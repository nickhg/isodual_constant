#!/bin/bash
# Instance side of stage 2 of the eta_5 = 6 run: for each job (a chunk of merged collected boxes), the producer stage2
# writes one certificate log per box, check_lift (an independent program) checks them all, and the job is packed for
# the local poller (pull_phi65.sh with PACK="logs.tar.xz txt out").
# A box whose log check_lift rejects is re-certified alone with a larger producer margin (1e-4, then 1e-3) and checked
# again (not a box that ran out of time: a margin cannot help it).  A box still rejected, or out of time (TCAP s, default
# 1800), is bisected at the midpoint of its widest coordinate and the halves proved in its place (logs <i>.0, <i>.1,
# recursively to depth 6; cover_s2.py recomputes the same bisection).  A box that fails even then sends the job to
# failed/<job>.out (not packed).  Resumable: jobs with done/<job>.sha256 are skipped.
# Usage: run_s2.sh JOBS [P]      (JOBS: lines "job_name chunk_file"; the chunk files are in chunks/)
set -u
JOBS=$1; P=${2:-$(nproc)}
mkdir -p results logs done failed
[ -s VERSION ] || { echo "missing VERSION"; exit 1; }
REV=$(cat VERSION)
gcc -O2 -DGITREV="\"$REV\"" -o stage2 stage2.c -lglpk -lm || exit 1
gcc -O2 -frounding-math -DGITREV="\"$REV\"" -o check_lift check_lift.c -lm || exit 1
TCAP=${TCAP:-1800}
# the bisection (cover_s2.py uses the same rule): widest coordinate, ties to the lowest index, at the float midpoint
BISECT='import sys
v = [float.fromhex(x) for x in sys.argv[1].split()[:18]]; lo, hi = v[0::2], v[1::2]
t = max(range(9), key=lambda i: (hi[i] - lo[i], -i)); m = (lo[t] + hi[t]) / 2
if not lo[t] < m < hi[t]: sys.exit(1)
for a, b in ((lo[t], m), (m, hi[t])):
    l, h = list(lo), list(hi); l[t], h[t] = a, b
    print(" ".join("%s %s" % (x.hex(), y.hex()) for x, y in zip(l, h)))'
checkall() { ./check_lift $(ls logs/$1/*.log | sort -V) > results/$1.out 2>&1; }
badboxes() { grep ' FAILED$\|incomplete\|not complete' results/$1.out | sed 's|^logs/[^/]*/\([0-9.]*\)\.log.*|\1|' | sort -uV; }
timedout() { awk -v i="$2" '$1 == i {r = $4} END {exit r != "TIME"}' results/$1.txt; }
prove_split() {     # prove_split JOB ID BOX DEPTH: prove the two halves of BOX as logs ID.0 and ID.1
    local name=$1 id=$2 d=$4 t h ok=0
    t=$(mktemp -d)
    python3 -c "$BISECT" "$3" > $t/halves || { rm -rf $t; return 1; }
    echo "split box $id (depth $d)" >> results/$name.err
    for h in 0 1; do
        sed -n "$((h + 1))p" $t/halves > $t/in.col; rm -rf $t/l; mkdir -p $t/l
        LOGDIR=$t/l TCAP=$TCAP BOXCAP=100000000 ./stage2 $t/in.col 2>> results/$name.err | sed -n "s/^0 /$id.$h /p" >> results/$name.txt
        mv $t/l/0.log logs/$name/$id.$h.log 2>/dev/null
        if ! ./check_lift logs/$name/$id.$h.log 2>&1 | grep -q "FAILURES 0"; then
            rm -f logs/$name/$id.$h.log
            [ $d -lt 6 ] && prove_split $name $id.$h "$(cat $t/in.col)" $((d + 1)) || { ok=1; break; }
        fi
    done
    rm -rf $t; return $ok
}
job() {
    local name=$1 chunk=chunks/$2
    [ -f done/$name.sha256 ] && return
    rm -rf logs/$name; mkdir -p logs/$name
    LOGDIR=logs/$name TCAP=$TCAP BOXCAP=100000000 ./stage2 $chunk > results/$name.txt 2> results/$name.err
    local i m
    checkall $name
    for m in 1e-4 1e-3; do
        grep -q "FAILURES 0" results/$name.out && break
        for i in $(badboxes $name); do
            timedout $name $i && continue
            echo "rerun box $i with MARGIN=$m" >> results/$name.err
            MARGIN=$m LOGDIR=logs/$name TCAP=$TCAP BOXCAP=100000000 ./stage2 $chunk $i 1 >> results/$name.txt 2>> results/$name.err
        done
        checkall $name
    done
    if ! grep -q "FAILURES 0" results/$name.out; then
        local unproved=""
        for i in $(badboxes $name); do
            rm -f logs/$name/$i.log
            prove_split $name $i "$(sed -n "$((i + 1))p" $chunk)" 1 || { unproved=$i; break; }
        done
        checkall $name
        # a bisection that fails at depth 6 leaves no log for the failed half, so check_lift (which sees only the
        # logs present) cannot notice the gap: fail the job here (cover_s2.py also catches it)
        [ -n "$unproved" ] && { cp results/$name.out failed/$name.out; echo "box $unproved: bisection failed at depth 6" >> failed/$name.out; return; }
    fi
    grep -q "FAILURES 0" results/$name.out || { cp results/$name.out failed/$name.out; return; }
    tar -cJf done/$name.logs.tar.xz.tmp -C logs/$name . && xz -t done/$name.logs.tar.xz.tmp && mv done/$name.logs.tar.xz.tmp done/$name.logs.tar.xz &&
        cp results/$name.txt done/$name.txt && cp results/$name.out done/$name.out &&
        (cd done && sha256sum $name.logs.tar.xz $name.txt $name.out > $name.sha256.tmp) &&
        s3put $name && mv done/$name.sha256.tmp done/$name.sha256 && rm -rf logs/$name &&
        { [ -z "${S3:-}" ] || rm -f done/$name.logs.tar.xz; }
}
# with S3 (file S3 holding "s3://bucket/prefix"): upload the logs with S3's SHA-256 checksum before the job counts as
# packed (.sha256 is written last), then drop the local copy; the poller checks S3's checksum against .sha256
[ -s S3 ] && S3=$(cat S3)
s3put() {
    [ -z "${S3:-}" ] && return 0
    local b k try
    b=${S3#s3://}; b=${b%%/*}; k=${S3#s3://*/}
    for try in 1 2 3 4 5; do
        aws s3api put-object --bucket "$b" --key "$k/$1.logs.tar.xz" --body "done/$1.logs.tar.xz" \
            --checksum-algorithm SHA256 > /dev/null 2>> results/$1.err && return 0
        sleep $((try * 30))
    done
    echo "S3 upload of $1 failed" >> results/$1.err; echo "S3 upload failed (FAILURES 1)" > failed/$1.out; return 1
}
export S3 TCAP BISECT; export -f s3put checkall badboxes timedout prove_split job
xargs -P "$P" -L 1 bash -c 'job "$@"' _ < "$JOBS"
echo "finished $(date -u +%FT%TZ) packed $(ls done/*.sha256 2>/dev/null | wc -l) failed $(ls failed/*.out 2>/dev/null | wc -l)" > done/FINISHED
