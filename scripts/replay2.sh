#!/bin/bash
# Replay 2: re-check EVERY certificate log of the paper with one static build of check_phi and check_lift whose exp and
# log are CORE-MATH's correctly rounded cr_exp/cr_log (third_party/core-math), widened by one ulp.
# Job names (one per line of JOBS):
#   p7_<tag>    Phi_5 < 7 (App. phi7), runs/cloud_phi7/certlogs/slice_<tag>.log          EXTRA 24 (28, 36, 48 if needed)
#   p65_<tag>   Phi_5 <= P (App. phi65), the 4,757 logs of final/selected.txt, MSTAR       EXTRA 24, 28, 36, 48
#   s1_<tag>    eta_5 = 6 stage 1 (App. s2), PVD5 + STAGE2; the collected boxes must be byte-identical to the stored ones
#   s2_<job>    stage 2, S3 eta6_s2/full/<job>.logs.tar.xz, check_lift on all logs in run_s2.sh's order
#   s2f_<job>   stage-2 fix-up, S3 eta6_s2/fix/<job>.fix.logs.tar.xz
# For each job the SHA-256 of what was fetched is compared with the one recorded when the log was produced (s2: recorded
# in done/<name>.arch and compared on the Chromebook against verified/<job>.sha256, as in replay_s2.sh).  Output:
# done/<name>.out (checker output), done/<name>.arch (hash record), done/<name>.sha256 (over both); a job that fails
# any test goes to failed/<name>.out.  done/BINARY: the sha256 of both binaries and the gcc version.  Resumable.
# Usage: replay2.sh JOBS [P]       (file S3: s3://bucket, holding replay2/ and eta6_s2/)
set -u
JOBS=$1; P=${2:-$(nproc)}
[ -s VERSION ] && [ -s S3 ] || { echo "missing VERSION or S3"; exit 1; }
REV=$(cat VERSION); S3=$(cat S3); B=${S3#s3://}; B=${B%%/*}
mkdir -p done failed tmp logs src
gcc -O2 -mfma -ffp-contract=off -c exp.c -o cr_exp.o && gcc -O2 -mfma -ffp-contract=off -I. -c log.c -o cr_log.o || exit 1
for c in check_phi check_lift; do
    gcc -O2 -static -frounding-math -ffp-contract=off -DGITREV="\"$REV\"" -o $c $c.c cr_exp.o cr_log.o -lm || exit 1
done
{ sha256sum check_phi check_lift; gcc --version | head -1; } > done/BINARY
get() { local try; for try in 1 2 3 4 5; do aws s3api get-object --bucket "$B" --key "$1" "$2" > /dev/null 2>&1 && return 0; rm -f "$2"; sleep $((try * 20)); done; return 1; }
# one-time inputs of the check_phi runs (skipped if already unpacked)
if [ ! -f src/READY ]; then
    for f in p7_certlogs.sha256 p65_selected.txt pivd5_U.txt mstar_U.txt; do get replay2/$f src/$f || exit 1; done
    get replay2/p7.tar src/p7.tar && tar -xf src/p7.tar -C src && rm src/p7.tar || exit 1
    get replay2/s1.tar src/s1.tar && mkdir -p src/s1 && tar -xf src/s1.tar -C src/s1 && rm src/s1.tar || exit 1
    get replay2/p65_incoming.tar src/p65.tar && mkdir -p src/p65 &&
        tar -xf src/p65.tar -C src/p65 $(awk '{print "incoming/" $2 "/" $1 ".log.gz"}' src/p65_selected.txt) && rm src/p65.tar || exit 1
    touch src/READY
fi
fail() { echo "$2 (FAILURES 1)" >> tmp/$1.out; mv tmp/$1.out failed/$1.out; rm -rf tmp/$1.*; }
cphi() {   # cphi NAME LOG ENV... : check_phi at EXTRA 24, 28, 36, 48 until FAILURES 0
    local n=$1 log=$2 x; shift 2
    for x in 24 28 36 48; do
        env "$@" EXTRA=$x ./check_phi "$log" > tmp/$n.out 2>&1
        grep -q "FAILURES 0" tmp/$n.out && { echo "EXTRA=$x" >> tmp/$n.out; return 0; }
    done
    return 1
}
job() {
    local n=$1 t h want a
    [ -f done/$n.sha256 ] && return
    case $n in
    p7_*) t=${n#p7_}; a=src/certlogs/slice_$t.log
        h=$(sha256sum < $a | cut -d' ' -f1); want=$(awk -v f=slice_$t.log '$2 == f {print $1}' src/p7_certlogs.sha256)
        echo "$h slice_$t.log recorded $want" > tmp/$n.arch
        [ -n "$want" ] && [ "$h" = "$want" ] || { fail $n "hash mismatch"; return; }
        cphi $n $a || { fail $n "check failed"; return; } ;;
    p65_*) t=${n#p65_}; read -r _ d want < <(awk -v t=$t '$1 == t' src/p65_selected.txt); a=src/p65/incoming/$d/$t.log.gz
        h=$(sha256sum < $a | cut -d' ' -f1); echo "$h $t.log.gz recorded $want" > tmp/$n.arch
        [ -n "$want" ] && [ "$h" = "$want" ] || { fail $n "hash mismatch"; return; }
        gzip -dc $a > tmp/$n.log && cphi $n tmp/$n.log MSTAR=src/mstar_U.txt || { fail $n "check failed"; return; } ;;
    s1_*) t=${n#s1_}; : > tmp/$n.arch
        for f in $t.log.gz $t.col.xz; do
            h=$(sha256sum < src/s1/$f | cut -d' ' -f1); want=$(awk -v f=$f '$2 == f {print $1}' src/s1/$t.sha256)
            echo "$h $f recorded $want" >> tmp/$n.arch
            [ -n "$want" ] && [ "$h" = "$want" ] || { fail $n "hash mismatch $f"; return; }
        done
        gzip -dc src/s1/$t.log.gz > tmp/$n.log && cphi $n tmp/$n.log PVD5=src/pivd5_U.txt STAGE2=tmp/$n.col ||
            { fail $n "check failed"; return; }
        xz -dc src/s1/$t.col.xz | cmp -s - tmp/$n.col || { fail $n "collected boxes differ"; return; }
        echo "collected boxes identical to the stored ones" >> tmp/$n.out ;;
    s2_*|s2f_*) if [ "${n#s2f_}" != "$n" ]; then t=${n#s2f_}; a=eta6_s2/fix/$t.fix.logs.tar.xz; else t=${n#s2_}; a=eta6_s2/full/$t.logs.tar.xz; fi
        get $a tmp/$n.tar.xz || { fail $n "download failed"; return; }
        echo "$(sha256sum < tmp/$n.tar.xz | cut -d' ' -f1) $(basename $a)" > tmp/$n.arch
        rm -rf logs/$t; mkdir -p logs/$t && tar -xJf tmp/$n.tar.xz -C logs/$t && rm -f tmp/$n.tar.xz || { fail $n "extract failed"; return; }
        ./check_lift $(ls logs/$t/*.log | sort -V) > tmp/$n.out 2>&1; rm -rf logs/$t
        grep -q "FAILURES 0" tmp/$n.out || { fail $n "check failed"; return; } ;;
    *) fail $n "unknown job type"; return ;;
    esac
    rm -f tmp/$n.log tmp/$n.col
    mv tmp/$n.arch done/$n.arch; mv tmp/$n.out done/$n.out
    (cd done && sha256sum $n.arch $n.out > $n.sha256.tmp && mv $n.sha256.tmp $n.sha256)
}
export B S3; export -f job cphi fail get
awk '{print $1}' "$JOBS" | xargs -P "$P" -n 1 bash -c 'job "$@"' _
echo "finished $(date -u +%FT%TZ) packed $(ls done/*.sha256 2>/dev/null | wc -l) failed $(ls failed/*.out 2>/dev/null | wc -l)" > done/FINISHED
