#!/bin/bash
# Instance side of stage 1 of the eta_5 = 6 run (a copy of run_phi65.sh): Phi_5 < 6 outside the pi_v D5 balls and the
# collected set C.  Independent jobs, each certified, checked and packed for the local poller (pull_phi65.sh with
# PACK="log.gz txt out col.xz"), which keeps the authoritative ledger on the Chromebook.
# A job is only "done" once the poller has verified its packed log locally; done/ here is just the outbox.
# Per job: certify (PHI = 6, c = 1.68749999 < 27/16, pi_v D5 balls of radius 1/20, COLLECT = W = 0.12 with HOPELESS = 0)
#   -> check_phi (PVD5, EXTRA=24, retried at 28, 36 and 48; if still failing, certify again with NOCOMB=1 and re-check),
#   whose STAGE2 output (the C leaves of the checked log, exactly) is the stage-2 input
#   -> gzip the log, xz the C boxes, sha256 -> done/<tag>.{log.gz,col.xz,sha256,txt,out}.  Resumable: a job with done/<tag>.sha256
#   is skipped.  A job whose check still fails at EXTRA=28 goes to failed/<tag>.out and is not packed.
# Usage: run_phi65.sh JOBS [P]      (JOBS: one k/K or k/K/j/J per line, hardest first; a sub-job of a split piece,
#   k/K/j/J/s<path>, is followed on its line by its exact root box, see split_phi65.py)
set -u
JOBS=$1; P=${2:-$(nproc)}
PHI=6; C=1.68749999; W=0.12     # c < 6 / (32/9) = 27/16, so the root box is valid (kappa^2 <= 32/9)
mkdir -p results certlogs done failed
# VERSION: the git commit these sources were deployed from (launch_phi65.sh), stamped into every log and check
[ -s VERSION ] || { echo "missing VERSION"; exit 1; }
REV=$(cat VERSION)
gcc -O2 -frounding-math -DGITREV="\"$REV\"" -o certify certify.c -lglpk -lm || exit 1
gcc -O2 -frounding-math -DGITREV="\"$REV\"" -o check_phi check_phi.c -lm || exit 1
job() {
    # a line is "spec" or "spec box": a sub-job of a split piece carries its exact root box (18 hex numbers).
    # xargs -I substitutes the line unquoted, so it arrives as separate words: take all of them
    local line="$*"; local spec=${line%% *} box=""; [ "$spec" != "$line" ] && box=${line#* }
    set -- "$spec"
    # only set for sub-jobs: certify treats any ROOTBOX, even an empty one, as the root
    if [ -n "$box" ]; then export ROOTBOX="$box"; else unset ROOTBOX; fi
    t=$(echo "$1" | tr / _)
    [ -f done/$t.sha256 ] && return
    grep -q "CERTIFIED\|COMPLETE" results/$t.txt 2>/dev/null ||
        INVBALL=1 NOELEMCG=1 HOPELESS=0 COLLECT=$W COLLECT_FILE=results/$t.certcol CERTLOG=certlogs PHI=$PHI PROGRESS=100000 \
            ./certify 4 $C 100000000000 pivd5_balls.txt $1 > results/$t.txt 2> results/$t.err
    grep -q "CERTIFIED\|COMPLETE" results/$t.txt || { echo "$1 not certified" >> failed/certify.txt; return; }
    ok=0
    for x in 24 28 36 48; do
        PVD5=pivd5_U.txt STAGE2=results/$t.col EXTRA=$x ./check_phi certlogs/slice_$t.log > results/$t.out 2>&1 && grep -q "FAILURES 0" results/$t.out && { ok=1; break; }
    done
    # still failing: a leaf closed by the LP combination of several flags, which the single-flag checker cannot follow
    # near M* (465/480/4/8/s100).  Re-certify without it (NOCOMB: same certify, every leaf then closes with one
    # flag, which is its hint) and check again.
    if [ $ok = 0 ]; then
        NOCOMB=1 INVBALL=1 NOELEMCG=1 HOPELESS=0 COLLECT=$W COLLECT_FILE=results/$t.certcol CERTLOG=certlogs PHI=$PHI PROGRESS=100000 \
            ./certify 4 $C 100000000000 pivd5_balls.txt $1 > results/$t.txt 2> results/$t.err
        grep -q "CERTIFIED\|COMPLETE" results/$t.txt &&
            for x in 24 28 36 48; do
                PVD5=pivd5_U.txt STAGE2=results/$t.col EXTRA=$x ./check_phi certlogs/slice_$t.log > results/$t.out 2>&1 && grep -q "FAILURES 0" results/$t.out && { ok=1; break; }
            done
    fi
    [ $ok = 1 ] || { cp results/$t.out failed/$t.out; return; }
    gzip -9 -c certlogs/slice_$t.log > done/$t.log.gz.tmp && gzip -t done/$t.log.gz.tmp && mv done/$t.log.gz.tmp done/$t.log.gz &&
        xz -T1 -c results/$t.col > done/$t.col.xz.tmp && xz -t done/$t.col.xz.tmp && mv done/$t.col.xz.tmp done/$t.col.xz &&
        cp results/$t.txt done/$t.txt && cp results/$t.out done/$t.out &&
        (cd done && sha256sum $t.log.gz $t.col.xz $t.txt $t.out > $t.sha256.tmp && mv $t.sha256.tmp $t.sha256) &&
        rm -f certlogs/slice_$t.log results/$t.col results/$t.certcol
}
export -f job; export PHI C W
xargs -P "$P" -I{} bash -c 'job {}' < "$JOBS"
echo "finished $(date -u +%FT%TZ) packed $(ls done/*.sha256 2>/dev/null | wc -l) failed $(ls failed/*.out 2>/dev/null | wc -l)" > done/FINISHED
