#!/bin/bash
# Tier 4, coverage: that the leaves of each certificate tile its search space.  Run after check_one.sh has re-checked
# every log, with OUT = the directory of its outputs (for stage 2).  The outputs are compared with the ones recorded
# in manifests/coverage/.
#   1. Phi_5 < 7: the 575 roots tile one box, and it contains the box B of the appendix (root_box_check.py);
#   2. Phi_5 <= P and 3. eta_5 = 6 stage 1: the roots, with the split manifests, tile the reference root box;
#   4. stage 2: every merged root is an exact union of stage-1 collected boxes (merge_check.py), and every merged root
#      of all 742 pieces is proved by an accepted stage-2 log (cover_s2.py, on the outputs in OUT).
# Usage: OUT=check_out scripts/coverage_all.sh DATA
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd); DATA=$1; OUT=${OUT:-check_out}; M=$R/manifests/coverage; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
heads() {    # heads DIR PATTERN: the first 10 lines (the ROOT header and box) of each log, uncompressed, as
             # slice_<tag>.log (the names cover_check.py expects, also for the children of split pieces)
    mkdir -p $T/$1; find "$DATA/$1" -name "$2" | while read -r f; do b=$(basename "$f"); b=${b%.gz}; b=${b%.log}; b=slice_${b#slice_}
        if [ "${f%.gz}" != "$f" ]; then { gzip -dc "$f" || true; } | head -n 10 > $T/$1/$b.log; else head -n 10 "$f" > $T/$1/$b.log; fi; done; }
echo "== 1. Phi_5 < 7"
heads phi5_lt7 'slice_*.log'        # includes the 96 pieces slice_469_480_<i>_96.log
python3 $R/scripts/cover_check.py $T/phi5_lt7/*.log | tee $T/p7.out | python3 $R/scripts/root_box_check.py
cmp <(head -n 10 $T/p7.out) <(head -n 10 $M/phi5_lt7/cover_check.out) && echo "matches the recorded output"
for c in phi5_le_P eta5_stage1; do
    echo "== $c"
    heads $c '*.log.gz'
    python3 $R/scripts/cover_check.py --ref $M/$c/ref_root.log $(printf -- '--manifest %s ' $M/$c/manifests/*.json) $T/$c/*.log | tee $T/$c.out | head -n 1
    ref=$M/$c/cover_check.out; [ -f $ref ] || ref=$M/$c/cover_check.txt
    cmp <(head -n 10 $T/$c.out) <(head -n 10 $ref) && echo "matches the recorded output"
done
echo "== stage 2"
mkdir -p $T/s2; tar -xf $DATA/eta5_stage1_merge.tar -C $T/s2; tar -xf $DATA/eta5_stage2/chunks.tar -C $T/s2; cp $DATA/eta5_stage2/jobs.txt $T/s2/
n=0; while read -r p; do      # merge_check per piece, on the collected boxes of its 32 stage-1 jobs (in the order prep.sh used)
    python3 $R/scripts/merge_check.py $T/s2/merged/$p.col $T/s2/merged/$p.trees $(ls $DATA/eta5_stage1/${p%%_*}/${p}_s*.col.xz) > /dev/null; n=$((n + 1))
done < $M/eta5_stage2/pieces.txt; echo "merge check OK on $n pieces"
CHECKVERS=b5e57c27466195192cfc8a46dc8bcd87a004de77 python3 $R/scripts/cover_s2.py $T/s2 $OUT \
    9d7c51a8414a68bcb445c8b8debacfeed4f8cf15,59fdb6acb6e41b8d621173c92b6c6ab89a91c2e2 $(cat $M/eta5_stage2/pieces.txt) | tail -n 1
