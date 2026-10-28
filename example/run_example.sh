#!/bin/bash
# The worked example of VERIFY.md (tier 2): every kind of check in the final re-check, on a few real files of the
# published certificate of eta_5 = 6, plus the checker's own tests.  About 5 minutes on one core.
# Data (example/data/, 3.9 MB; each file also listed with its SHA-256 in manifests/):
#   eta5_stage1/391/391_480_4_8_s01110.{log.gz,col.xz}       a stage-1 log and its collected boxes
#   eta5_stage1/477/477_480_4_8_s*.col.xz                     the collected boxes of the 32 stage-1 jobs of piece 477_480_4_8
#   eta5_stage1_merge/477_480_4_8.{col,trees}             that piece's merged stage-2 roots and merge trees (in the dataset: eta5_stage1_merge.tar)
#   eta5_stage2/chunks/477_480_4_8_c000.col, jobs.txt     its one stage-2 job
#   eta5_stage2/logs/409/409_480_1_8_c023.logs.tar.xz         a stage-2 archive (97 logs, 114,074 leaves, one of them (B))
#   reference_outputs/                                    the final re-check's check_lift outputs for both stage-2 jobs
# The stage-2 archive of 477_480_4_8_c000 is 67 MB and not included: step 4 uses the final re-check's output for it,
# unless ARCHIVE_477=<path to 477_480_4_8_c000.logs.tar.xz> is set, in which case it is re-checked here first.
# Usage: example/run_example.sh       (needs gcc, python3; python-flint for step 5 and mpmath for step 7, else skipped)
set -u
cd "$(dirname "$0")/.."; D=example/data; W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; export OUT=$W/out
PV=9d7c51a8414a68bcb445c8b8debacfeed4f8cf15,59fdb6acb6e41b8d621173c92b6c6ab89a91c2e2   # stage-2 producer versions
bad=0; step() { echo; echo "== $*"; }; chk() { if "$@"; then echo "   -> passed"; else echo "   -> FAILED"; bad=$((bad + 1)); fi; }

step "1. build the checkers (as in the final re-check)"
chk scripts/build.sh

step "2. stage 1: check_phi re-proves a log; its collected boxes must equal the stored ones (about 1 minute)"
chk scripts/check_one.sh $D/eta5_stage1/391/391_480_4_8_s01110.log.gz

step "3. stage 2: check_lift re-proves an archive; its output must equal the final re-check's, byte for byte"
chk scripts/check_one.sh $D/eta5_stage2/logs/409/409_480_1_8_c023.logs.tar.xz
chk cmp $OUT/409_480_1_8_c023.out $D/reference_outputs/409_480_1_8_c023.out

step "4. coverage of piece 477_480_4_8: merged roots = exact unions of all stage-1 C leaves; every root proved"
hashes() { (cd $D && grep -E ' eta5_stage1/477/477_480_4_8_s[01]+\.col\.xz$' ../../manifests/eta5_stage1_logs.sha256 | sha256sum --quiet -c -); }
chk hashes
chk python3 scripts/merge_check.py $D/eta5_stage1_merge/477_480_4_8.col $D/eta5_stage1_merge/477_480_4_8.trees $(ls $D/eta5_stage1/477/477_480_4_8_s*.col.xz)
mkdir -p $W/s2/chunks $W/ver; ln -s "$(realpath $D/eta5_stage1_merge)" $W/s2/merged
cp $D/eta5_stage2/chunks/477_480_4_8_c000.col $W/s2/chunks/; grep '^477_480_4_8_' $D/eta5_stage2/jobs.txt > $W/s2/jobs.txt
if [ -n "${ARCHIVE_477:-}" ]; then
    chk scripts/check_one.sh "$ARCHIVE_477"; cp $OUT/477_480_4_8_c000.out $W/ver/
else
    echo "   (using the final re-check's output for 477_480_4_8_c000; set ARCHIVE_477 to re-check its 67 MB archive)"
    cp $D/reference_outputs/477_480_4_8_c000.out $W/ver/
fi
chk env CHECKVERS=b5e57c27466195192cfc8a46dc8bcd87a004de77 python3 scripts/cover_s2.py $W/s2 $W/ver $PV 477_480_4_8

mkdir -p $W/a409; tar -xJf $D/eta5_stage2/logs/409/409_480_1_8_c023.logs.tar.xz -C $W/a409   # 67.log has N, F, P and B leaves
step "5. the independent checker (Python + Arb): the (B) leaf, and a 0.2% sample of the other leaves of 67.log"
if python3 -c "import flint" 2>/dev/null; then
    chk python3 independent/check_lift2.py --types B $W/a409/67.log
    chk python3 independent/check_lift2.py --types NFP --sample 0.002 $W/a409/67.log
else echo "   (skipped: pip install python-flint)"; fi

step "6. 23 deliberately broken certificates: check_lift must reject each one"
chk scripts/check_lift_tests.sh build/check_lift $W/a409/67.log

step "7. CORE-MATH exp and log against mpmath at 60 digits, rounding up and to nearest (20,000 points)"
if python3 -c "import mpmath" 2>/dev/null; then chk scripts/test_crexplog.sh 20000; else echo "   (skipped: pip install mpmath)"; fi

step "8. 12 negative controls at M* and pi_v D5: without the balls, check_phi must reject a box containing the maximizer"
chk scripts/check_phi_controls.sh build/check_phi

step "9. 30 mutation tests of the coverage scripts: each damaged input must be rejected"
chk scripts/coverage_mutation_tests.sh

echo; [ $bad = 0 ] && echo "EXAMPLE PASSED" || { echo "EXAMPLE FAILED: $bad checks"; exit 1; }
