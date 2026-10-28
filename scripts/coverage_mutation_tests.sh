#!/bin/bash
# Mutation tests for the coverage layer (cover_check.py, root_box_check.py, merge_check.py, cover_s2.py).
# Each script is first run on correct inputs (it must accept), then on deliberately damaged copies (it must reject):
# a deleted, duplicated or truncated file or line, an endpoint shifted by one ulp, a wrong version or header, and an
# output taken from another job.  Inputs: the worked example's data (piece 477_480_4_8) and small synthetic root
# logs built from the reference root box of Phi_5 <= P.  About 30 seconds; needs python3.
# Usage: scripts/coverage_mutation_tests.sh
set -u
cd "$(dirname "$0")/.."; R=$PWD; D=example/data; W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
pass=0; fail=0
ok()  { if "$@" > $W/last.txt 2>&1; then pass=$((pass + 1)); else fail=$((fail + 1)); echo "  MUST ACCEPT, rejected: $*"; tail -n 2 $W/last.txt; fi; }
bad() { local what=$1; shift; if "$@" > $W/last.txt 2>&1; then fail=$((fail + 1)); echo "  MUST REJECT, accepted: $what"; else pass=$((pass + 1)); echo "  rejected (as required): $what"; fi; }
ulp() { python3 -c "import sys,math; x=float.fromhex(sys.argv[1]); print(math.nextafter(x, math.inf).hex())" "$1"; }

echo "== cover_check.py: synthetic root logs tiling the reference root box of Phi_5 <= P"
REF=manifests/coverage/phi5_le_P/ref_root.log
mk() {   # mk DIR NAME LO HI [HEADER]: a root log whose box is REF's with coordinate 1 replaced by [LO, HI]
    { if [ -n "${5:-}" ]; then echo "$5"; else head -n 1 $REF; fi
      echo "$3 $4"; sed -n '3,10p' $REF; echo "END incomplete"; } > $1/$2; }
lo=$(sed -n 2p $REF | cut -d' ' -f1); hi=$(sed -n 2p $REF | cut -d' ' -f2)
m1=$(python3 -c "print((float.fromhex('$lo')/2 + float.fromhex('$hi')/2).hex())"); m2=$(python3 -c "print((float.fromhex('$m1')/2 + float.fromhex('$hi')/2).hex())")
good() { rm -rf $W/cc; mkdir -p $W/cc; mk $W/cc slice_a.log $lo $m1; mk $W/cc slice_b.log $m1 $m2; mk $W/cc slice_c.log $m2 $hi; }
good; ok python3 scripts/cover_check.py --ref $REF $W/cc/slice_a.log $W/cc/slice_b.log $W/cc/slice_c.log
good; rm $W/cc/slice_b.log;                       bad "a missing log (gap)"                 python3 scripts/cover_check.py --ref $REF $W/cc/*.log
good; mk $W/cc slice_b.log $(ulp $m1) $m2;        bad "an endpoint shifted by one ulp"      python3 scripts/cover_check.py --ref $REF $W/cc/*.log
good; cp $W/cc/slice_b.log $W/cc/slice_b2.log;    bad "a duplicated log (overlap)"          python3 scripts/cover_check.py --ref $REF $W/cc/*.log
good; mk $W/cc slice_c.log $m2 $hi "$(head -n 1 $REF | sed 's/c=1.8204/c=1.8205/')"
                                                  bad "a log with a different header (c)"   python3 scripts/cover_check.py --ref $REF $W/cc/*.log
good; sed -i '3s/^\(\S*\) \(\S*\)$/\1 0x1.27p-1/' $W/cc/slice_a.log
                                                  bad "coordinate 2 differs in one log"     python3 scripts/cover_check.py --ref $REF $W/cc/*.log
good; mk $W/cc slice_c.log $m2 $(python3 -c "import math; print(math.nextafter(float.fromhex('$hi'), -math.inf).hex())")
                                                  bad "union one ulp short of the reference" python3 scripts/cover_check.py --ref $REF $W/cc/*.log
good; mk $W/cc slice_b_s0.log $m1 $m2;            bad "a sub-log without its manifest"      python3 scripts/cover_check.py --ref $REF $W/cc/*.log

echo "== root_box_check.py: the recorded union of each certificate must contain B, exactly"
for x in "phi5_lt7/cover_check.out 63/32" "phi5_le_P/cover_check.out 2048/1125" "eta5_stage1/cover_check.txt 27/16"; do
    set -- $x; ok python3 scripts/root_box_check.py $2 < manifests/coverage/$1
done
bad "Phi_5 < 7 with c below the program's 1.96874999"   python3 scripts/root_box_check.py 196874998/100000000 < manifests/coverage/phi5_lt7/cover_check.out
bad "stage 1 with c below the program's 1.68749999"     python3 scripts/root_box_check.py 168749998/100000000 < manifests/coverage/eta5_stage1/cover_check.txt
python3 - "$W/short.out" <<'EOF'
import sys, re
s = open("manifests/coverage/phi5_lt7/cover_check.out").read()
m = re.search(r"coordinate 2: \[(\S+), (\S+)\]", s); hi = float.fromhex(m.group(2))
import math; s = s.replace(m.group(2), math.nextafter(hi, -math.inf).hex(), 1)   # y2 upper end one ulp below ybar_2
open(sys.argv[1], "w").write(s)
EOF
bad "a union whose y2 range stops one ulp short"         python3 scripts/root_box_check.py 63/32 < $W/short.out

echo "== merge_check.py: piece 477_480_4_8 (merged roots = exact unions of the stage-1 collected boxes)"
MC="$D/eta5_stage1_merge/477_480_4_8.col"; MT="$D/eta5_stage1_merge/477_480_4_8.trees"; COLS=$(ls $D/eta5_stage1/477/477_480_4_8_s*.col.xz)
ok python3 scripts/merge_check.py $MC $MT $COLS
mkdir -p $W/mc
python3 - "$MC" "$W/mc/shift.col" <<'EOF'
import sys, math
L = open(sys.argv[1]).read().split("\n"); v = L[0].split("#")[0].split()
v[3] = math.nextafter(float.fromhex(v[3]), math.inf).hex()          # first merged box, coordinate 2 upper end + 1 ulp
L[0] = " ".join(v); open(sys.argv[2], "w").write("\n".join(L))
EOF
bad "a merged box shifted by one ulp"                   python3 scripts/merge_check.py $W/mc/shift.col $MT $COLS
bad "one stage-1 collected-box file left out"            python3 scripts/merge_check.py $MC $MT $(echo $COLS | cut -d' ' -f2-)
NE=$(for f in $COLS; do [ -n "$(xz -dc $f | head -c 1)" ] && { echo $f; break; }; done)   # first non-empty file
with() { for f in $COLS; do [ $f = $NE ] && echo $1 || echo $f; done; }                   # COLS with NE replaced
xz -dc $NE | sed '$d' > $W/mc/del.col
bad "one collected box deleted"                          python3 scripts/merge_check.py $MC $MT $(with $W/mc/del.col)
{ xz -dc $NE; xz -dc $NE | head -n 1; } > $W/mc/dup.col
bad "one collected box duplicated (left unused)"         python3 scripts/merge_check.py $MC $MT $(with $W/mc/dup.col)
sed '$d' $MT > $W/mc/short.trees
bad "the last merge tree deleted"                        python3 scripts/merge_check.py $MC $W/mc/short.trees $COLS

echo "== cover_s2.py: piece 477_480_4_8 (every merged root proved by an accepted stage-2 log)"
PV=9d7c51a8414a68bcb445c8b8debacfeed4f8cf15,59fdb6acb6e41b8d621173c92b6c6ab89a91c2e2; CV=b5e57c27466195192cfc8a46dc8bcd87a004de77
setup() { rm -rf $W/s2 $W/ver; mkdir -p $W/s2/chunks $W/s2/merged $W/ver
    cp $MC $W/s2/merged/; cp $D/eta5_stage2/chunks/477_480_4_8_c000.col $W/s2/chunks/
    grep '^477_480_4_8_' $D/eta5_stage2/jobs.txt > $W/s2/jobs.txt; cp $D/reference_outputs/477_480_4_8_c000.out $W/ver/; }
s2() { CHECKVERS=$CV python3 scripts/cover_s2.py $W/s2 $W/ver $PV 477_480_4_8; }
O=$W/ver/477_480_4_8_c000.out
setup; ok s2
setup; n=$(grep -n ' OK$' $O | head -n 1 | cut -d: -f1); sed -i "${n}d" $O
                                                  bad "one accepted-log line deleted"           s2
setup; sed -i 's/FAILURES 0$/FAILURES 1/' $O;     bad "FAILURES 0 changed to FAILURES 1"        s2
setup; grep -v 'FAILURES' $O > $O.t && mv $O.t $O
                                                  bad "output truncated (no FAILURES line)"     s2
setup; n=$(grep -n ' OK$' $O | head -n 1 | cut -d: -f1); sed -n "${n}p" $O >> $O
                                                  bad "one accepted-log line duplicated"        s2
setup; n=$(grep -n ' OK$' $O | head -n 1 | cut -d: -f1); sed -i "${n}s/check_lift ver=$CV/check_lift ver=0000000/" $O
                                                  bad "a log checked by an unknown checker version" s2
setup; sed -i '$d' $W/s2/chunks/477_480_4_8_c000.col
                                                  bad "a chunk with its last box deleted"       s2
setup; python3 - $W/s2/chunks/477_480_4_8_c000.col <<'EOF'
import sys, math
L = open(sys.argv[1]).read().split("\n"); v = L[0].split("#")[0].split()
v[0] = math.nextafter(float.fromhex(v[0]), math.inf).hex(); L[0] = " ".join(v); open(sys.argv[1], "w").write("\n".join(L))
EOF
                                                  bad "a chunk box shifted by one ulp"          s2
setup; sed "s#logs/409_480_1_8_c023/#logs/477_480_4_8_c000/#" $D/reference_outputs/409_480_1_8_c023.out > $O
                                                  bad "the output of another job, relabelled"   s2
setup; rm $O;                                     bad "the job's output missing"                s2

echo; echo "$pass passed, $fail failed"
[ $fail = 0 ] && echo "COVERAGE MUTATION TESTS PASSED" || { echo "COVERAGE MUTATION TESTS FAILED"; exit 1; }
