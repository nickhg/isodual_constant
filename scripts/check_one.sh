#!/bin/bash
# Re-check one file of the published dataset with the checkers from scripts/build.sh, as the final re-check did:
#   1. its SHA-256 must equal the one recorded at production (manifests/*.sha256);
#   2. the checker must report FAILURES 0 (check_phi is retried at EXTRA = 24, 28, 36, 48 as in replay2.sh: EXTRA only
#      sets how hard check_phi tries before giving up on a leaf, never what it accepts);
#   3. for a stage-1 log, the collected boxes it outputs must be byte-identical to the stored <tag>.col.xz next to it.
# FILE is a path into a copy of the dataset that keeps its layout, e.g. DATA/eta5_stage2/logs/<slice>/<job>.logs.tar.xz,
# DATA/eta5_stage1/<slice>/<tag>.log.gz, DATA/phi5_le_P/<tag>.log.gz or DATA/phi5_lt7/slice_<k>_480.log; the part of
# the path from the top-level directory on is the manifest entry.
# Prints the checker's summary line and "OK" or "FAILED: reason"; the full output goes to OUT (default ./check_out/).
# Usage: scripts/check_one.sh FILE...
set -u
R=$(cd "$(dirname "$0")/.." && pwd); B=$R/build; OUT=${OUT:-check_out}; mkdir -p "$OUT"
[ -x $B/check_phi ] && [ -x $B/check_lift ] || { echo "run scripts/build.sh first"; exit 1; }
st=0
for f in "$@"; do
    base=$(basename "$f")
    rel=$(printf '%s\n' "$f" | grep -oE '(phi5_lt7|phi5_le_P|eta5_stage1|eta5_stage2)/.*$')
    want=$(awk -v r="$rel" '$2 == r {print $1}' $R/manifests/*.sha256 | head -1)
    have=$(sha256sum < "$f" | cut -d' ' -f1)
    o=$OUT/${base%%.*}.out; T=$(mktemp -d)
    res=OK
    if [ -z "$want" ]; then res="FAILED: $rel not in the manifests"
    elif [ "$have" != "$want" ]; then res="FAILED: SHA-256 differs from the manifest"
    else
        cphi() { local x; for x in 24 28 36 48; do env "$@" EXTRA=$x $B/check_phi $T/log > "$o" 2>&1; grep -q "FAILURES 0" "$o" && return 0; done; return 1; }
        case $rel in
        eta5_stage2/logs/*|eta5_stage2/fix/*)
            j=${base%.logs.tar.xz}; j=${j%.fix}; [ "${rel#eta5_stage2/fix/}" != "$rel" ] && o=$OUT/$j.fix.out   # cover_s2.py's names
            mkdir -p $T/logs/$j; o=$(realpath "$o")
            tar -xJf "$f" -C $T/logs/$j && (cd $T && $B/check_lift $(ls logs/$j/*.log | sort -V)) > "$o" 2>&1
            grep -q "FAILURES 0$" "$o" || res="FAILED: check_lift" ;;
        eta5_stage1/*.log.gz)
            col=${f%.log.gz}.col.xz; wcol=$(awk -v r="${rel%.log.gz}.col.xz" '$2 == r {print $1}' $R/manifests/eta5_stage1_logs.sha256)
            if [ ! -f "$col" ] || [ "$(sha256sum < "$col" | cut -d' ' -f1)" != "$wcol" ]; then res="FAILED: missing or altered $(basename "$col")"
            else
                gzip -dc "$f" > $T/log; cphi PVD5=$R/inputs/pivd5_U.txt STAGE2=$T/col || res="FAILED: check_phi"
                [ "$res" = OK ] && { xz -dc "$col" | cmp -s - $T/col || res="FAILED: collected boxes differ from $(basename "$col")"; }
                [ "$res" = OK ] && echo "collected boxes identical to the stored ones" >> "$o"
            fi ;;
        phi5_le_P/*) gzip -dc "$f" > $T/log; cphi MSTAR=$R/inputs/mstar_U.txt || res="FAILED: check_phi" ;;
        phi5_lt7/*) cp "$f" $T/log; cphi || res="FAILED: check_phi" ;;
        *) res="FAILED: unknown file type" ;;
        esac
    fi
    rm -rf $T
    echo "$rel: $(grep -h 'FAILURES' "$o" 2>/dev/null | tail -n 1) $res"
    [ "$res" = OK ] || st=1
done
exit $st
