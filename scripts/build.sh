#!/bin/bash
# Build check_phi and check_lift exactly as in the paper's final re-check (replay 2): CORE-MATH's exp and log compiled
# with -mfma, the checkers with -frounding-math -ffp-contract=off, statically linked.  Output in build/.
# Usage: scripts/build.sh      (from the repository root; needs gcc and a static libc/libm)
set -eu
cd "$(dirname "$0")/.."; C=src/third_party/core-math; mkdir -p build
REV=b5e57c27466195192cfc8a46dc8bcd87a004de77      # the version stamp the paper's checker outputs carry
# -mfma only for CORE-MATH, which calls fma explicitly: the build used for every check in the paper needs FMA hardware
grep -qw fma /proc/cpuinfo 2>/dev/null || { echo "this CPU has no FMA instruction; the paper's build requires one (x86-64, Intel Haswell / AMD Piledriver or later)"; exit 1; }
gcc -O2 -mfma -ffp-contract=off -c $C/exp.c -o build/cr_exp.o
gcc -O2 -mfma -ffp-contract=off -I$C -c $C/log.c -o build/cr_log.o
for c in check_phi check_lift; do
    gcc -O2 -static -frounding-math -ffp-contract=off -I$C -DGITREV="\"$REV\"" -o build/$c src/$c.c build/cr_exp.o build/cr_log.o -lm
done
{ sha256sum build/check_phi build/check_lift; gcc --version | head -1; } | tee build/BINARY
echo "(the paper's build, gcc 11.4.0 on Ubuntu 22.04: check_phi 2ae354dc..., check_lift 94439b73...; another compiler gives other hashes, which is fine)"
