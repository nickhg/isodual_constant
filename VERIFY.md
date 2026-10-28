# Verifying the results

This guide is for referees and readers who want to check the computer-assisted parts of

> Nick Howgrave-Graham, *The isodual constant of the LLL potential: exact values and structure* (2026).

The main result, η₅ = 6 (attained only by D₅ and D₅\*), has two kinds of ingredients:

1. **Proofs on paper:** the reduction to lattices with a shortest vector, the root box, the meaning of each leaf type, and the results the leaves invoke.
2. **Machine-checked certificates:** about 9·10⁹ leaves of logged branch-and-bound searches, each re-proved by an interval checker, plus coverage scripts showing that the leaves tile the search space.

Checking can be done at four levels of effort. Each level says what it establishes, what it costs, and what it still leaves you to trust. Stop at whichever level convinces you.

| Tier | What you check | Time | Machine | Download |
|---|---|---|---|---|
| 1 | every **exact** (non-search) certificate | ~3.5 h total (one core) | laptop with SageMath | none (in this repo) |
| 2 | a complete worked example: checker → coverage, plus the checker's own tests | ~3 min | laptop | none (in this repo) |
| 3 | a **sample** of the published certificate, against the hash manifests | ~1–2 h | laptop / workstation | a few GB |
| 4 | **every** log, with one static checker build, plus full coverage | ~620 vCPU-hours | cloud or cluster | 263 GB |

## The trusted base

Whatever tier you choose, the following are trusted. They are listed in the paper's paragraph "What has to be trusted" (end of Appendix E).

- **The checkers:** `src/check_phi.c` (696 lines) and `src/check_lift.c` (557 lines).
  - Their exp and log come from **CORE-MATH** (`src/third_party/core-math`), which is correctly rounded in every rounding mode, widened by one ulp.
  - They rely on IEEE-754 arithmetic under directed rounding, and on the C compiler.
- **For the exact certificates (tier 1):** SageMath, through its exact rational arithmetic and, where intervals are used, its interval arithmetic (`RealIntervalField`). Theorem 13 is one of these, and the stage-2 (B) leaves rest on it.
- **The coverage and merge scripts:** `scripts/cover_check.py`, `root_box_check.py`, `merge_check.py` and `cover_s2.py`, about 240 lines of Python comparing hexadecimal floats exactly.
- **The reductions on paper.**

**Leaf-type names.** The paper calls the two leaf types of `check_phi` (D) and (A) (Appendix C). The header comment of `src/check_phi.c` calls the second one (C); the source is left unchanged so that it stays byte-identical to the version stamped in the checker outputs (`manifests/sources.sha256`).

**Not trusted:** the producers (`src/producers/*/`) and all orchestration scripts. A log is only a *candidate* until a checker has proved every leaf in it.

**An independent second checker** (`independent/check_lift2.py`, in Python with Arb ball arithmetic, written from the log specification alone) re-proves a large sample (see tier 3). It gives evidence against a common-mode bug in `check_lift`.

---

## Tier 1: the exact certificates (laptop, SageMath)

```
cd certificates && ./run_all.sh
```

| script | claim in the paper | expected output (see `certificates/expected/`) | time |
|---|---|---|---|
| `d5_lemma.sage 1/10 5` | **Theorem 13** (used in the η₅ proof): η ≤ 6 within distance 1/10 of D₅, equality only at D₅; likewise for D₅\* | `LEMMA HOLDS: eta(G) <= 6 for d(G, G0) <= 1/10 …` | 3 min |
| `radius_d4d5.sage` | Appendix A: the D₄ radius r = 2/3 exactly (all 105,184 facets, via Normaliz), so η ≤ 3 within distance 33/100 of D₄; and r ≤ 4/(3√5) at D₅ | `LEMMA HOLDS: eta(G) <= 3 for d(G, G0) <= 33/100 …` | 1 min (needs `normaliz`) |
| `rad_produce.py cert 8`, then `rad_check.py cert` | Appendix A: the D₅ radius is exactly 4/(3√5) (lower bound; the upper bound is in `radius_d4d5.sage`). The producer (numpy, scipy; untrusted) writes a cover of 448 logs; the checker re-proves every box in exact rationals. Inputs come from `rad_inputs.sage` | `THEOREM HOLDS: … r(D5) >= 4/(3 sqrt5)`, with 199,619 boxes | 10 min on 8 cores + 1 min |
| `pivd5_lemma.sage 1/20 11/20` | Prop. 21: π_vD₅ is a sharp local maximum of Φ₅ (landscape only, not used for η₅) | interior and exterior margins 0.0277, 0.1276 | 62 min |
| `mstar_lemma.sage 1/40 9/10` | Prop. 22: M\* is a sharp local maximum of Φ₅ (strictness in Thm. 20; not used for η₅) | margins 0.0691, 0.0460 | 117 min |
| `mstar_unique.sage`, `lstar_check.sage` | κ(M\*) attained only at ±φ₀; η(L\*) = 486/125 | as stated | minutes |
| `mstar_invariants.sage`, `family_exact.py`, `geodesic_exact.py` | \|Aut(M\*)\| = 72, not strongly modular; M\* lies on the isodual geodesic D₄ → A₂⊗A₂ at parameter log₂(6/5) | exact | minutes |
| `a2a2.sage`, `shadow_iso.sage` | the exact entries of Remark 24 | exact | minutes |
| `mstar_genus.sage`, `mstar_mass.sage` | Prop. 17 (7): the genus of M\* has class number 6 (mass check Σ1/\|Aut\| = 1/3), and M\* is its unique lattice of largest minimum | `class number of the genus: 6`; minima 4, 16, 10, 4, 12, 10 | 1 min |
| `mstar_eisenstein.sage`, `mstar_herm_check.sage`, `mstar_herm_genus.sage` | Prop. 17 (8) and §6: M\* is Hermitian of rank 2 over ℤ[ω] (Gram matrix (8, 3−2ω; 3−2ω̄, 8), det 45); the four structured classes of its genus | exact | 1 min |
| `d4_kappa_classes.py` | Prop. 12: exact search over the torsion grids finds exactly 48 κ-maximizer classes at D₄ and 24 at A₄ | `classes attaining 1/h: 48` (D₄), `24` (A₄) | 4 min |
| `e6_flags.sage`, `e6_cross.sage`, `e7_design.sage` | Prop. 26: the Ψ-minimal flags of E₆ (19,440) and E₇ (544,320) and the design identities, exactly; the LP certificate r(E₆) ≥ 3/10 | exact | 6 min |
| `schur_check.sage 6`, `7`, `8` | Prop. 26: W(E_n) acts irreducibly on traceless symmetric forms (commutant of dimension 1); for E₈ the interval check of the radius 1/20 | `commutant dim: 1` | 1 min |

**What this establishes:** all exact claims that do not rest on a search.
- **For η₅:** this covers Theorem 13, which is the equality case and the (B) leaves of stage 2.
- **Arithmetic:** these computations use exact rational arithmetic, or interval arithmetic in Sage with rigorous rounding.

## Tier 2: a complete worked example (laptop, minutes)

```
example/run_example.sh
```
It needs gcc and python3, plus python-flint and mpmath for steps 5 and 7 (skipped with a message if missing). The data is in `example/data/` (3.9 MB), and every log in it is listed in the manifests. The example:
1. builds both checkers exactly as in the final re-check (`scripts/build.sh`: `gcc -O2 -frounding-math -ffp-contract=off`, static, CORE-MATH with `-mfma`; see "Build requirements" below);
2. **stage 1:** re-checks one stage-1 log with `check_phi` and requires its collected boxes (the inputs of stage 2) to be byte-identical to the stored ones;
3. **stage 2:** re-checks one stage-2 archive (97 logs, 114,074 leaves, one of them (B)) with `check_lift`, and requires the output to equal the final re-check's byte for byte;
4. **coverage**, for the stage-1 piece 477_480_4_8:
   - its 32 collected-box files must match the manifest;
   - `merge_check.py` rebuilds its 165 merged stage-2 roots as exact unions of all 554 C leaves;
   - `cover_s2.py` checks that every root is proved by an accepted stage-2 log. That job's archive is 67 MB and not included, so by default this step reads the final re-check's output for it; set `ARCHIVE_477=<file>` to re-check the archive itself;
5. runs the independent checker `check_lift2.py` on the (B) leaf and a 0.2% sample of one log;
6. runs `scripts/check_lift_tests.sh`: 23 deliberately broken certificates, each of which the checker must reject;
7. runs `scripts/test_crexplog.sh`: CORE-MATH's exp and log against mpmath at 60 digits, in both rounding modes (20,000 points).

Every step prints `-> passed`, and the script ends with `EXAMPLE PASSED`, after about 3 minutes on one laptop core. We ran it with gcc 12.2 (Debian 12) on an Intel Core i3-1215U. The paper's final re-check used gcc 11.4 (Ubuntu 22.04), and the stage-2 checker output was still byte-identical to the final re-check's. This is a robustness check, not a requirement: any C99 compiler with IEEE-754 doubles and directed rounding (`-frounding-math`) should do.

**Build requirements.** An x86-64 CPU with the FMA instruction: Intel from Haswell (2013) on, or AMD from Piledriver (2012) on. `scripts/build.sh` checks this and stops if it is missing. `-mfma` is used only for the two CORE-MATH files (`exp.c`, `log.c`), which call `fma` explicitly. The checkers themselves are compiled without `-mfma` and with `-ffp-contract=off`, so their own arithmetic contains no fused operations. This is the build used for every check reported in the paper.

**What this establishes:**
- **It does:** the pipeline works as described; the log format is what the paper says; the checker rejects tampered logs; the arithmetic library is correctly rounded on the test points.
- **It does not:** say anything about the rest of the certificate.

## Tier 3: a sample of the published certificate (workstation, a few GB)

**Get the data** (see `DATA.md`): download individual files or folders from Hugging Face, free and with no account needed.

**Check any files you like.** Here DATA is where you put the dataset, keeping its layout. For each file, `check_one.sh` does three things:
- it compares the file's SHA-256 with the one recorded at production (`manifests/`);
- it runs the checker and requires `FAILURES 0`;
- for a stage-1 log, it also requires the collected boxes to match.

It prints one line per file, ending `OK` or `FAILED: <reason>`.
```
scripts/build.sh
scripts/check_one.sh DATA/eta5_stage2/logs/<slice>/<job>.logs.tar.xz   # check_lift on all its logs
scripts/check_one.sh DATA/eta5_stage1/<slice>/<tag>.log.gz             # check_phi (stage 1); needs <tag>.col.xz beside it
scripts/check_one.sh DATA/phi5_le_P/<tag>.log.gz               # check_phi (Phi_5 <= P, with the M* balls)
scripts/check_one.sh DATA/phi5_lt7/slice_<k>_480.log           # check_phi (Phi_5 < 7)
```

**Indicative cost:**

| check | download per item | time per item | 100 random items |
|---|---|---|---|
| stage-2 archive | 16 MB on average (3.7 MB in the example) | ~1.5 min on average (21 s for the example's) | 1.6 GB, ~2.5 h on one core |
| stage-1 log (+ its collected boxes) | ~50 KB | ~0.5–1 min (47 s for the example's) | 5 MB, ~1–1.5 h |

Times (here and in tier 1) are measured on one core of a laptop (Intel Core i3-1215U); a cloud core is somewhat faster. They scale with the number of leaves, which the checker prints.

**The independent checker:**
```
python3 independent/check_lift2.py --types B <logs>                 # every (B) leaf
python3 independent/check_lift2.py --types NFP --sample 0.01 <logs>   # 1% of the other leaves
```
This needs python-flint (`pip install python-flint`), and is about 0.2 s per leaf. For the paper we ran it on **every** (B) leaf (47,437 leaves, in 304 archives), and on a seeded 0.5% sample of all other leaves in those archives and in 300 random others. The results are in `manifests/check_lift2_summary.txt`. `independent/validate_gz.sage` checks the coordinate maps from first principles: at 100,000 random points it builds explicit bases and compares their Gram matrices, agreeing to 10⁻⁵⁷.

**What this establishes:**
- **For the files you check:** they are the ones produced and checked (the hashes), and every leaf in them is proved.
- **For a uniform random sample of k leaves of the certificate:** if a fraction p of all leaves were wrong, the sample would miss every wrong one with probability (1 − p)ᵏ.
- **The paper's independent-checker sample is stratified, not uniform.** It supports three statements:
  - **(B) leaves:** all 47,437 were checked.
  - **Archives:** 300 of the 15,537 stage-2 jobs without (B) leaves were drawn at random. If errors that check_lift2 would detect were confined to a fraction q of those jobs, the chance that none of the 300 is among them is (1 − q)³⁰⁰: about 5% for q = 1%, and 0.2% for q = 2%.
  - **Leaves within a drawn archive:** each non-(B) leaf was checked with probability 0.005. In those 300 archives, 736,494 leaves were checked. If a fraction p of the leaves of an archive with n leaves were wrong, all of them would be missed with probability about e^(−0.005·p·n).

  The second checker is therefore evidence against a systematic (common-mode) error in `check_lift`, not a re-proof of every leaf. That is `check_lift`'s job, and tier 4 repeats it.

## Tier 4: everything (cloud or cluster, ~620 vCPU-hours, 263 GB)

This repeats the paper's final re-check: every log of the four certificates (Φ₅ < 7, Φ₅ ≤ P, and stages 1 and 2 of η₅ = 6), 45,307 jobs in all.

**Compute:**
- **Measured:** about 620 vCPU-hours on c7i/c6i-class cloud cores. Of this, about 260 is stage 2 and the rest is the three check_phi runs.
- **Indicative cost:** at 2026 spot prices, roughly $15–30 of compute.

**Data:** 263 GB from Hugging Face, free. 257 GB of it is the stage-2 logs; the rest is the other runs' logs (about 3.5 GB), the merge and chunk data (about 3 GB) and the checker outputs. Download into the cloud region you compute in; most clouds don't charge for inbound transfer. Hugging Face may rate-limit very large anonymous downloads, so a logged-in `hf` session (a free account, `hf auth login`) may be faster for all 263 GB.

**Requirements:** gcc, Python 3 with `mpmath` (`pip install mpmath`; `root_box_check.py` computes the box B with it), and the `hf` tool (`pip install huggingface_hub`).

**Procedure:**
1. **Download and build:**
   ```
   hf download nickhg1971/isodual-constant-certificates --repo-type dataset --local-dir DATA
   (cd DATA && sha256sum --quiet -c SHA256SUMS)
   scripts/build.sh
   ```
   The paper's build has SHA-256 `2ae354dc…` (check_phi) and `94439b73…` (check_lift), with gcc 11.4.0 on Ubuntu 22.04. A different compiler gives a different hash, which is fine.
2. **Re-check every log** (45,307 of them), in parallel:
   ```
   find DATA \( -path '*/phi5_lt7/*.log' -o -path '*/phi5_le_P/*.log.gz' -o -path '*/eta5_stage1/*.log.gz' \
        -o -name '*.logs.tar.xz' \) | xargs -P $(nproc) -n 1 scripts/check_one.sh > results.txt
   grep -c ' OK$' results.txt          # must be 45307
   ```
   For each log, `check_one.sh` compares its SHA-256 with the manifest recorded at production, runs the checker and requires `FAILURES 0`. For stage 1 it also requires the collected boxes to be byte-identical to the stored ones.
3. **Coverage:** `scripts/coverage_all.sh DATA`. It runs:
   - `cover_check.py` with the split manifests, for Φ₅ < 7, Φ₅ ≤ P and stage 1, plus `root_box_check.py`;
   - `merge_check.py`, which rebuilds every merged stage-2 root from the stage-1 collected boxes;
   - `cover_s2.py` on all 742 pieces, using your outputs from step 2.

   Expected: "union contains B", "matches the recorded output" three times (the recorded outputs are in `manifests/coverage/`), "merge check OK on 742 pieces" and "stage-2 coverage OK: 742 pieces". This takes about 7 minutes on 8 cores.

`scripts/replay2.sh` is the script the paper's own final re-check used, reading from our storage. It is included as a record; the steps above do the same from the public dataset.

**What this establishes:** every leaf of every certificate has been proved by a checker you built, and the leaves tile the search spaces. Together with the reductions on paper and the tier-1 certificates, that is the whole proof of Theorems 19, 20 and 14.
