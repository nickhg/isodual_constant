# The isodual constant of the LLL potential: code and certificates

Code, certificates and verification instructions for

> **Nick Howgrave-Graham, *The isodual constant of the LLL potential: exact values, structure and reduction guarantees* (2026).** [arXiv link to be added]

**Status:** preprint. DOIs: this repository, release v1.0 (the version of the paper), [10.5281/zenodo.23650217](https://doi.org/10.5281/zenodo.23650217) (all versions: [10.5281/zenodo.23650216](https://doi.org/10.5281/zenodo.23650216)); the dataset, [10.57967/hf/11260](https://doi.org/10.57967/hf/11260). To report a problem, or a failed check, please open an issue on this repository.

## The paper

The constant η_n is the supremum over n-dimensional lattices of the LLL potential Ψ² = ∏ det Gram(b₁,…,b_k), minimized over bases and normalized. It is isodual: η(L) = η(L\*).

The paper proves:
- **η₅ = 6**, attained only by D₅ and D₅\*;
- **η₈ = 2304**, attained by E₈;
- **the uniqueness of D₄** for η₄ = 3;
- **η of every root lattice**;
- **log η_n = (¼ + o(1)) n² log n**;
- **reduction guarantees** for potential-based lattice reduction, which for blocks of size m ≤ 5 match slide reduction's Hermite factor.

The proof of η₅ = 6 is computer-assisted. This repository contains everything needed to check that part.

**To verify the results, start with [`VERIFY.md`](VERIFY.md).** It describes four tiers, from a laptop-scale check of every exact certificate up to a full re-check of all logs. For each it gives the compute, data transfer and cost needed, and what exactly is established.

## Layout

| directory | contents |
|---|---|
| `src/` | **The trusted checkers** `check_phi.c` and `check_lift.c`, with CORE-MATH's correctly rounded exp and log (`src/third_party/core-math`, MIT). **The untrusted producers** `certify.c` and `stage2.c`, at each version stamped in the logs (`src/producers/<stamp>/`). |
| `scripts/` | coverage and merge checks; the instance scripts used to produce and re-check the certificates; the broken-certificate tests; the CORE-MATH accuracy test |
| `certificates/` | the exact (non-search) certificates in SageMath: Theorem D5, Propositions pivd5 and mstar, M\*, L\*, the isodual geodesic, A₂⊗A₂, π_vA₄. `run_all.sh` runs them all; expected outputs are in `certificates/expected/` |
| `inputs/` | the two small fixed inputs of `check_phi`: the bases defining the balls about π_vD₅ (`pivd5_U.txt`, stage 1) and M\* (`mstar_U.txt`, Φ₅ ≤ P) |
| `independent/` | `check_lift2.py`: an independent stage-2 checker (Python + Arb ball arithmetic). `validate_gz.sage`: the coordinate maps checked from first principles |
| `manifests/` | the SHA-256 of every certificate log as recorded when it was produced; the split manifests and reference root boxes for coverage; coverage outputs; checker binary hashes |
| `example/` | a complete worked example (`example/run_example.sh`, about 3 minutes, data included) |

**The certificate logs themselves** (about 263 GB) are not in this repository. They are in the Hugging Face dataset [`nickhg1971/isodual-constant-certificates`](https://huggingface.co/datasets/nickhg1971/isodual-constant-certificates) ([doi:10.57967/hf/11260](https://doi.org/10.57967/hf/11260), revision `96e9321`), free to download with no account; see [`DATA.md`](DATA.md).

## Version stamps

Every log records the version (a git commit of the development repository) of the program that wrote it. Every checker output records the checker's version. The stamps and the corresponding files here:

| stamp | file | used for |
|---|---|---|
| (none) | an earlier certify.c (not archived; the producer is untrusted) | Φ₅ < 7 (its logs predate stamping) |
| `1399bb7` | `src/producers/1399bb7/certify.c` | Φ₅ ≤ 2¹⁶/(3⁴·5³) |
| `d8f5c29` | `src/producers/d8f5c29/certify.c` | η₅ = 6, stage 1 |
| `9d7c51a`, `59fdb6a` | `src/producers/<stamp>/stage2.c` | η₅ = 6, stage 2 (the producer was fixed mid-run; both are accepted, since soundness rests on the checker) |
| `b5e57c2` | `src/check_phi.c`, `src/check_lift.c` | the final re-check of every log of all four certificates (CORE-MATH exp/log), so the trusted side is versioned for each, including Φ₅ < 7 |

The SHA-256 of each source file is in `manifests/sources.sha256`.

## Citing

Please cite the paper (see `CITATION.cff`).

## Licences

- **Code:** MIT (see `LICENSE`).
- **Data and manifests:** CC BY 4.0 (see `LICENSE-DATA`).
- **CORE-MATH:** MIT, by its authors.
