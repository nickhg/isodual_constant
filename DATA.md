# The certificate data

The certificate logs total about 263 GB in about 70,000 files. They are published as one dataset with the layout below. Every log is listed, with the SHA-256 recorded when it was produced, in this repository's `manifests/`, so a copy from any source can be checked against them. A `SHA256SUMS` file at the top of the dataset covers every file in it. It is mirrored in this repository as `manifests/dataset_SHA256SUMS`, and so is archived with each Zenodo release, where it can't be changed.

## Where to get it

The dataset is on Hugging Face: [`nickhg1971/isodual-constant-certificates`](https://huggingface.co/datasets/nickhg1971/isodual-constant-certificates). Downloads are free and need no account. You can take any single file, any folder, or everything.

- **One file:** every file has a direct link:
  ```
  wget https://huggingface.co/datasets/nickhg1971/isodual-constant-certificates/resolve/main/eta5_stage2/logs/409/409_480_1_8_c023.logs.tar.xz
  ```
- **Folders or everything:** use the `hf` tool (`pip install huggingface_hub`), which resumes interrupted downloads:
  ```
  hf download nickhg1971/isodual-constant-certificates --repo-type dataset --include "eta5_stage1/391/*" --local-dir DATA
  hf download nickhg1971/isodual-constant-certificates --repo-type dataset --local-dir DATA     # all 263 GB
  ```

## Layout

`<slice>` is the first number of a job's name: stage-1 job `391_480_4_8_s01110` is in `eta5_stage1/391/`.

| path | files | size | contents |
|---|---|---|---|
| `phi5_lt7/` | 575 | 0.6 GB | the Φ₅ < 7 logs (`slice_<k>_480.log`) |
| `phi5_le_P/` | 4,757 | 1.7 GB | the Φ₅ ≤ 2¹⁶/(3⁴·5³) logs (`<tag>.log.gz`) |
| `eta5_stage1/<slice>/` | 24,133 × 2 | 1.1 GB | stage-1 logs (`<tag>.log.gz`) and their collected boxes (`<tag>.col.xz`) |
| `eta5_stage1_merge.tar` | 1 | 1.5 GB | for each of the 742 pieces, the merged stage-2 roots (`merged/<piece>.col`) and merge trees (`merged/<piece>.trees`) |
| `eta5_stage2/jobs.txt` | 1 | | the 15,841 stage-2 jobs, "job chunk" per line |
| `eta5_stage2/chunks.tar` | 1 | 1.4 GB | the merged boxes of each job (`chunks/<job>.col`) |
| `eta5_stage2/logs/<slice>/` | 15,839 | 257 GB | stage-2 logs, one `<job>.logs.tar.xz` per job (the largest 78 MB) |
| `eta5_stage2/fix/` | 3 | 70 MB | the box-by-box proofs of Appendix E (`<job>.fix.logs.tar.xz`) |
| `checker_outputs/` | 4 | 41 MB | the outputs of the paper's final re-check, one archive per certificate (`phi5_lt7.tar.xz` 45 KB, `phi5_le_P.tar.xz` 0.4 MB, `eta5_stage1.tar.xz` 2.4 MB, `eta5_stage2.tar.xz` 38 MB) |

To run the coverage scripts, unpack the two tars beside `jobs.txt`: `cover_s2.py` expects `merged/`, `chunks/` and `jobs.txt` in one directory (see `scripts/cover_s2.py`).

## What each part is for

- **The logs:** each is a search tree whose leaves a checker re-proves (see the paper, Appendices C–E).
- **The collected boxes, merge trees and chunks:** the bookkeeping that ties stage 1 to stage 2 exactly.
  - `scripts/merge_check.py` rebuilds every merged root from the collected boxes.
  - `scripts/cover_s2.py` checks that the chunks of each piece concatenate to its merged roots, and that every box has accepted logs.
- **The checker outputs:** not needed to verify anything, since you produce your own. They let you compare your outputs line by line with ours. For the stage-2 archive of the worked example, our output with a different compiler (gcc 12) was byte-identical to the stored one.

## Licence

CC BY 4.0 (see `LICENSE-DATA`).
