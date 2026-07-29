# ColorMIP-search pipeline (Deng/Kondo LM → BANC neurons)

End-to-end pipeline for scoring every Deng/Kondo LM sample against the
BANC neuron colour-MIP library and emitting a single ranked CSV.

This runs **in parallel to** the layer-registration pipeline
(`{kondo,deng_brain,deng_vnc}_to_banc.R`), which targets BANC voxel
space. ColorMIP search needs the NeuronBridge template grids
(JRC2018_UNISEX_20x_HR for brain, JRC2018_VNC_UNISEX_HR for VNC), so
we re-warp the LM sources to those grids first.

## Stage 0 — One-time asset sync

Keep the BANC colormip library off Dropbox (~60 GB, no need to sync
between machines) — symlink from data-root:

```bash
mkdir -p "$HOME/banc_colormips/JRC2018_UNISEX_20x_HR"  \
         "$HOME/banc_colormips/JRC2018_VNC_UNISEX_461"
ln -s   "$HOME/banc_colormips/JRC2018_UNISEX_20x_HR"    \
        "$DATA/banc_colormips/JRC2018_UNISEX_20x_HR"
ln -s   "$HOME/banc_colormips/JRC2018_VNC_UNISEX_461"   \
        "$DATA/banc_colormips/JRC2018_VNC_UNISEX_461"
```

Assets (source paths canonicalise to `<data-root>/`):

| Path | Source | Files |
|---|---|---|
| `banc_colormips/JRC2018_UNISEX_20x_HR/` | `gs://lee-lab.../neuron_colormips/template_alignment_240721/JRC2018_UNISEX_20x_HR/` | 128,588 PNGs (~800 MB) |
| `banc_colormips/JRC2018_VNC_UNISEX_461/` | same GCS prefix, `JRC2018_VNC_UNISEX_461/` | 25,531 PNGs (~280 MB) |
| `banc_colormips/{brain,VNC}_cache_thr100.rds` | `Rscript build_colormip_index.R <region>` | one-shot sparse index (463 MB brain / 241 MB VNC) |
| `banc_meta/banc_888_meta.feather` | `gs://lee-lab.../compiled_data/banc_888/banc_888_meta.feather` | root_626 ↔ root_888 + cell metadata |
| `banc_meta/banc_banc_space_swc/<root_888>_l2.swc` | `gs://lee-lab.../compiled_data/banc_888/banc_banc_space_swc/` | L2 skeletons in BANC nm (only ~10% of hit ids have SWCs; expected) |

**Cache build:**

```bash
Rscript inst/scripts/build_colormip_index.R brain  --threshold 100 --mc.cores 8
Rscript inst/scripts/build_colormip_index.R VNC    --threshold 100 --mc.cores 8
```

Takes ~15 min for brain, ~5 min for VNC. Re-run only if you change the
search threshold or the library changes upstream.

## Stage 1 — Query MIPs

`inst/scripts/make_colormip_library.R DATASET REGION [test|all] [--force] [--threshold X]`

Re-warps each LSM (Deng) or NRRD (Kondo) to the NB target grid and
renders the search-ready colour-MIP:

- brain → `JRC2018_UNISEX_20x_HR` (1210 × 566 × 174 @ 0.519/0.519/1.0 µm)
- VNC   → `JRC2018_VNC_UNISEX_HR` (573 × 1119 × 219 @ 0.461/0.461/0.7 µm; `nrrd_to_mip(target_space="VNC")` adds the standard 90-px header)

**`--threshold` matters a lot.** Default `0.999` (top 0.1% of nonzero
voxels) makes query MIPs sparse enough to match the library sparsity —
dense receptor stains at Triangle-threshold produce 40-80% foreground
vs library MIPs' ~0.001%, drowning the scorer. Verified: top-5 hits
stable, cm_scores ~20× higher, per-query search ~10× faster.

Outputs per sample:
- `<dataset>_to_banc/colormip_library/<region>/<name>.nrrd` — the warped LM volume (needed downstream for NBLAST + voxel_attr)
- `<dataset>_to_banc/colormip_library/<region>/<name>.png`  — the query MIP

```bash
Rscript inst/scripts/make_colormip_library.R deng  brain all
Rscript inst/scripts/make_colormip_library.R deng  VNC   all
Rscript inst/scripts/make_colormip_library.R kondo brain all
```

## Stage 2 — Per-query search + augment

`inst/scripts/colormip_search_against_banc.R QUERY_DIR REGION WARPS_ROOT OUT_DIR [TOP_K] [--force] [--cache <rds>]`

Per query MIP, in-process:

1. **`colormip_search()`** — query PNG × BANC library
   (`threshold=100`, `z_tolerance=2`, `xy_shift=2`, `mirror=TRUE`,
   `mc.cores=8`). Keeps top-K (default 25) BANC hits with `cm_score`,
   `cm_n_match`, `cm_dx/dy`, `cm_mirror`. Pass `--cache <rds>` for a
   ~5× speedup by loading the pre-built sparse library index instead of
   decoding + indexing each PNG.
2. **v888 join** — left-join `banc_id` (which is `root_626`) against
   `banc_888_meta.feather` for `root_888`, `super_cluster`,
   `cell_function`, `neurotransmitter_predicted`.
   `preferred_id = root_888 %||% banc_id`.
3. **NBLAST vs local SWC** — LM dotprops from the warped NRRD (voxels
   > 20 intensity, subsample 50k, k=5). Per hit: read
   `<preferred_id>_l2.swc`, warp BANC nm → target template via
   `bancr::banc_to_JRC2018F(method="tpsreg", region="brain"|"vnc")`
   + `xform_brain(JRC2018F → JRC2018U)`, dotprops, `nblast(normalised=TRUE)`.
   Skips neurons with missing SWC or <6 dotprops points (no CAVE
   backend hit). Column: `nblast`. Note: only root_888 IDs have SWCs
   on GCS, so root_626-only hits stay `NA`. Expected coverage ~10%
   until v888 catches up.
4. **Voxel attribution** — forward-warp every LM voxel µm → BANC nm
   via `bancr::{jrc2018f,jrcvnc2018f}_to_banc_tpsreg` (`Morpho::tps3d`),
   then FNN k=1 NN against pooled L2 point cloud of the top-K
   skeletons. Per-hit vote fraction = `voxel_attr`.
5. **Combined ranking** — `mean_rank = mean(rank(cm_score),
   rank(nblast), rank(voxel_attr))`. NA scores get −Inf so
   NBLAST/attribution failures rank last, not drop.

Writes `<OUT_DIR>/<name>_hits.csv` (19 columns; character-typed
`banc_id`, `root_888`, `preferred_id` to preserve 18-digit precision).

## Stage 3 — Master CSV

The same script's tail: `bind_rows(list.files("*_hits.csv"))` →
`<OUT_DIR>/lm_to_banc_colormip_hits.csv`, sorted by `mean_rank`.

## Stage 4 — NBLAST + voxel_attr fill-in (optional but recommended)

`inst/scripts/nblast_fill.R REGION HITS_DIR WARPS_ROOT [--data-root P]`

The Stage 2 search leaves `nblast` = NA for any hit whose L2 SWC isn't
already cached under `banc_meta/banc_banc_space_swc/`. This script:

1. reads every `_hits.csv` under `HITS_DIR`
2. collects the union of `preferred_id`s
3. bulk-syncs missing `<id>_l2.swc` from GCS (missing IDs — root_626-only —
   are logged and continued past; not fatal)
4. rewrites each per-sample CSV with NBLAST + voxel_attr filled where
   possible, using the exact same code path as Stage 2
5. refreshes `lm_to_banc_colormip_hits.csv`

Run once per region:

```bash
Rscript inst/scripts/nblast_fill.R brain "$DATA/deng_to_banc/colormip_hits"  \
        "$DATA/deng_to_banc/colormip_library/brain"
Rscript inst/scripts/nblast_fill.R VNC   "$DATA/deng_to_banc/colormip_hits"  \
        "$DATA/deng_to_banc/colormip_library/VNC"
```

The brain run also handles VNC CSVs' SWC sync (single union of IDs), so
the VNC run is fast (fill only).

## End-to-end example

```bash
DATA="$HOME/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"

# Stage 0 — one-time cache build (after gsutil rsync of the library)
Rscript inst/scripts/build_colormip_index.R brain --threshold 100 --mc.cores 8
Rscript inst/scripts/build_colormip_index.R VNC   --threshold 100 --mc.cores 8

# Stage 1 — query MIPs (threshold=0.999 default; --force to regenerate)
Rscript inst/scripts/make_colormip_library.R deng brain all
Rscript inst/scripts/make_colormip_library.R deng VNC   all

# Stage 2 + 3 — search + master CSV (cached)
Rscript inst/scripts/colormip_search_against_banc.R  \
  "$DATA/deng_to_banc/colormip_library/brain"  brain \
  "$DATA/deng_to_banc/colormip_library/brain"        \
  "$DATA/deng_to_banc/colormip_hits"  25             \
  --cache "$DATA/banc_colormips/brain_cache_thr100.rds"

Rscript inst/scripts/colormip_search_against_banc.R  \
  "$DATA/deng_to_banc/colormip_library/VNC"    VNC   \
  "$DATA/deng_to_banc/colormip_library/VNC"          \
  "$DATA/deng_to_banc/colormip_hits"  25             \
  --cache "$DATA/banc_colormips/VNC_cache_thr100.rds"

# Stage 4 — NBLAST + voxel_attr fill-in + master refresh
Rscript inst/scripts/nblast_fill.R brain  \
  "$DATA/deng_to_banc/colormip_hits"      \
  "$DATA/deng_to_banc/colormip_library/brain"
Rscript inst/scripts/nblast_fill.R VNC    \
  "$DATA/deng_to_banc/colormip_hits"      \
  "$DATA/deng_to_banc/colormip_library/VNC"
```

## Scaling + cost (as measured 2026-07-28 on Apple silicon, 16 cores)

- Cache build: 15 min brain (128k), 5 min VNC (25k). One-shot.
- Query MIP (Stage 1): ~40 s/sample re-MIP-only, ~5 min/sample if
  re-warping is needed too.
- Search per query (Stage 2 with cache): ~5 min brain, ~1 min VNC.
- NBLAST + voxel_attr fill per sample: seconds.
- **End-to-end for Deng cohort (54 brain + 39 VNC)**: ~5 hr brain
  search + ~40 min VNC search + ~10 min fill. Compared to the
  pre-speedup estimate of ~5-6 days for just the search.

Segfaults under `mclapply` at `mc.cores=8` are transient (macOS +
Accelerate.framework fork issue) — the outer tryCatch reports FAILED
and moves on. If a specific sample fails repeatedly, retry it single-
threaded with BLAS thread-limit env vars:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
MKL_NUM_THREADS=1  Rscript ... mc.cores=1
```
