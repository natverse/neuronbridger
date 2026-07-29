# ColorMIP-search pipeline (Deng/Kondo LM → BANC neurons)

End-to-end pipeline for scoring every Deng/Kondo LM sample against the
BANC neuron colour-MIP library and emitting a ranked master CSV.

This runs **in parallel to** the layer-registration pipeline
(`{kondo,deng_brain,deng_vnc}_to_banc.R`), which targets BANC voxel
space. ColorMIP search needs the NeuronBridge template grids
(JRC2018_UNISEX_20x_HR for brain, JRC2018_VNC_UNISEX_HR for VNC), so
we re-warp the LM sources to those grids first.

## Stage 0 — One-time asset sync

All under `<data-root>/`:

| Path | Source | Rows/files |
|---|---|---|
| `banc_colormips/JRC2018_UNISEX_20x_HR/` | `gs://lee-lab.../neuron_colormips/template_alignment_240721/JRC2018_UNISEX_20x_HR/` | 93,724 PNGs |
| `banc_colormips/JRC2018_VNC_UNISEX_461/` | same GCS prefix, `JRC2018_VNC_UNISEX_461/` | 25,034 PNGs |
| `banc_meta/banc_888_meta.feather` | `gs://lee-lab.../compiled_data/banc_888/banc_888_meta.feather` | root_626 ↔ root_888 + cell metadata |
| `banc_meta/banc_banc_space_swc/<root_888>_l2.swc` | `gs://lee-lab.../compiled_data/banc_888/banc_banc_space_swc/` | L2 skeletons in BANC nm |

## Stage 1 — Query MIPs

`inst/scripts/make_colormip_library.R DATASET REGION [test|all] [--force]`

Re-warps each LSM (Deng) or NRRD (Kondo) to the NB target grid and
renders the search-ready colour-MIP:

- brain → `JRC2018_UNISEX_20x_HR` (1210 × 566 × 174 @ 0.519/0.519/1.0 µm)
- VNC   → `JRC2018_VNC_UNISEX_HR` (573 × 1119 × 219 @ 0.461/0.461/0.7 µm; `nrrd_to_mip(target_space="VNC")` adds the standard 90-px header)

Outputs per sample:
- `<dataset>_to_banc/colormip_library/<region>/<name>.nrrd` — the warped LM volume (needed downstream for NBLAST + voxel_attr)
- `<dataset>_to_banc/colormip_library/<region>/<name>.png`  — the query MIP

```bash
Rscript inst/scripts/make_colormip_library.R deng  brain all
Rscript inst/scripts/make_colormip_library.R deng  VNC   all
Rscript inst/scripts/make_colormip_library.R kondo brain all
```

## Stage 2 — Per-query search + augment

`inst/scripts/colormip_search_against_banc.R QUERY_DIR REGION WARPS_ROOT OUT_DIR [TOP_K] [--force]`

Per query MIP, in-process:

1. **`colormip_search()`** — query PNG × BANC library
   (`threshold=100`, `z_tolerance=2`, `xy_shift=2`, `mirror=TRUE`,
   `mc.cores=8`). Keeps top-K (default 25) BANC hits with `cm_score`,
   `cm_n_match`, `cm_dx/dy`, `cm_mirror`.
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
   backend hit). Column: `nblast`.
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

## End-to-end example

```bash
DATA="$HOME/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"

# Stage 1 (Deng brain)
Rscript inst/scripts/make_colormip_library.R deng brain all

# Stage 2 + 3 (Deng brain)
Rscript inst/scripts/colormip_search_against_banc.R \
  "$DATA/deng_to_banc/colormip_library/brain"  brain \
  "$DATA/deng_to_banc/colormip_library/brain"  \
  "$DATA/deng_to_banc/colormip_hits"  25
```

## Scaling + cost

- Search alone: ~140 min per brain sample (93,724-MIP library),
  ~45 min per VNC sample (25,034 lib), `mc.cores=8` on Apple silicon.
- NBLAST per hit + voxel_attr are seconds each once SWCs are cached.
- Master CSV write is trivial.
- Speedup levers still open (per `tasks.md`): smaller `xy_shift`
  (2→1 or omit), tighter pre-filter threshold, or pre-computing a
  pixel→depth-LUT index per library MIP (cache).
