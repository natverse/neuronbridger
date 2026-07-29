# `inst/scripts/`

Each script is a self-contained reproducer. The first comment block in
every file documents inputs/outputs/usage in detail; this README is the
top-level map.

## Pipeline scripts (LM → BANC + downstream)

These are the primary reproducible workflows. They use the `R/`
helpers (`lsm_to_nrrd`, `lm_to_jrc2018u_elastix`, `lsm_to_banc_layer`,
`lsm_vnc_to_banc_layer`, `nrrd_to_precomputed`, `colormip_search`,
etc.) and ship as the practical entry points for re-running each
data-import step.

| Script | What it does |
|---|---|
| `kondo_to_banc.R`         | Bridge **Kondo et al. 2020** receptor LM stacks (IS2 → FCWB → JRC2018F → BANC voxel) and upload `<gene>_<sample>.ng` to GCS. |
| `deng_brain_to_banc.R`    | Bridge **Deng et al. 2019 BRAIN** LSMs (native → JRC2018U → JRC2018F → BANC) and upload to GCS. |
| `deng_vnc_to_banc.R`      | Bridge **Deng et al. 2019 VNC** LSMs (native → JRC2018VNCF → BANC via the bancr tpsreg + Python inverse-warp) and upload to GCS. |
| `make_colormip_library.R` | Re-warp Deng/Kondo source volumes to the **NeuronBridge target template** (`JRC2018U_HR` for brain, `JRC2018VNCU_HR` for VNC) and render colour-MIP PNGs (default `--threshold 0.999`: top 0.1% of nonzero voxels for library-matching sparsity). |
| `build_colormip_index.R`  | One-shot: read every BANC library PNG, write a sparse `(fg_idx, z)` RDS keyed at `<data-root>/banc_colormips/<region>_cache_thr<T>.rds`. Consumed by `colormip_search_against_banc.R --cache` for a ~5× per-query speedup. |
| `colormip_search_against_banc.R` | Search a directory of query colour-MIPs against the BANC neuron MIP library (via `--cache` when the index exists), NBLAST top-K vs locally-cached L2 SWCs, voxel-attribute LM signal to nearest skeleton, emit `<name>_hits.csv` per query + a single `lm_to_banc_colormip_hits.csv` master. |
| `nblast_fill.R`           | Post-search: bulk-sync missing L2 SWCs from GCS, walk each `_hits.csv`, and fill in `nblast` + `voxel_attr` for whichever rows now have SWCs available. Refreshes `lm_to_banc_colormip_hits.csv`. |

See `colormip_pipeline.md` for the end-to-end walkthrough of the
colormip-search half (Stages 0–4).

### Common arguments

Every pipeline script accepts:

* `--data-root <path>`   root for every local source + output. Defaults
  to `~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat`
  (or the `$NEURONBRIDGER_DATA_ROOT` env var if set).
* `--force`              recompute even if the per-sample output
  already exists (default skips and prints a SKIP line).

The three `*_to_banc.R` registration scripts also accept:

* `--upload`     on completion, `gsutil -m cp -Z -r` each `*.ng/`
  directory into `gs://lee-lab.../light_level/<dataset>/`.
* `--register`   on completion, merge the per-batch
  `registry_entries.json` into the master `registry.json` on GCS,
  then re-mint `bancr::banc_lm_links` (runs
  `bancr/data-raw/make_banc_lm_links.R`). Implies `--upload`.
* `--one <name>` worker mode: process a single sample and write per-sample
  RDS stubs. The outer loop invokes this via `Rscript` so each sample
  runs in a fresh R heap, sidestepping the mid-batch macOS memory-
  pressure slowdown. Not for hand invocation; use `test`/`all` modes.

`colormip_search_against_banc.R` also accepts `--cache <rds>` (path
to a `build_colormip_index.R` output). `make_colormip_library.R` also
accepts `--threshold X` (passed through to `nrrd_to_mip()`; default
0.999).

### Typical run order

```bash
DATA="$HOME/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat"

# (one-time) sync the BANC neuron MIP libraries to a real dir OFF Dropbox
# (~1 GB total, avoids Dropbox sync churn), then symlink from data-root.
mkdir -p "$HOME/banc_colormips/JRC2018_UNISEX_20x_HR" \
         "$HOME/banc_colormips/JRC2018_VNC_UNISEX_461"
gsutil -m rsync -r \
  gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_colormips/template_alignment_240721/JRC2018_UNISEX_20x_HR/ \
  "$HOME/banc_colormips/JRC2018_UNISEX_20x_HR/"
gsutil -m rsync -r \
  gs://lee-lab_brain-and-nerve-cord-fly-connectome/neuron_colormips/template_alignment_240721/JRC2018_VNC_UNISEX_461/ \
  "$HOME/banc_colormips/JRC2018_VNC_UNISEX_461/"
mkdir -p "$DATA/banc_colormips"
ln -s "$HOME/banc_colormips/JRC2018_UNISEX_20x_HR"  "$DATA/banc_colormips/"
ln -s "$HOME/banc_colormips/JRC2018_VNC_UNISEX_461" "$DATA/banc_colormips/"

# (one-time) sync v888 metadata
mkdir -p "$DATA/banc_meta/banc_banc_space_swc"
gsutil cp \
  gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/banc_888_meta.feather \
  "$DATA/banc_meta/"
# L2 SWCs are synced on demand by `nblast_fill.R` at Stage 4.

# (one-time) build the library index cache for fast search
Rscript inst/scripts/build_colormip_index.R brain --threshold 100 --mc.cores 8
Rscript inst/scripts/build_colormip_index.R VNC   --threshold 100 --mc.cores 8

# 1. Register LM data + ship to GCS (skips per-sample outputs that exist)
Rscript inst/scripts/kondo_to_banc.R       all   --register
Rscript inst/scripts/deng_brain_to_banc.R  all   --register
Rscript inst/scripts/deng_vnc_to_banc.R    all   --register

# 2. Make NeuronBridge-compatible MIPs (--threshold 0.999 is the default)
Rscript inst/scripts/make_colormip_library.R deng  brain   all
Rscript inst/scripts/make_colormip_library.R deng  VNC     all
Rscript inst/scripts/make_colormip_library.R kondo brain   all

# 3. Search against BANC + emit master ranked CSV (with cache)
Rscript inst/scripts/colormip_search_against_banc.R \
  "$DATA/deng_to_banc/colormip_library/brain"  brain \
  "$DATA/deng_to_banc/colormip_library/brain"        \
  "$DATA/deng_to_banc/colormip_hits"  25             \
  --cache "$DATA/banc_colormips/brain_cache_thr100.rds"

# 4. NBLAST + voxel_attr fill-in (bulk-syncs L2 SWCs, refreshes master)
Rscript inst/scripts/nblast_fill.R brain "$DATA/deng_to_banc/colormip_hits" \
  "$DATA/deng_to_banc/colormip_library/brain"
Rscript inst/scripts/nblast_fill.R VNC   "$DATA/deng_to_banc/colormip_hits" \
  "$DATA/deng_to_banc/colormip_library/VNC"
```

### Data paths

All local paths are derived from `--data-root` (or the
`$NEURONBRIDGER_DATA_ROOT` env var; default
`~/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/`).

**Local sources** (read; relative to `<data-root>/`):

| What | Path |
|---|---|
| Deng raw LSMs (Bowen / Lee lab mirror) | `imaging-CCT-Bowen/<gene>-<...>-(BRAIN|VNC)-<sex>.lsm` |
| Kondo IS2 NRRDs (G-Node mirror)        | `kondo_et_al_2020/nrrd/IS2_<gene>_no<n>_02_warp_*.nrrd` |
| Templates                              | `templates/JRC2018_UNISEX_20x_HR.nrrd`, `JRC2018_VNC_UNISEX_HR.nrrd`, `JRC2018_VNC_FEMALE_461.nrrd` |
| BANC v888 metadata (synced once)       | `banc_meta/banc_888_meta.feather` |
| BANC v888 L2 SWCs (synced once)        | `banc_meta/banc_banc_space_swc/<root_888>_l2.swc` |
| BANC neuron colorMIP libraries (synced once, symlink from `~/banc_colormips/`) | `banc_colormips/{JRC2018_UNISEX_20x_HR,JRC2018_VNC_UNISEX_461}/` |
| Sparse library-index cache (`build_colormip_index.R`) | `banc_colormips/{brain,VNC}_cache_thr<T>.rds` |

**Local outputs** (write; relative to `<data-root>/`):

| What | Path |
|---|---|
| BANC-aligned `.ng` (Deng)              | `deng_to_banc/deng_et_al_2019/<gene>_<region>-<sex>.ng/` |
| BANC-aligned `.ng` (Kondo)             | `kondo_to_banc/kondo_et_al_2020/<gene>_no<n>.ng/` |
| Per-batch registry stub + timings      | `<dataset_out>/registry_entries.json`, `timings.csv` |
| NB-template MIPs                       | `{deng,kondo}_to_banc/colormip_library/{brain,VNC}/<name>.png` |
| Search hits per sample                 | `deng_to_banc/colormip_hits/<name>_hits.csv` |
| Master ranked CSV                      | `deng_to_banc/colormip_hits/lm_to_banc_colormip_hits.csv` |

**GCS** (`gs://lee-lab_brain-and-nerve-cord-fly-connectome/`):

| What | Path |
|---|---|
| Master LM registry                                   | `light_level/registry.json` |
| Per-dataset LM layers                                | `light_level/{kondo_et_al_2020,deng_et_al_2019}/<name>.ng/` |
| BANC neuron colour-MIP library                       | `neuron_colormips/template_alignment_240721/{JRC2018_UNISEX_20x_HR,JRC2018_VNC_UNISEX_461}/` |
| BANC v888 metadata (root_626 ↔ root_888 ↔ cell_type) | `compiled_data/banc_888/banc_888_meta.feather` |
| BANC v888 L2 skeletons (BANC nm)                     | `compiled_data/banc_888/banc_banc_space_swc/<root_888>_l2.swc` |

**R packages**:

| What | Where |
|---|---|
| Master LM-link table (213 rows = 120 Kondo + 54 Deng BRAIN + 39 Deng VNC as of 2026-07-08) | `bancr::banc_lm_links` |
| `bancr` repo (for re-minting links)                  | `~/Projects/flyconnectome/bancr` |
| Saalfeldlab JRC2018 H5 bridges (brain + VNC)         | `~/Library/Application Support/R/nat.jrcbrains/{JRC2018U_JRC2018F,JRCVNC2018U_JRCVNC2018F,...}/` (VNC ones symlinked from `~/flybrain-data/` where Python `flybrains.download_jrc_vnc_transforms()` writes them) |

## Vignette reproducers

These exist solely to regenerate the figures, tables, and cached
artefacts referenced from a specific vignette. None of them are part
of the LM-to-BANC pipeline.

| Script | Backs vignette | Notes |
|---|---|---|
| `run_asta_sez.R`              | `match_em_to_line.Rmd` | Drives the full SS32423 → AstA-SEZ workflow; caches NB hits in `inst/extdata/asta_sez/`. |
| `run_asta_sez_figs.R`         | `match_em_to_line.Rmd` | Re-renders SEZ figures from the cached RDS. |
| `asta_sez_4cell.R`            | `match_em_to_line.Rmd` | The canonical 4-cell SEZ-AstA panel (CB0602 + CB0108). |
| `asta_sez_meshes.R`           | `match_em_to_line.Rmd` | Mesh-based 2-D MIP renders for the candidate panels. |
| `asta_sez_mip_panel.R`        | `match_em_to_line.Rmd` | Side-by-side colour-MIP comparison: NB vs FlyWire renders. |
| `asta_sez_natggplot.R`        | `match_em_to_line.Rmd` | nat.ggplot brain-figure variant of the candidate panel. |
| `asta_sez_ss32423_montage.R`  | `match_em_to_line.Rmd` | All 16 SS32423 brain MIPs tiled into a montage. |
| `asta_other_neurons.R`        | `match_em_to_line.Rmd` | Non-SEZ AstA candidates panel. |
| `run_abdominal_peripheral.R`  | `find_peripheral_targets.Rmd` | MANC abdominal motor + endocrine neurons → NB hits → KDRC. |
| `colormip_methods_panel.R`    | `mip_from_em.Rmd` | Three-back-end colour-MIP comparison panel for AstA1. |
| `colormip_param_sweep.R`      | `search_banc_with_lm.Rmd` | Tune `colormip_search()` parameters against an NBLAST ground-truth ranking. |
| `banc_colormip_sren.R`        | `search_banc_with_lm.Rmd` | End-to-end SREN-vs-BANC colour-MIP search reproducer. |
| `lm_capar_to_precomputed.R`   | `publish_lm_layer.Rmd` | Convert one IS2 NRRD to a Neuroglancer precomputed layer. |

## Helpers

`R/lm_pipeline.R` exposes the per-sample orchestrators
`lsm_to_banc_layer()` (brain) and `lsm_vnc_to_banc_layer()` (VNC).
`inst/python/lm_vnc_inverse_warp.py` is the Python helper the VNC
orchestrator drives via `system2()` to inverse-warp + max-pool the
JRC2018VNCF GFP volume into BANC voxel space using the bancr tpsreg.
