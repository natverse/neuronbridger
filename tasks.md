# neuronbridger — work-in-progress tasks

A running record of the multi-session work on `nrrd_to_mip()` /
`nrrd_to_precomputed()` and the LM-on-Neuroglancer pipeline. The most
recent task is at the bottom.

## Done

- **Pure-R `nrrd_to_mip(method = "direct")`.** Replaces the FIJI launcher
  as the default code path; the FIJI macro is still callable via
  `method = "fiji"`. R port is byte-equivalent to Stephan Gerhard's
  Python port in `fanc.render_neurons.make_colormip` (verified
  numerically: 99.97 % pixels exactly identical, max diff = 1/255 RGB
  unit). Embedded the verbatim 256-entry Janelia depth LUT.
- **`nrrd_to_mip(method = "python")`.** Calls the BANC Python
  implementation directly via `reticulate` for byte-level validation.
  Auto-imports `banc.render_neurons` (or `fanc.render_neurons` from the
  source repo); falls back to the bundled LUT if the BANC python pkg
  is not installed.
- **`vignettes/colormip_direct_vs_fiji.Rmd`.** Three-back-end
  walkthrough; uses **BANC v888 AstA1** (`cell_type == "AstA1"`,
  `root_888` IDs `720575941506055874` and `720575941541909965`) as
  the example neuron. Cached sub-sampled point cloud at
  `inst/extdata/asta1/banc_asta1_points_nm.rds` (2.1 MB) so the
  reproducer runs without re-fetching draco meshes.
- **`inst/scripts/colormip_methods_panel.R` reproducer** for the
  three-row panel saved at `inst/images/colormip_methods_panel.png`.
- **`asta_sez.Rmd` follow-ups** (per user requests): removed
  *Iteration* section + BANC-colormips Python-fallback line; replaced
  the FIJI MIP-rendering example with `nrrd_to_mip_direct(...)` /
  `nrrd_to_mip(method="direct")` against `JRC2018U_HR`; promoted
  `CB0108` from "non-SEZ" to the canonical SEZ AstA pair alongside
  `CB0602`; updated tables / captions / render lists accordingly.
- **`nrrd_to_precomputed()`.** Pure-R wrapper around `cloud-volume`
  (via `reticulate`) that writes any 3-D image volume in
  Neuroglancer's "precomputed" format (info JSON + chunked, gzip'd
  raw bricks). Tested round-trip via `cv$CloudVolume(...)[]` reads.
- **`vignettes/lm_layer_neuroglancer.Rmd`.** Four-stage pipeline:
  IS2 → JRC2018F (CMTK + `nat.jrcbrains`; documented but `eval=FALSE`
  — Stage 1 install pending), JRC2018F → BANC voxel via Elastix
  (cited from BANC repo `template_to_BANC.txt`), precomputed write,
  scene assembly. Self-contained `eval=TRUE` sanity check writes a
  synthetic volume to local precomputed and round-trips it.
- **R CMD check**: status 0, 0 errors, 0 warnings, 3 pre-existing
  NOTEs (`.github` dir, unused `dplyr`, `split_halves.Rd` `\itemize`
  formatting). Tests: 22 / 22 pass (network test skipped).
  `pkgdown::build_site()` clean — 5 vignettes + 3 new reference pages.

## Done (2026-05-03, late evening)

- **URL bug diagnosed and fixed.** The first short URL paired a
  daf-apis state ID with `ng.banc.community/view/` (which redirects
  to GitHub-committed states, not dynamic POSTs). Refactored
  `bancr::banc_lm_scene` to use `bancr::banc_shorturl` (same as
  `bancsee()`) → spelunker URL form. Documented the
  `ng.banc.community/view/<state-name>` GitHub-PR publishing path
  separately.
- **The deeper layer-alignment bug**: my CapaR layer was in
  JRC2018U_HR voxel space (`1210 × 566 × 174 @ 519 nm`) while BANC EM
  expects layers in BANC voxel space (`2400 × 924 × 789 @ 400 nm`).
  Vignettes 2–3 (color-MIP) target JRC2018U_HR; vignette 4 (NG layer
  overlay) needs BANC voxel space — different downstream targets.
- **Elastix 5.3.1** downloaded from official GitHub release; BANC
  Elastix transform parameters fetched from
  `the-BANC-fly-connectome/fanc/transforms/transform_parameters/brain_240721/`.
  Note the naming inversion: `BANC_to_template.txt` is the file with
  FIXED = BANC and is the one used to land output on the BANC grid.
- **Vignette 4 reworked into 4 stages**: source → JRC2018F (CMTK +
  H5 points-mode), JRC2018F → BANC voxel space (transformix), BANC
  NRRD → precomputed, scene assembly. Working pipeline run on the
  CapaR no1_02 stack:
  1. `nat::xform_brain` (IS2 → JRC2018F via FCWB)
  2. `transformix -tp BANC_to_template.txt` (JRC2018F → BANC voxel)
  3. `nrrd_to_precomputed` (400 nm BANC grid, 64 chunk size)
  4. `gsutil cp -r` to
     `gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/kondo_et_al_2020/CapaR_no1_02_aligned240721_to_BANC.ng/CapaR_BANC_pc/`
  5. `bancr::banc_lm_scene` → spelunker short URL.
- **Working URL pinned**:
  [5028046288453632](https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/5028046288453632)
  — same naming convention + grid as the public
  `JRC2018F_aligned240721_to_BANC.ng` atlas.

## Done (2026-05-03, end of day)

- **Saalfeldlab transforms downloaded** (all 8, ~10 GB total, via the
  built-in `nat.jrcbrains::download_saalfeldlab_registrations()`
  after the figshare rate limit eased).
- **OpenJDK 25 installed** via brew + `JAVA_HOME` persisted in
  `~/.Rprofile`; `nat.h5reg::dr_h5reg()` reports OK.
- **Stage 1 — IS2 → JRC2018U_HR — done for real on the Kondo CapaR
  no1_02 NRRD.** The H5-based image-mode `xformimage` path is not
  yet supported in `nat.h5reg`, so the working pipeline is the
  **points-mode**: threshold + denoise (using the same defaults
  `nrrd_to_mip` would apply) → take above-threshold voxel centres as
  3-D points in IS2 microns → `xform_brain(points, sample = "IS2",
  reference = "JRC2018U")` (CMTK then H5) → voxelise into the
  `JRC2018U_HR` grid. ~448k IS2 fg voxels → ~308k JRC2018U_HR fg
  voxels (with 121 untransformable points; ~0.03 %). Total wall
  time ~3 min.
- **`inst/images/lm_capar_colormip.png`** regenerated from the
  bridged JRC2018U_HR mask — distinct bilateral cell clusters, SEZ
  signal, posterior magenta blob, all consistent with CapaR
  expression.
- **Precomputed conversion + upload**: `/tmp/CapaR_no1_02_pc/` (4.1
  MB) → `gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/kondo_et_al_2020/CapaR_no1_02_JRC2018U_HR/CapaR_no1_02_pc/`.
  1003 chunks at 519 × 519 × 1000 nm, info JSON valid, layer is
  publicly readable.
- **Working spelunker URL** for the layer overlaid on the canonical
  public BANC scene
  ([5028046288453632](https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/5028046288453632))
  built via `bancr::banc_lm_scene(shorten = TRUE)` — the same URL
  form `bancr::bancsee()` produces. Pinned into the LM-NG vignette.
  Note: the public `ng.banc.community/view/` viewer doesn't accept
  ad-hoc dynamically-POSTed states; it loads named states from
  GitHub-committed JSON files under
  `the-BANC-fly-connectome/neuroglancer_states/view/`. The vignette
  documents that publishing path separately.
- **Vignette 3 (`colormip_lm.Rmd`)** rewritten to show the
  points-mode bridging pipeline as the default and document the
  `nat.h5reg` image-mode caveat. Vignette 4 (`lm_layer_neuroglancer.Rmd`)
  pinned the working URL.

## Done (2026-05-03, evening)

- **`threshold` + `denoise` options on `nrrd_to_mip()`.** Default
  `"auto"` applies Triangle thresholding (Zack et al. 1977) + 3×3×3
  median filter on grayscale input; skips both for binary input.
  Power-user knobs: `"otsu"`, `"none"`, raw cutoff, or quantile in
  `[0, 1]`. `mmand` added to Suggests for the median filter. CapaR
  preview image regenerated — receptor cell bodies + neuropil now
  clearly visible (FG fraction 86.9% → 12.2%).
- **`bancr::banc_lm_scene(viewer = ...)`** now defaults to
  `"ng.banc.community/view"` (public BANC viewer); also accepts
  `"ng.banc.community"` (private CAVE-authenticated), `"spelunker"`,
  or any custom URL ending in `/`.
- **Vignettes reorganised into 4** per user direction:
  1. `asta_sez` — AstA search in SEZ (existing).
  2. `colormip_direct_vs_fiji` — colour-MIP from a connectome neuron
     (BANC AstA1) in JRC2018U_HR. Explicit IS2 → JRC2018F → JRC2018U
     bridging chain (now that saalfeldlab transforms are arriving;
     the old JRC2018F-bbox hack is gone). Includes back-end
     comparison section.
  3. `colormip_lm` — *(NEW)* colour-MIP from registered LM data
     (Kondo 2020 CapaR) in JRC2018U_HR. Documents the G-Node data
     source (https://doi.gin.g-node.org/10.12751/g-node.10246f/),
     the Kondo et al. 2020 citation
     (https://pubmed.ncbi.nlm.nih.gov/31914394/), and notes that the
     example NRRD path is a private Dropbox cache — readers download
     their own from G-Node. Includes a parallel back-end comparison
     section.
  4. `lm_layer_neuroglancer` — trimmed to focus on `.ng` write +
     scene assembly. Cross-links to `colormip_lm` for upstream
     bridging instead of duplicating it. Updated viewer docs to
     cover all `viewer =` options.
- **`_pkgdown.yml`** lists the four in the user's order.

## Done (2026-05-03)

- **Refactor: `banc_lm_scene()` moved into `bancr`** (per user
  direction: this repo is general; BANC-specific helpers live in
  `bancr`). Now exported from
  `flyconnectome/bancr/R/banc-lm-scene.R`. The neuronbridger
  vignette + README + reference index now call `bancr::banc_lm_scene`;
  the local copy + Rd + `_pkgdown.yml` entry are deleted. Refactored
  function adds `viewer = c("spelunker", "ng.banc.community")` so
  callers can pick the BANC-community URL prefix.
- **Switched LM example to `IS2_CapaR_no1_02_warp_m0g40c4e1e-1x16r3.nrrd`**
  in vignette + script + preview images
  (`inst/images/lm_capar_grayscale.png`,
  `inst/images/lm_capar_colormip.png` regenerated from no1_02).
- **Documented the upload + viewer policy** in the vignette + README:
  canonical bucket is
  `gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/kondo_et_al_2020/`
  but **not public-write** — readers need their own GCS / S3 / static
  HTTP host. View at `ng.banc.community/view` (default) or
  `spelunker.cave-explorer.org`.
- **CMTK 3.4.0 installed** at
  `/opt/homebrew/Cellar/cmtk-contrib-3.4.0/`; `~/.Rprofile` sets
  `options(cmtk.bindir)` + `PATH` so `nat::cmtk.version()` returns
  `‘3.4.0’` from any R session. `reformatx`, `gregxform`, etc.
  available.
- **R CMD check** (neuronbridger) post-refactor: status 0, 0 errors,
  0 warnings, 3 pre-existing NOTEs.

## In progress / pending

### Finish `nat.jrcbrains` saalfeldlab download (blocked on figshare rate-limiting)

The targeted file Stage 1 needs is `JRC2018F_JFRC2010.h5` (~1.65 GB).
Today figshare is heavily rate-limiting:

| Attempt | Method | Rate | Result |
|---|---|---|---|
| 1 | `nat.jrcbrains::download_saalfeldlab_registrations()` (single connection, default) | ~5 MB/min | Killed at 226 MB / file 1 of 8 |
| 2 | `aria2c --split=8` parallel | ~50 MB/sec | Got 99.996% (1 648 355 301 / 1 648 355 301 bytes) — actually complete; aria2's expected-size from header was wrong by 69 KB and it kept retrying |
| 3 | `aria2c --continue` resume | hit 502 | figshare started rejecting parallel connections |
| 4 | `curl -L` single connection | ~7 MB/min | Stalled |
| 5 | `aria2c --split=16 --file-allocation=none` | got 2 MB then stalled | Killed |

Net: corrupted partial file. **Cleared the partial.** User to re-run
when figshare cooperates:

```r
nat.jrcbrains::download_saalfeldlab_registrations(filenames = "JRC2018F_JFRC2010.h5")
```

or with aria2 fast path (when figshare is not rate-limiting):

```bash
cd ~/Library/Application\ Support/R/nat.jrcbrains/
aria2c --max-connection-per-server=8 --split=8 --max-tries=10 --retry-wait=10 \
       --out=JRC2018F_JFRC2010.h5 \
       'https://ndownloader.figshare.com/files/14368358?private_link=b29e25b6e47ccf9187a8'
```

Once the .h5 file is in place,
`nat.jrcbrains::register_saalfeldlab_registrations()` registers the
JFRC2 ↔ JRC2018F bridge and Stage 1 of the LM-on-NG vignette becomes
runnable.

### Stage 1 — IS2 → JRC2018F image bridge — pending the JFRC2010 download

Once `JRC2018F_JFRC2010.h5` is in place, run:

```r
suppressMessages({
  library(nat); library(nat.flybrains); library(nat.templatebrains)
  library(nat.jrcbrains)
})
nat.jrcbrains::register_saalfeldlab_registrations()

NRRD_IN  <- "/Users/asbates/Library/CloudStorage/Dropbox-HMS/Alexander Bates/neuroanat/kondo_et_al_2020/nrrd/IS2_CapaR_no1_02_warp_m0g40c4e1e-1x16r3.nrrd"
NRRD_OUT <- "~/CapaR_in_JRC2018F.nrrd"
v_is2 <- nat::read.nrrd(NRRD_IN)
v_jrc <- nat.templatebrains::xform_brain(v_is2, sample = "IS2",
                                         reference = "JRC2018F")
nat::write.nrrd(v_jrc, NRRD_OUT)
```

This goes IS2 → JFRC2 (CMTK) → JRC2018F (saalfeldlab H5). The
output volume sits at native JRC2018F dimensions (1652 × 768 × 479,
380 nm isotropic) — pass `resolution = c(380, 380, 380)` to
`nrrd_to_precomputed()`.

### Stage 4 — upload + working ng.banc.community URL — pending Stage 1

Once Stage 1 produces `CapaR_in_JRC2018F.nrrd`:

```r
neuronbridger::nrrd_to_precomputed(
  NRRD_OUT,
  output     = "/tmp/CapaR_no1_02_pc",
  resolution = c(380, 380, 380),
  data_type  = "uint8"
)
system(paste("gsutil -m cp -r /tmp/CapaR_no1_02_pc",
             "gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/kondo_et_al_2020/CapaR_no1_02/"))
u <- bancr::banc_lm_scene(
  "gs://lee-lab_brain-and-nerve-cord-fly-connectome/light_level/kondo_et_al_2020/CapaR_no1_02",
  layer_name = "Kondo 2020 — CapaR (no1_02)",
  viewer     = "ng.banc.community",
  shorten    = TRUE
)
# pin u into vignettes/lm_layer_neuroglancer.Rmd
```

### Tracking
- `inst/scripts/lm_capar_to_precomputed.R` runs end-to-end on the
  user's Dropbox path and produces a 1.5 MB precomputed dir at
  `/tmp/CapaR_IS2_pc/`. Round-trip verified.
- `bancr` lives at `/Users/asbates/Projects/flyconnectome/bancr`,
  branch `main`. Working tree dirty (new
  `R/banc-lm-scene.R` + `man/banc_lm_scene.Rd` + `NAMESPACE` export).
- CMTK on PATH via `~/.Rprofile`: `nat::cmtk.version()` → `‘3.4.0’`.
- jrcbrains transform folder cleared
  (`~/Library/Application Support/R/nat.jrcbrains/`); will re-fill
  with `download_saalfeldlab_registrations(filenames = "JRC2018F_JFRC2010.h5")`
  once figshare cooperates.
