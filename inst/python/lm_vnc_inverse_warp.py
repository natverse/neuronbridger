"""Inverse-warp a JRC2018VNCF-aligned LM volume into BANC voxel space using
the BANC -> JRC2018VNCF tpsreg coefficients (exported by R as JSON).

For each output BANC voxel center, applies the BANC->JRC2018VNCF TPS to find
its source coord in JRC2018VNCF microns, then samples the LM volume at that
coord using max-pool over the output voxel's physical footprint. No
threshold; every LM voxel's intensity is preserved (the same way transformix
preserves intensity for the brain pipeline).

Usage:
  python lm_vnc_inverse_warp.py LM_NRRD TPS_JSON OUT_NRRD [DS_NM_X DS_NM_Y DS_NM_Z]

  LM_NRRD     JRC2018VNCF-aligned LM volume (the GFP channel after Stage A')
  TPS_JSON    bancr::banc_to_jrcvnc2018f_tpsreg exported as JSON; see the R
              helper export_banc_to_jrcvnc2018f_tpsreg() in lm_pipeline.R
  OUT_NRRD    output BANC-voxel-grid NRRD
  DS_NM_*     output voxel size in nm (default 1600 1600 800).

The output spatial extent is the bbox of the BANC VNC neuropil mesh + 10%
padding, so the work scales with cord volume rather than full BANC volume.
"""
import os, sys, json, numpy as np, SimpleITK as sitk

if len(sys.argv) < 4:
    sys.exit(__doc__)

LM_NRRD  = sys.argv[1]
TPS_JSON = sys.argv[2]
OUT_NRRD = sys.argv[3]
DS_NM    = (int(sys.argv[4]) if len(sys.argv) > 4 else 1600,
            int(sys.argv[5]) if len(sys.argv) > 5 else 1600,
            int(sys.argv[6]) if len(sys.argv) > 6 else 800)

BANC_DIM = (2622, 2950, 789)   # at 400 nm
OUT_DIM = (BANC_DIM[0] * 400 // DS_NM[0],
           BANC_DIM[1] * 400 // DS_NM[1],
           BANC_DIM[2] * 400 // DS_NM[2])
print(f"output grid: {OUT_DIM} @ {DS_NM} nm")

# Read LM volume (JRC2018VNCF)
img = sitk.ReadImage(LM_NRRD)
arr_lm = sitk.GetArrayFromImage(img)             # (z, y, x) float32
arr_lm = np.clip(arr_lm, 0, 255).astype(np.float32)
sx_lm, sy_lm, sz_lm = img.GetSpacing()           # um per voxel
print(f"LM (JRC2018VNCF): shape {arr_lm.shape}  spacing {(sx_lm, sy_lm, sz_lm)} um")

# Load TPS landmarks (BANC nm refmat -> JRC2018VNCF um tarmat)
tps = json.load(open(TPS_JSON))
refmat = np.asarray(tps["refmat"], dtype=np.float64)
tarmat = np.asarray(tps["tarmat"], dtype=np.float64)
n = refmat.shape[0]
print(f"tpsreg landmarks: {n}")

# Solve TPS once.
# K_ij = -|src_i - src_j|; L = [[K, P], [P^T, 0]]; coef = L^-1 [tar; 0]
print("solving TPS...")
diff = refmat[:, None, :] - refmat[None, :, :]
K = -np.linalg.norm(diff, axis=-1)
P = np.hstack([np.ones((n, 1)), refmat])
L = np.zeros((n + 4, n + 4))
L[:n, :n] = K
L[:n, n:] = P
L[n:, :n] = P.T
Y = np.vstack([tarmat, np.zeros((4, 3))])
coef = np.linalg.solve(L, Y)
print("TPS solved")

# Restrict iteration to the BANC VNC neuropil bbox + 10% padding. The
# bbox below is hard-coded from `bancr::banc_vnc_neuropil.surf`'s vertex
# extent; widening it keeps any nearby off-cord signal visible.
mb = {"x0": 300_000, "x1":   700_000,
      "y0": 480_000, "y1": 1_080_000,
      "z0":  30_000, "z1":   330_000}
ix0, ix1 = mb["x0"] // DS_NM[0], mb["x1"] // DS_NM[0] + 1
iy0, iy1 = mb["y0"] // DS_NM[1], mb["y1"] // DS_NM[1] + 1
iz0, iz1 = mb["z0"] // DS_NM[2], mb["z1"] // DS_NM[2] + 1
ix1, iy1, iz1 = (min(ix1, OUT_DIM[0]),
                 min(iy1, OUT_DIM[1]),
                 min(iz1, OUT_DIM[2]))
sub_dim = (ix1 - ix0, iy1 - iy0, iz1 - iz0)
print(f"bbox sub-grid: {sub_dim}  ({sub_dim[0]*sub_dim[1]*sub_dim[2]:,} voxels)")

out_arr = np.zeros((OUT_DIM[2], OUT_DIM[1], OUT_DIM[0]), dtype=np.uint8)

# Output voxel footprint in source-voxel space (for max-pool sampling)
half_x = (DS_NM[0] / 1000.0) / sx_lm / 2.0
half_y = (DS_NM[1] / 1000.0) / sy_lm / 2.0
half_z = (DS_NM[2] / 1000.0) / sz_lm / 2.0

CHUNK = 200_000
for kk in range(iz0, iz1):
    bz_nm = kk * DS_NM[2]
    ii = np.arange(ix0, ix1)
    jj = np.arange(iy0, iy1)
    II, JJ = np.meshgrid(ii, jj, indexing="ij")
    bx_nm = II.ravel() * DS_NM[0]
    by_nm = JJ.ravel() * DS_NM[1]
    q = np.column_stack([bx_nm, by_nm,
                         np.full(bx_nm.shape, bz_nm)]).astype(np.float64)

    out_um = np.empty_like(q)
    for s in range(0, q.shape[0], CHUNK):
        e = min(s + CHUNK, q.shape[0])
        d = q[s:e, None, :] - refmat[None, :, :]
        Kc = -np.linalg.norm(d, axis=-1)
        Pq = np.hstack([np.ones((e - s, 1)), q[s:e]])
        out_um[s:e] = Kc @ coef[:n] + Pq @ coef[n:]

    vx = out_um[:, 0] / sx_lm
    vy = out_um[:, 1] / sy_lm
    vz = out_um[:, 2] / sz_lm

    inb = ((vx >= 0) & (vx < arr_lm.shape[2] - 1) &
           (vy >= 0) & (vy < arr_lm.shape[1] - 1) &
           (vz >= 0) & (vz < arr_lm.shape[0] - 1))

    vals = np.zeros(q.shape[0], dtype=np.float32)
    if inb.any():
        # max-pool over each output voxel's physical footprint
        xs0 = np.maximum(0, np.floor(vx[inb] - half_x)).astype(np.int32)
        xs1 = np.minimum(arr_lm.shape[2] - 1,
                          np.ceil(vx[inb] + half_x)).astype(np.int32)
        ys0 = np.maximum(0, np.floor(vy[inb] - half_y)).astype(np.int32)
        ys1 = np.minimum(arr_lm.shape[1] - 1,
                          np.ceil(vy[inb] + half_y)).astype(np.int32)
        zs0 = np.maximum(0, np.floor(vz[inb] - half_z)).astype(np.int32)
        zs1 = np.minimum(arr_lm.shape[0] - 1,
                          np.ceil(vz[inb] + half_z)).astype(np.int32)
        tmp = np.zeros(int(inb.sum()), dtype=np.float32)
        for kk2 in range(len(tmp)):
            tmp[kk2] = arr_lm[zs0[kk2]:zs1[kk2] + 1,
                              ys0[kk2]:ys1[kk2] + 1,
                              xs0[kk2]:xs1[kk2] + 1].max()
        vals[inb] = tmp
    vals = np.clip(np.round(vals), 0, 255).astype(np.uint8)

    out_slice = vals.reshape(sub_dim[0], sub_dim[1])
    out_arr[kk, iy0:iy1, ix0:ix1] = out_slice.T   # (Y, X)

    if (kk - iz0) % 20 == 0:
        nz_count = (out_arr[kk] > 0).sum()
        print(f"  z {kk-iz0}/{sub_dim[2]}  filled-this-slice={nz_count}")

print(f"total non-zero voxels: {(out_arr > 0).sum():,}  "
      f"max={out_arr.max()}  mean={out_arr.mean():.2f}")

img_out = sitk.GetImageFromArray(out_arr)
img_out.SetSpacing([v / 1000.0 for v in DS_NM])
sitk.WriteImage(img_out, OUT_NRRD, useCompression=True)
print("wrote", OUT_NRRD)
