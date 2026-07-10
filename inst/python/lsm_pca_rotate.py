#!/usr/bin/env python3
"""PCA-based pre-rotation of a partial-FOV LSM onto a target template's
axes, with sub-voxel placement at a user-chosen anchor.

Designed for cases where the LSM captures only a fragment of the brain
or VNC (e.g. just the abdominal tip) at an arbitrary orientation. The
NC82 mask's principal axes are computed in physical microns and mapped
to the template's (X, Y, Z) axes; the longest PC -> Y, medium -> X,
shortest -> Z. Sign disambiguation puts the GFP CoG at +Y. The result
is a JRC2018VNCF-aligned NRRD with identity Direction, ready to feed
into downstream image-space tools.

Inputs are three NRRDs from the SAME LSM at native voxdims:
  NC82_NRRD       structural reference channel
  GFP_NRRD        label channel
  MASK_NRRD       binary mask of NC82 tissue (Otsu + close + largest CC)

Plus a fixed reference for placement:
  TEMPLATE_NRRD       output-grid reference (its sizes / Direction are not
                      used here, only voxdims via the target args).
  FIXED_MASK_NRRD     binary mask in template coords; its CoG is the
                      default placement target. Override with --target-um
                      "x,y,z" to anchor the rotated source at a known
                      reference point (e.g. a previous registration's
                      GFP CoG).

Output (written to OUT_DIR):
  nc82_rot.nrrd, gfp_rot.nrrd, mask_rot.nrrd   isotropic 0.4/0.4/0.7 um
                                               grid with identity Direction
                                               and an origin set so the
                                               alignment source lands on
                                               the target.
  pca_transform.json                            R, signs, target, origin.

Usage:
  python lsm_pca_rotate.py NC82 GFP MASK TEMPLATE FIXED_MASK OUT_DIR \\
      [--align-source nc82_mask|gfp_cog]  [--target-um "x,y,z"] \\
      [--pc2-sign +1|-1] [--out-spacing 0.4] [--out-z-spacing 0.7]
"""
import os, json, argparse
import numpy as np, SimpleITK as sitk
from scipy import ndimage as ndi

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("nc82_nrrd")
ap.add_argument("gfp_nrrd")
ap.add_argument("mask_nrrd")
ap.add_argument("template_nrrd")
ap.add_argument("fixed_mask_nrrd")
ap.add_argument("out_dir")
ap.add_argument("--pc2-sign", type=int, default=+1, choices=[+1, -1])
ap.add_argument("--out-spacing", type=float, default=0.4)
ap.add_argument("--out-z-spacing", type=float, default=0.7)
ap.add_argument("--align-source", default="nc82_mask",
                choices=["nc82_mask", "gfp_cog"])
ap.add_argument("--target-um", default=None,
                help='explicit "x,y,z" um target in template coords; '
                     'overrides fixed_mask_nrrd CoG')
ap.add_argument("--prefix", default="",
                help="optional output filename prefix")
args = ap.parse_args()

os.makedirs(args.out_dir, exist_ok=True)

def read(p):
    im = sitk.ReadImage(p)
    arr = sitk.GetArrayFromImage(im)
    return arr, np.asarray(im.GetSpacing()), np.asarray(im.GetOrigin()), im

nc82, sN, _, _ = read(args.nc82_nrrd)
gfp,  sG, _, _ = read(args.gfp_nrrd)
mask, sM, _, _ = read(args.mask_nrrd)
print(f"NC82 shape {nc82.shape} spacing {sN}")
print(f"GFP  shape {gfp.shape}  spacing {sG}")
print(f"mask shape {mask.shape} spacing {sM}")

if args.target_um:
    fcog_world = np.array([float(x) for x in args.target_um.split(",")])
    print(f"target (override) (x,y,z um): {fcog_world}")
else:
    fix, sF, oF, _ = read(args.fixed_mask_nrrd)
    fix_idx = np.argwhere(fix.astype(bool))
    fc = fix_idx.mean(axis=0)
    fcog_world = np.array([fc[2]*sF[0], fc[1]*sF[1], fc[0]*sF[2]]) + oF
    print(f"fixed mask CoG (x,y,z um): {fcog_world}")

# --- PCA on the NC82 mask ---
m_idx = np.argwhere(mask.astype(bool))
pts = np.stack([m_idx[:,2]*sM[0], m_idx[:,1]*sM[1], m_idx[:,0]*sM[2]], axis=1)
mu = pts.mean(axis=0)
X  = pts - mu
print(f"moving mask CoG (x,y,z um): {mu}  (N {len(pts):,})")
cov = (X.T @ X) / len(X)
w, V = np.linalg.eigh(cov)
order = np.argsort(-w); w = w[order]; V = V[:, order]
print(f"eigenvalues: {w}")

# Sign disambiguation: PC1 toward GFP (template posterior).
g_idx = np.argwhere(gfp > 0)
if len(g_idx) == 0:
    raise SystemExit("GFP channel is empty -- nothing to align.")
g_world = np.stack([g_idx[:,2]*sG[0], g_idx[:,1]*sG[1],
                    g_idx[:,0]*sG[2]], axis=1).mean(axis=0)
proj1 = (g_world - mu) @ V[:, 0]
sign1 = +1 if proj1 > 0 else -1
V[:, 0] *= sign1
V[:, 1] *= args.pc2_sign
print(f"GFP offset on PC1: {proj1:.3f} -> v1 sign={sign1:+d}")

# Map PC1->Y, PC2->X, PC3->Z. Enforce det(R)=+1 by flipping PC3 if needed.
M = np.column_stack([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
R = M @ V.T
if np.linalg.det(R) < 0:
    V[:, 2] *= -1
    R = M @ V.T
    print("flipped PC3 sign to enforce det(R)=+1")
print(f"final R:\n{R}\ndet(R): {np.linalg.det(R):+.4f}")

# --- Output grid: axis-aligned bbox of rotated FOV ---
out_sp = np.array([args.out_spacing, args.out_spacing, args.out_z_spacing])
corners = np.array([(z, y, x) for z in [0, mask.shape[0]-1]
                              for y in [0, mask.shape[1]-1]
                              for x in [0, mask.shape[2]-1]])
cw = np.stack([corners[:,2]*sM[0], corners[:,1]*sM[1],
               corners[:,0]*sM[2]], axis=1)
rc = (R @ (cw - mu).T).T
bbox_lo, bbox_hi = rc.min(axis=0), rc.max(axis=0)
out_size = np.ceil((bbox_hi - bbox_lo) / out_sp).astype(int)
print(f"rotated bbox: {bbox_hi - bbox_lo}  out_size: {out_size}")

# --- Resample via scipy.ndimage.affine_transform ---
# For output voxel idx_zyx, find source idx_zyx via:
#   src_world_xyz = R^T @ (origin_out + diag(out_sp) @ P @ idx_zyx) + mu
#   src_idx_xyz   = diag(1/src_sp) @ src_world_xyz
#   src_idx_zyx   = P @ src_idx_xyz
# where P permutes zyx<->xyz.
P = np.array([[0,0,1],[0,1,0],[1,0,0]], float)
def resample(arr, src_sp):
    Dout = np.diag(out_sp); Din = np.diag(src_sp); Dinv = np.linalg.inv(Din)
    A = P @ Dinv @ R.T @ Dout @ P
    b = P @ (Dinv @ R.T @ bbox_lo + Dinv @ mu)
    return ndi.affine_transform(arr, A, offset=b,
                                output_shape=(out_size[2], out_size[1], out_size[0]),
                                order=1, mode="constant", cval=0)

nc82_rot = resample(nc82, sN)
gfp_rot  = resample(gfp,  sG)
mask_rot = resample(mask.astype(np.float32), sM) >= 0.5

# --- Place final origin so alignment-source lands on target ---
mr_idx = np.argwhere(mask_rot)
mr_world = np.stack([mr_idx[:,2]*out_sp[0], mr_idx[:,1]*out_sp[1],
                     mr_idx[:,0]*out_sp[2]], axis=1).mean(axis=0)
gw = gfp_rot.astype(np.float32)
gi = np.argwhere(gw > 0)
if len(gi):
    gw_world = np.average(np.stack([gi[:,2]*out_sp[0], gi[:,1]*out_sp[1],
                                    gi[:,0]*out_sp[2]], axis=1),
                          axis=0, weights=gw[gw > 0])
else:
    gw_world = None
source_in_rot = mr_world if args.align_source == "nc82_mask" else gw_world
if source_in_rot is None:
    raise SystemExit("--align-source=gfp_cog but no GFP signal post-rotation.")
final_origin = fcog_world - source_in_rot
print(f"aligning {args.align_source} -> target {fcog_world}")
print(f"final origin (x,y,z um): {final_origin}")

def write(arr, path):
    im = sitk.GetImageFromArray(arr.astype(np.uint8))
    im.SetSpacing(tuple(out_sp.astype(float)))
    im.SetOrigin(tuple(final_origin.astype(float)))
    sitk.WriteImage(im, path, useCompression=True)
    print("wrote", path)

px = args.prefix
write(nc82_rot, os.path.join(args.out_dir, f"{px}nc82_rot.nrrd"))
write(gfp_rot,  os.path.join(args.out_dir, f"{px}gfp_rot.nrrd"))
write(mask_rot.astype(np.uint8),
      os.path.join(args.out_dir, f"{px}nc82_mask_rot.nrrd"))

with open(os.path.join(args.out_dir, f"{px}pca_transform.json"), "w") as f:
    json.dump({
        "R": R.tolist(),
        "mu_moving_um": mu.tolist(),
        "target_um": fcog_world.tolist(),
        "out_spacing": out_sp.tolist(),
        "final_origin_um": final_origin.tolist(),
        "v1_sign": int(sign1),
        "v2_sign": int(args.pc2_sign),
    }, f, indent=2)
print("dumped pca_transform.json")
