#!/usr/bin/env python3
"""Build a neuropil mask from an NC82 volume and apply it to the
matching label channel.

Used when only the part of the label signal that sits inside the
neuropil should contribute to downstream matching (e.g. NBLAST, colour-
MIP search) -- background, cuticle bleed-through and signal that falls
outside the cord neuropil are zeroed.

The mask is built from the NC82 channel by Otsu thresholding followed
by a binary closing (3x3x3, 2 iterations) and a small dilation (2 vox)
to fully cover the neuropil edge. The largest connected component is
kept so isolated background blobs are dropped.

Usage:
  python lsm_neuropil_mask.py NC82_NRRD GFP_NRRD OUT_NRRD \\
      [--mask-out NEUROPIL_MASK_NRRD]
"""
import sys, os, argparse, numpy as np, SimpleITK as sitk
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("nc82_nrrd")
ap.add_argument("gfp_nrrd")
ap.add_argument("out_nrrd")
ap.add_argument("--mask-out", default=None,
                help="optional path to also write the neuropil mask NRRD")
args = ap.parse_args()

nc = sitk.ReadImage(args.nc82_nrrd)
gf = sitk.ReadImage(args.gfp_nrrd)
nc_arr = sitk.GetArrayFromImage(nc).astype(np.float32)
gf_arr = sitk.GetArrayFromImage(gf).astype(np.float32)
if nc_arr.shape != gf_arr.shape:
    raise SystemExit(f"shape mismatch: NC82 {nc_arr.shape} vs GFP {gf_arr.shape}")
print(f"NC82 shape {nc_arr.shape}  nonzero {(nc_arr>0).sum():,}")
print(f"GFP  shape {gf_arr.shape}  nonzero {(gf_arr>0).sum():,}")

nz = nc_arr[nc_arr > 0]
if nz.size == 0:
    raise SystemExit("NC82 channel is empty -- cannot build neuropil mask.")
t = threshold_otsu(nz)
print(f"NC82 Otsu threshold: {t:.2f}")

m = nc_arr > t
m = ndi.binary_closing(m, iterations=2)
lab, n = ndi.label(m)
if n > 0:
    sizes = ndi.sum(m, lab, range(1, n + 1))
    keep = int(np.argmax(sizes)) + 1
    m = (lab == keep)
m = ndi.binary_dilation(m, iterations=2)
print(f"neuropil mask voxels: {int(m.sum()):,}")

gf_out = np.where(m, gf_arr, 0).astype(np.uint8)
print(f"GFP after neuropil mask: nonzero {(gf_out>0).sum():,}  "
      f"({100*(gf_out>0).sum()/max(1,(gf_arr>0).sum()):.1f}% of pre-mask)")

out = sitk.GetImageFromArray(gf_out)
out.CopyInformation(gf)
sitk.WriteImage(out, args.out_nrrd, useCompression=True)
print("wrote", args.out_nrrd)

if args.mask_out:
    mout = sitk.GetImageFromArray(m.astype(np.uint8))
    mout.CopyInformation(nc)
    sitk.WriteImage(mout, args.mask_out, useCompression=True)
    print("wrote", args.mask_out)
