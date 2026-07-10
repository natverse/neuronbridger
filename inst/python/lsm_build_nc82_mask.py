#!/usr/bin/env python3
"""Build a tissue mask from an NC82 LSM channel via Otsu + close +
largest CC + small dilation. Used as the moving-image mask for
`lsm_pca_rotate.py` (which needs a binary mask to compute PCA axes).

Usage:
  python lsm_build_nc82_mask.py NC82_NRRD MASK_OUT_NRRD
"""
import sys, numpy as np, SimpleITK as sitk
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu

if len(sys.argv) != 3:
    sys.exit(__doc__)
NC82_IN, MASK_OUT = sys.argv[1], sys.argv[2]

im = sitk.ReadImage(NC82_IN)
arr = sitk.GetArrayFromImage(im)
t = threshold_otsu(arr[arr > 0]) if (arr > 0).any() else 0
print(f"NC82 Otsu threshold = {t} (min={arr.min()} max={arr.max()})")

bw = arr > t
bw = ndi.binary_closing(bw, iterations=2)
lab, n = ndi.label(bw)
if n:
    sizes = ndi.sum(bw, lab, range(1, n + 1))
    bw = (lab == int(np.argmax(sizes)) + 1)
bw = ndi.binary_dilation(bw, iterations=2).astype(np.uint8)
print(f"NC82 mask voxels: {int(bw.sum()):,}  ({100*bw.mean():.1f}% of FOV)")

out = sitk.GetImageFromArray(bw)
out.CopyInformation(im)
sitk.WriteImage(out, MASK_OUT, useCompression=True)
print("wrote", MASK_OUT)
