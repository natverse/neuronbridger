#!/usr/bin/env python3
"""Extract NC82 (channel 1) and GFP (channel 0) from a Carl Zeiss
LSM file to two NRRDs with proper voxdims metadata.

Used by the Deng et al. 2019 LM -> JRC2018 -> BANC pipeline (see
neuronbridger::lsm_to_banc_layer). Stand-alone so it can be called
from R via system2() without rJava / scyjava plumbing.

Usage:
    python3 lsm_to_nrrd.py <input.lsm> <output_dir> [--prefix STEM]

Writes <output_dir>/<stem>_nc82.nrrd and <stem>_gfp.nrrd, both uint8,
gzip-compressed, voxdims pulled from LSM metadata (XYZ in microns).

Channel convention assumed (verified for Deng et al. 2019 stacks):
    ch 0 = GFP signal (sparse)
    ch 1 = NC82 anti-Brp reference (broad neuropil)

If a future LSM cohort swaps the channel order, this script will
silently swap labels --- inspect the resulting NRRDs and re-extract
with --gfp-channel / --nc82-channel overrides.
"""
import argparse
import os
import sys

import SimpleITK as sitk
import tifffile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output_dir")
    ap.add_argument("--prefix", default=None,
                    help="output stem; default = input file basename")
    ap.add_argument("--gfp-channel", type=int, default=0)
    ap.add_argument("--nc82-channel", type=int, default=1)
    args = ap.parse_args()

    src = args.input
    out_dir = args.output_dir
    stem = args.prefix or os.path.splitext(os.path.basename(src))[0]
    os.makedirs(out_dir, exist_ok=True)

    with tifffile.TiffFile(src) as t:
        if not t.is_lsm:
            print(f"WARNING: {src} is not an LSM-format TIFF", file=sys.stderr)
        arr = t.series[0].asarray()  # (Z, C, Y, X) for typical LSM
        meta = t.lsm_metadata or {}

    if arr.ndim != 4 or arr.shape[1] < 2:
        sys.exit(f"unexpected LSM shape {arr.shape}; expected (Z, C>=2, Y, X)")

    vx_um = float(meta.get("VoxelSizeX", 1e-6) * 1e6)
    vy_um = float(meta.get("VoxelSizeY", 1e-6) * 1e6)
    vz_um = float(meta.get("VoxelSizeZ", 1e-6) * 1e6)

    def write_channel(ch_idx, label):
        if ch_idx >= arr.shape[1]:
            sys.exit(f"channel index {ch_idx} out of range "
                     f"(LSM has {arr.shape[1]} channels)")
        # SimpleITK expects (Z, Y, X) for 3D; the LSM array already has that
        # order at axes 0, 2, 3, with channel at axis 1.
        vol = arr[:, ch_idx, :, :]
        img = sitk.GetImageFromArray(vol)
        img.SetSpacing([vx_um, vy_um, vz_um])
        img.SetOrigin([0.0, 0.0, 0.0])
        out = os.path.join(out_dir, f"{stem}_{label}.nrrd")
        sitk.WriteImage(img, out, useCompression=True)
        print(f"  wrote {out}  shape(Z,Y,X)={vol.shape}  "
              f"voxdims_um={vx_um:.4f},{vy_um:.4f},{vz_um:.4f}")
        return out

    write_channel(args.nc82_channel, "nc82")
    write_channel(args.gfp_channel,  "gfp")


if __name__ == "__main__":
    main()
