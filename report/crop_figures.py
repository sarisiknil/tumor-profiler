#!/usr/bin/env python3
"""Crop the dashboard screenshots to the region that carries the point, for use on slides.

A full page screenshot is 16:10; a slide's content area is roughly 2.2:1. Fitting the whole page therefore
leaves it small enough to be unreadable from the back of a room. Each crop below keeps the part of the page the
slide is actually about, and the uncropped originals remain for the report's appendix.

Fractions are (left, top, right, bottom) of the full image.
"""
from pathlib import Path

from PIL import Image

REPO = Path(__file__).resolve().parents[1]
FIG = REPO / "report" / "figures"

# the sidebar occupies the left ~19 % of every page and carries no information on a slide
CROPS = {
    "02_qc_library.png":  (0.187, 0.145, 1.0, 0.66),   # header, library facts, base-composition chart
    "03_variants.png":    (0.187, 0.215, 1.0, 0.72),   # counts, the agreement callout, the table
    "04_rna_fusions.png": (0.187, 0.265, 1.0, 0.78),   # the artefact banner and the flagged calls
    "05_pathways.png":    (0.187, 0.245, 1.0, 0.99),   # the pathway bar chart
    "06_therapies.png":   (0.187, 0.215, 1.0, 0.80),
}


def main():
    for name, (l, t, r, b) in CROPS.items():
        src = FIG / name
        if not src.exists():
            print(f"  missing {name}")
            continue
        im = Image.open(src)
        w, h = im.size
        box = (int(l * w), int(t * h), int(r * w), int(b * h))
        out = FIG / ("crop_" + name)
        im.crop(box).save(out)
        cw, ch = im.crop(box).size
        print(f"  {name:22s} -> {out.name:27s} {cw}x{ch}  aspect {cw/ch:.2f}")


if __name__ == "__main__":
    main()
