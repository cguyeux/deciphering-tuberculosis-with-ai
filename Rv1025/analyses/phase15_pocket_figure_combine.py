#!/usr/bin/env python3
"""phase15_pocket_figure_combine.py — assemble the fpocket druggability figure (P5.1).

Combines the two PyMOL panels rendered by phase15_pocket_figure.pml into a single
A/B supplementary figure with panel letters and a colour key, at 600 DPI.

Writes: article/supplementary_materials/figureS_druggability.{png,pdf}
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent.parent
SP = Path("/tmp/claude-1000/-home-christophe-docs-codes-mtbc-Rv1025/"
          "d9a98c91-4927-4b5f-9291-2769bad55c2b/scratchpad")
OUT = ROOT / "article/supplementary_materials/figureS_druggability"

A = Image.open(SP / "panelA.png").convert("RGB")
B = Image.open(SP / "panelB.png").convert("RGB")

# normalise heights
h = min(A.height, B.height)
A = A.resize((round(A.width * h / A.height), h))
B = B.resize((round(B.width * h / B.height), h))

pad, top, keyh = 30, 70, 90
W = A.width + B.width + 3 * pad
H = h + top + keyh
canvas = Image.new("RGB", (W, H), "white")
canvas.paste(A, (pad, top))
canvas.paste(B, (2 * pad + A.width, top))

draw = ImageDraw.Draw(canvas)


def font(sz, bold=True):
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    for p in (f"/usr/share/fonts/truetype/dejavu/{name}",
              f"/usr/share/fonts/TTF/{name}"):
        if Path(p).exists():
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()


fL, fK = font(52), font(34)
draw.text((pad + 8, 8), "A", fill="black", font=fL)
draw.text((2 * pad + A.width + 8, 8), "B", fill="black", font=fL)

# colour key strip at the bottom
y = top + h + 22
x = pad + 10
orange, marine = (255, 127, 0), (60, 110, 200)
draw.ellipse([x, y, x + 34, y + 34], fill=orange)
draw.text((x + 46, y - 2), "Metal-site pocket (Pocket 3, DS 0.24)", fill="black", font=fK)
x2 = x + 1040
draw.ellipse([x2, y, x2 + 34, y + 34], fill=marine)
draw.text((x2 + 46, y - 2), "Druggable pocket (Pocket 1, DS 0.65)", fill="black", font=fK)

canvas.save(f"{OUT}.png", dpi=(600, 600))
canvas.save(f"{OUT}.pdf", "PDF", resolution=600.0)
print(f"Written: {OUT}.png / .pdf  ({canvas.width}x{canvas.height})")
