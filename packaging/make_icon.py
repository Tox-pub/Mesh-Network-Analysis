# -*- coding: utf-8 -*-
"""
make_icon.py - draw the application icon.

The icon is a small concept network: a highlighted hub with satellites around
it, which is what the pipeline actually produces.

Each size is drawn natively rather than rendered once and downsampled. At 16 px
a scaled-down 256 px drawing turns into grey mush - strokes fall below one pixel
and the nodes merge. Drawing per size lets the stroke weights and the node radii
stay legible where it matters, which is the taskbar.

Run:  python packaging/make_icon.py
Writes src/mesh_workbench/assets/mesh_workbench.ico (and PNG previews).
"""

import math
import os
import sys

try:
    from PIL import Image, ImageDraw
except ImportError:
    sys.exit('Pillow is required: python -m pip install pillow')

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(HERE), 'src', 'mesh_workbench', 'assets')
SIZES = (16, 24, 32, 48, 64, 128, 256)

BG = (27, 42, 74)            # deep navy plate
EDGE = (127, 168, 208)       # steel blue links
NODE = (236, 242, 250)       # near-white satellites
HUB = (240, 168, 48)         # amber hub - the one thing the eye lands on

# Satellites on a circle. Five reads as a network without crowding at 16 px;
# more and they merge into a ring.
SATELLITES = 5
START_ANGLE = -90            # one node straight up, so the shape looks upright


def draw(size):
    # Supersample, but only by 4x and only for the geometry - the stroke widths
    # below are computed from the final size, so they survive the reduction.
    ss = 4 if size >= 24 else 8
    px = size * ss
    img = Image.new('RGBA', (px, px), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)

    radius = px * 0.16
    d.rounded_rectangle([0, 0, px - 1, px - 1], radius=radius, fill=BG + (255,))

    cx = cy = px / 2
    ring = px * 0.30
    hub_r = max(px * 0.105, ss * 1.5)
    node_r = max(px * 0.072, ss * 1.0)
    edge_w = max(int(px * 0.030), ss)

    points = []
    for i in range(SATELLITES):
        a = math.radians(START_ANGLE + i * 360 / SATELLITES)
        points.append((cx + ring * math.cos(a), cy + ring * math.sin(a)))

    # Spokes first, then one chord, so the graph reads as connected rather than
    # as a star. The chord is skipped at the smallest size where it just blurs.
    for x, y in points:
        d.line([cx, cy, x, y], fill=EDGE + (255,), width=edge_w)
    if size >= 24:
        d.line([points[1], points[2]], fill=EDGE + (255,), width=edge_w)
        d.line([points[3], points[4]], fill=EDGE + (255,), width=edge_w)

    for x, y in points:
        d.ellipse([x - node_r, y - node_r, x + node_r, y + node_r],
                  fill=NODE + (255,))
    d.ellipse([cx - hub_r, cy - hub_r, cx + hub_r, cy + hub_r], fill=HUB + (255,))

    return img.resize((size, size), Image.LANCZOS)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    frames = [draw(s) for s in SIZES]
    ico = os.path.join(OUT_DIR, 'mesh_workbench.ico')
    # Pillow writes every supplied size into the .ico when they are passed as
    # append_images; a single image would leave Windows to rescale.
    frames[-1].save(ico, format='ICO',
                    sizes=[(s, s) for s in SIZES],
                    append_images=frames[:-1])
    png = os.path.join(OUT_DIR, 'mesh_workbench.png')
    frames[-1].save(png, format='PNG')          # for Linux .desktop entries
    print(f'  wrote {ico}  ({os.path.getsize(ico):,} bytes, {len(SIZES)} sizes)')
    print(f'  wrote {png}  ({os.path.getsize(png):,} bytes)')

    preview = os.path.join(HERE, '_icon_preview.png')
    strip = Image.new('RGBA', (sum(s for s in SIZES) + 8 * len(SIZES), 256),
                      (255, 255, 255, 255))
    x = 0
    for s, f in zip(SIZES, frames):
        strip.paste(f, (x, (256 - s) // 2), f)
        x += s + 8
    strip.save(preview)
    print(f'  wrote {preview} (visual check, not shipped)')


if __name__ == '__main__':
    main()
