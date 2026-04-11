#!/usr/bin/env python3
"""Animate a 2D wave equation .tgm file.

Usage:
    python3 visualize_tgm.py <file> [skip] [--save] [--zoom r0:r1,c0:c1]

Requires: tensogram, numpy, matplotlib, Pillow
"""

import io
import os
import sys

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import tensogram
from PIL import Image


def decode_zoom_frame(f, idx, grid_size, zoom):
    """Decode one frame using range decode for a sub-region."""
    (r0, r1), (c0, c1) = zoom
    ranges = [(r * grid_size + c0, c1 - c0) for r in range(r0, r1)]
    return f.file_decode_range(idx, 0, ranges, join=True).reshape(r1 - r0, c1 - c0)


def parse_zoom(s):
    row_part, col_part = s.split(",")
    r0, r1 = map(int, row_part.split(":"))
    c0, c1 = map(int, col_part.split(":"))
    return (r0, r1), (c0, c1)


def main():
    tgm_path = sys.argv[1] if len(sys.argv) > 1 else "grid.tgm"
    skip = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 10
    save = "--save" in sys.argv
    zoom = None
    for i, arg in enumerate(sys.argv):
        if arg == "--zoom" and i + 1 < len(sys.argv):
            zoom = parse_zoom(sys.argv[i + 1])

    print(f"Loading {tgm_path}...")
    f = tensogram.TensogramFile.open(tgm_path)
    count = f.message_count()

    info = f.file_decode_descriptors(0)
    grid_size = info["descriptors"][0].shape[0]
    dt_val = info["metadata"]["wave"]["dt"]
    print(f"  {count} messages, {grid_size}x{grid_size}, dt={dt_val}")

    if zoom:
        (r0, r1), (c0, c1) = zoom
        print(f"  zoom: rows [{r0}:{r1}], cols [{c0}:{c1}] (range decode)")

    indices = list(range(0, count, skip))
    print(f"Decoding {len(indices)} frames (every {skip}th)...")

    if zoom:
        frames = []
        for i, idx in enumerate(indices):
            frames.append(decode_zoom_frame(f, idx, grid_size, zoom))
            if (i + 1) % 50 == 0:
                print(f"  decoded {i + 1}/{len(indices)}")
    else:
        frames = [msg.objects[0][1][::4, ::4] for msg in f[::skip]]

    vmax = max(np.abs(fr).max() for fr in frames)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=50)
    im = ax.imshow(frames[0], cmap="seismic", origin="lower", vmin=-vmax, vmax=vmax)
    plt.colorbar(im, ax=ax, shrink=0.8)
    title = ax.set_title("t = 0.000")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    def update(frame_num):
        im.set_data(frames[frame_num])
        title.set_text(f"t = {indices[frame_num] * dt_val:.3f}")
        return [im, title]

    ani = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=200, blit=True
    )

    if save:
        out = tgm_path.replace(".tgm", ".gif")
        print(f"Rendering {len(frames)} frames to GIF (PIL optimized)...")
        pil_frames = []
        for i in range(len(frames)):
            update(i)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0.3)
            buf.seek(0)
            pil_frames.append(
                Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=128)
            )
        pil_frames[0].save(
            out,
            save_all=True,
            append_images=pil_frames[1:],
            duration=33,
            loop=0,
            optimize=True,
        )
        print(f"Saved {out} ({os.path.getsize(out) / 1024 / 1024:.1f} MiB)")
    else:
        plt.show()


if __name__ == "__main__":
    main()
