#!/usr/bin/env python3
"""Visualize a 2D wave equation .tgm file produced by the MPI solver.

Requires: tensogram (maturin develop), numpy, matplotlib, Pillow

Usage:
    python3 visualize_tgm.py <file> [skip] [--save] [--zoom row_start:row_end,col_start:col_end]
"""

import io
import os
import sys

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import tensogram
from PIL import Image


def decode_zoom(f, idx, row_range, col_range, grid_size):
    """Decode a rectangular sub-region using range decode (tensorjump).

    Maps 2D row/column ranges to 1D element offsets in the flattened row-major
    array: row i starts at element i*grid_size.
    """
    r0, r1 = row_range
    c0, c1 = col_range
    raw = f.read_message(idx)
    ranges = [(r * grid_size + c0, c1 - c0) for r in range(r0, r1)]
    flat = tensogram.decode_range(raw, 0, ranges)
    return flat.reshape(r1 - r0, c1 - c0)


def parse_zoom(zoom_str):
    row_part, col_part = zoom_str.split(",")
    r0, r1 = map(int, row_part.split(":"))
    c0, c1 = map(int, col_part.split(":"))
    return (r0, r1), (c0, c1)


def main():
    tgm_path = sys.argv[1] if len(sys.argv) > 1 else "grid.tgm"
    skip = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    zoom = None
    for i, arg in enumerate(sys.argv):
        if arg == "--zoom" and i + 1 < len(sys.argv):
            zoom = parse_zoom(sys.argv[i + 1])

    print(f"Loading {tgm_path}...")
    f = tensogram.TensogramFile.open(tgm_path)
    count = f.message_count()

    meta0, _ = f.decode_message(0)
    dt = meta0["wave"]["dt"]
    grid_size = int(meta0.objects[0].shape[0])
    compression = meta0.payload[0].compression
    encoding = meta0.payload[0].encoding
    print(
        f"  {count} messages, {grid_size}x{grid_size}, encoding={encoding}, compression={compression}"
    )
    print(f"  dt={dt}")

    if zoom:
        row_range, col_range = zoom
        print(
            f"  zoom: rows [{row_range[0]}:{row_range[1]}], cols [{col_range[0]}:{col_range[1]}]"
        )
        print(f"  using range decode (tensorjump)")

    indices = list(range(0, count, skip))
    print(f"Decoding {len(indices)} frames (every {skip}th)...")

    frames = []
    for i, idx in enumerate(indices):
        if zoom:
            frame = decode_zoom(f, idx, row_range, col_range, grid_size)
        else:
            _, arrays = f.decode_message(idx)
            frame = arrays[0][::4, ::4]
        frames.append(frame)
        if (i + 1) % 50 == 0:
            print(f"  decoded {i + 1}/{len(indices)}")

    print("Starting animation...")
    vmax = max(np.abs(frame).max() for frame in frames)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=50)
    im = ax.imshow(frames[0], cmap="seismic", origin="lower", vmin=-vmax, vmax=vmax)
    plt.colorbar(im, ax=ax, shrink=0.8)
    title = ax.set_title("t = 0.000")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    def update(frame_num):
        im.set_data(frames[frame_num])
        t_val = indices[frame_num] * dt
        title.set_text(f"t = {t_val:.3f}")
        return [im, title]

    ani = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=200, blit=True
    )

    if "--save" in sys.argv:
        out = tgm_path.replace(".tgm", ".gif")
        print(f"Rendering {len(frames)} frames to optimized GIF...")
        pil_frames = []
        for i in range(len(frames)):
            update(i)
            frame_buf = io.BytesIO()
            fig.savefig(frame_buf, format="png", bbox_inches="tight", pad_inches=0.3)
            frame_buf.seek(0)
            img = Image.open(frame_buf).convert("P", palette=Image.ADAPTIVE, colors=128)
            pil_frames.append(img)
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
