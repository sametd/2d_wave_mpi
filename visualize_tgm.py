#!/usr/bin/env python3
"""Visualize a 2D wave equation .tgm file produced by the MPI solver.

Requires: tensogram (maturin develop), numpy, matplotlib, Pillow
"""

import io
import os
import sys

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import tensogram
from PIL import Image


def main():
    tgm_path = sys.argv[1] if len(sys.argv) > 1 else "grid.tgm"
    skip = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    print(f"Loading {tgm_path}...")
    f = tensogram.TensogramFile.open(tgm_path)
    count = f.message_count()
    print(f"Found {count} messages")

    meta0, _ = f.decode_message(0)
    dt = meta0["wave"]["dt"]
    print(f"dt = {dt} (from metadata)")

    indices = list(range(0, count, skip))
    print(f"Decoding {len(indices)} frames (every {skip}th)...")

    frames = []
    for i, idx in enumerate(indices):
        meta, arrays = f.decode_message(idx)
        frames.append(arrays[0][::4, ::4])
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
