#!/usr/bin/env python3
"""Benchmark report: static comparison + per-codec animated GIFs.

Usage:
    python3 plot_bench.py [frame_index] [--save] [--gifs]

--save   Write bench_comparison.png and bench_errors.png
--gifs   Write one animated GIF per codec (using matplotlib, not PIL)

Requires: tensogram, numpy, matplotlib
"""

import argparse
import glob
import os

import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import tensogram


def load_frame(path, frame_idx):
    """Decode a single frame from a .tgm file."""
    f = tensogram.TensogramFile.open(path)
    idx = min(frame_idx, f.message_count() - 1)
    _, objects = f.decode_message(idx)
    _desc, arr = objects[0]
    return arr


def load_all_frames(path, skip=5):
    """Decode all frames (with skip) for animation."""
    f = tensogram.TensogramFile.open(path)
    count = f.message_count()
    frames = []
    for idx in range(0, count, skip):
        _, objects = f.decode_message(idx)
        _desc, arr = objects[0]
        frames.append(arr)
    return frames


def plot_comparison(frames, frame_idx, raw_arr, save):
    """Plot all codecs side by side in a grid."""
    n = len(frames)
    ncols = 5
    nrows = (n + ncols - 1) // ncols
    vmax = np.abs(raw_arr).max() if raw_arr is not None else max(
        np.abs(a).max() for _, a in frames if a is not None
    )

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5), squeeze=False)
    fig.suptitle(f"Tensogram Codec Comparison \u2014 Frame {frame_idx}", fontsize=16, y=1.02)

    for i, (name, arr) in enumerate(frames):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        if arr is not None:
            ax.imshow(arr, cmap="seismic", origin="lower", vmin=-vmax, vmax=vmax)
            ax.set_title(name, fontsize=10, fontweight="bold")
        else:
            ax.text(0.5, 0.5, "FAILED", ha="center", va="center",
                    transform=ax.transAxes, fontsize=14, color="red")
            ax.set_title(name, fontsize=10, color="red")
        ax.set_xticks([])
        ax.set_yticks([])

    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r][c].set_visible(False)

    fig.tight_layout()
    if save:
        out = "bench_comparison.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved {out} ({os.path.getsize(out) / 1024:.0f} KB)")
    return fig


def plot_errors(frames, raw_arr, frame_idx, save):
    """Plot error heatmaps relative to the raw baseline."""
    if raw_arr is None:
        return None

    lossy = [(n, a) for n, a in frames if a is not None and n != "raw"]
    n = len(lossy)
    ncols = 5
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5), squeeze=False)
    fig.suptitle(f"Absolute Error vs Raw \u2014 Frame {frame_idx}", fontsize=16, y=1.02)

    for i, (name, arr) in enumerate(lossy):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        diff = np.abs(arr - raw_arr)
        max_err = diff.max()
        if max_err > 0:
            im = ax.imshow(diff, cmap="hot", origin="lower")
            plt.colorbar(im, ax=ax, shrink=0.7, format="%.1e")
        else:
            ax.imshow(np.zeros_like(diff), cmap="hot", origin="lower")
            ax.text(0.5, 0.5, "EXACT", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="lime", fontweight="bold")
        ax.set_title(f"{name}\nmax={max_err:.2e}", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r][c].set_visible(False)

    fig.tight_layout()
    if save:
        out = "bench_errors.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved {out} ({os.path.getsize(out) / 1024:.0f} KB)")
    return fig


def make_gif(path, name, skip=5):
    """Create an animated GIF for one codec using matplotlib's writer."""
    print(f"  {name}: loading frames...", end="", flush=True)
    all_frames = load_all_frames(path, skip)
    if not all_frames:
        print(" empty, skipping")
        return

    vmax = max(np.abs(fr).max() for fr in all_frames)
    fig, ax = plt.subplots(figsize=(6, 6), dpi=80)
    im = ax.imshow(all_frames[0], cmap="seismic", origin="lower", vmin=-vmax, vmax=vmax)
    plt.colorbar(im, ax=ax, shrink=0.8)
    title = ax.set_title(name, fontsize=14, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()

    def update(i):
        im.set_data(all_frames[i])
        title.set_text(f"{name}  (frame {i * skip})")
        return [im, title]

    ani = animation.FuncAnimation(fig, update, frames=len(all_frames), interval=100, blit=True)

    out = f"bench_{name}.gif"
    ani.save(out, writer="pillow", fps=15)
    plt.close(fig)
    size_kb = os.path.getsize(out) / 1024
    print(f" {len(all_frames)} frames -> {out} ({size_kb:.0f} KB)")


def main():
    parser = argparse.ArgumentParser(description="Benchmark visualization")
    parser.add_argument("frame", nargs="?", type=int, default=50, help="Frame index for static plots")
    parser.add_argument("--save", action="store_true", help="Save static PNG plots")
    parser.add_argument("--gifs", action="store_true", help="Generate per-codec animated GIFs")
    parser.add_argument("--skip", type=int, default=5, help="Frame skip for GIFs (default: 5)")
    args = parser.parse_args()

    if args.gifs:
        matplotlib.use("Agg")

    files = sorted(glob.glob("bench_*.tgm"))
    if not files:
        print("No bench_*.tgm files found. Run the benchmark first.")
        return

    names = [os.path.basename(f).replace("bench_", "").replace(".tgm", "") for f in files]
    print(f"Found {len(files)} benchmark files")

    # ── Static plots ──
    print(f"\nLoading frame {args.frame} from each codec...")
    frames = []
    raw_arr = None
    for path, name in zip(files, names):
        try:
            arr = load_frame(path, args.frame)
            frames.append((name, arr))
            if name == "raw":
                raw_arr = arr
            print(f"  {name}: {arr.shape}, range [{arr.min():.6f}, {arr.max():.6f}]")
        except Exception as e:
            print(f"  {name}: FAILED ({e})")
            frames.append((name, None))

    plot_comparison(frames, args.frame, raw_arr, args.save)
    plot_errors(frames, raw_arr, args.frame, args.save)

    # ── Animated GIFs ──
    if args.gifs:
        print(f"\nGenerating animated GIFs (skip={args.skip})...")
        for path, name in zip(files, names):
            try:
                make_gif(path, name, args.skip)
            except Exception as e:
                print(f"  {name}: FAILED ({e})")

    if not args.save and not args.gifs:
        plt.show()

    if args.save or args.gifs:
        print("\nDone.")


if __name__ == "__main__":
    main()
