# 2D Wave Equation Solver with MPI

Solves the 2D wave equation on a square domain with zero boundary conditions. Uses implicit time stepping with a Gauss-Seidel iterative solver, parallelized via MPI 2D Cartesian domain decomposition.

<img src="wave_simulation.gif?raw=true" width="400px">

## Quick Start

```bash
# 1. Install MPI (Ubuntu/Debian)
sudo apt install libopenmpi-dev
#    or (Arch)
sudo pacman -S openmpi

# 2. Install the Rust toolchain (once)
#    https://rustup.rs
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# 3. Build libtensogram_ffi.so + fetch tensogram.h (from crates.io)
make tensogram

# 4. Create the Python venv and install tensogram from PyPI
make py-setup

# 5. Build and run (16 MPI ranks by default)
make run
```

`make` alone will also build everything; the separate targets above just
make the two one-time installation steps visible.

## Build System

The project uses a Makefile with configurable variables:

| Variable | Default | Description |
|---|---|---|
| `NP` | `16` | Number of MPI ranks (must be a perfect square) |
| `VENV` | `.venv` | Python virtualenv directory (auto-used if present) |
| `PYTHON` | `$(VENV)/bin/python3` or `python3` | Python interpreter for visualization scripts |
| `SKIP` | `10` | Frame skip for GIF generation |
| `FRAME` | `150` | Frame index for static benchmark plots |
| `FILE` | `bench_raw.tgm` | Input file for `viz`/`viz-save` targets |
| `NETCDF` | `0` | Set to `1` to link NetCDF |

### Make Targets

```bash
make                 # Build wave.x (auto-builds tensogram if needed)
make tensogram       # Build libtensogram_ffi.so + copy tensogram.h
make py-setup        # Create .venv and pip install tensogram, numpy, matplotlib, Pillow
make run             # Build + run simulation
make bench           # Full pipeline: run + plots + gifs (sequential)
make plots           # Static benchmark PNGs (comparison + errors)
make gifs            # Animated GIF per codec
make viz             # Interactive visualizer (default: bench_raw.tgm)
make viz-save        # Save visualizer output as GIF
make print-tgm       # Debug: show resolved tensogram paths
make clean           # Remove wave.x
make distclean       # Remove wave.x + outputs + vendor/target + venv
```

Override any variable: `make run NP=4`, `make gifs SKIP=5`, `make viz FILE=bench_zstd-3.tgm`.

## Project Structure

```
Makefile            Build, run, and visualization targets
waveEq.c            Main solver (MPI topology, Gauss-Seidel, time loop)
config.h            Simulation parameters and output toggles
io_tensogram.h      Tensogram benchmark (included when write_tensogram=1)
io_netcdf.h         NetCDF output (included when write_netcdf=1)
tensogram.h         Tensogram C API header (gitignored, copied by `make tensogram`)
vendor/             Cargo helper that fetches tensogram-ffi from crates.io
  Cargo.toml          (depends on `tensogram-ffi = "*"`)
  src/lib.rs          (intentionally empty)
visualize_tgm.py    Animate a single .tgm file (PIL-optimized GIFs)
plot_bench.py       Benchmark report: comparison PNGs + per-codec GIFs
benchmark.md        Compression benchmark results (auto-generated)
```

## Output Formats

Controlled by flags in `config.h`. Set to `1` to enable, `0` to disable.

| Flag | Default | Output file | Requires |
|---|---|---|---|
| `write_posix` | 0 | `t_NNNN.txt` | nothing |
| `write_netcdf` | 0 | `grid.nc` | `libnetcdf` |
| `write_tensogram` | 1 | `bench_*.tgm` | `libtensogram_ffi` (from crates.io) |

### Tensogram (default)

[Tensogram](https://github.com/ecmwf/tensogram) is a binary format for
N-dimensional tensors with built-in compression, published as:

- `tensogram-ffi` on [crates.io](https://crates.io/crates/tensogram-ffi) — the C FFI library
- `tensogram` on [PyPI](https://pypi.org/project/tensogram/) — the Python bindings

The C side is fetched and built via a tiny helper crate in `vendor/`:

```bash
make tensogram
```

Under the hood this runs `cargo build --release` on `vendor/Cargo.toml`
(which declares `tensogram-ffi = "*"`), symlinks the produced cdylib into
`vendor/target/release/libtensogram_ffi.so`, and copies `tensogram.h`
from the Cargo registry cache into the project root. No manual version
management, no path environment variables, no local tensogram checkout.

To pin a specific version (e.g. for reproducible benchmarks), edit
`vendor/Cargo.toml`:

```toml
[dependencies]
tensogram-ffi = "=0.18.1"
```

### NetCDF

Set `write_netcdf` to `1` and `write_tensogram` to `0` in `config.h`, then:

```bash
make NETCDF=1
```

### No optional output

Set both `write_netcdf` and `write_tensogram` to `0` in `config.h`:

```bash
mpicc -O3 -march=native waveEq.c -o wave.x -lm
```

## Compression Benchmark

When `write_tensogram` is enabled, the solver captures simulation frames
and benchmarks all 21 compression codecs. Each codec writes a separate
`bench_<name>.tgm` file and results are printed as a Markdown table and
saved to `benchmark.md`.

The codecs tested include:

- **Lossless:** raw, lz4, zstd (levels 1/3/9), blosc2
- **Lossless + shuffle:** byte shuffle with zstd, lz4, blosc2
- **Lossy quantization:** simple_packing at 24/16/12 bits per value, combined with zstd, lz4, or szip
- **Floating-point lossy:** ZFP (fixed rate, fixed accuracy), SZ3 (absolute and relative error bounds)

See [benchmark.md](benchmark.md) for the full results table.

## Visualization

The Python scripts need the `tensogram` PyPI package plus numpy /
matplotlib / Pillow. `make py-setup` installs all of these into a local
`.venv/`:

```bash
make py-setup
```

Manual equivalent:

```bash
python3 -m venv .venv
.venv/bin/pip install tensogram numpy matplotlib Pillow
```

### Single file animation

```bash
make viz                              # Interactive window (bench_raw.tgm)
make viz FILE=bench_zstd-3.tgm        # Visualize a specific codec
make viz-save                         # Save as GIF
make viz-save SKIP=5                  # Save with custom frame skip

# Zoom into a sub-region using range decode
.venv/bin/python3 visualize_tgm.py bench_raw.tgm 10 --zoom 100:400,100:400 --save
```

### Benchmark comparison

```bash
make plots                            # Static PNGs (codec grid + error heatmaps)
make gifs                             # Animated GIF per codec
make bench                            # Full pipeline: run + plots + gifs
```

### NetCDF

```bash
ncview grid.nc
```

## Running the Solver

The number of MPI ranks must be a perfect square (1, 4, 9, 16, ...) and
the grid size (528 by default) must be divisible by `sqrt(ranks)`.

```bash
make run NP=1                         # single process
make run NP=4                         # 2x2 decomposition
make run NP=16                        # 4x4 decomposition (recommended)
```

## Configuration

All parameters are in `config.h`.

### Simulation

| Parameter | Default | Description |
|---|---|---|
| `dx` | 0.00594998 | Grid spacing (gives 528x528 grid) |
| `dt` | 0.001 | Time step |
| `f_to_go` | 3.14159 (pi) | Domain size |
| `t_to_go` | 3.0 | Total simulation time |
| `gerr_max` | 1e-6 | Gauss-Seidel convergence tolerance |

### Output

| Parameter | Default | Description |
|---|---|---|
| `write_posix` | 0 | Enable POSIX text file output |
| `write_netcdf` | 0 | Enable NetCDF output |
| `write_tensogram` | 1 | Enable Tensogram compression benchmark |
| `write_interval` | 10 | Write every Nth solver step (reduces output size) |
| `tgm_bench_steps` | 0 | Max frames to capture (0 = all) |

### Output size control

The solver runs 3001 steps for t=3.0. `write_interval` controls how many are captured:

| write_interval | Frames | Raw size | Covers |
|---|---|---|---|
| 1 | 3001 | 6.69 GB | every step |
| 10 | 301 | 671 MB | every 10th step |
| 20 | 151 | 337 MB | every 20th step |
