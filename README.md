# 2D Wave Equation Solver with MPI

Solves the 2D wave equation on a square domain with zero boundary conditions. Uses implicit time stepping with a Gauss-Seidel iterative solver, parallelized via MPI 2D Cartesian domain decomposition.

<img src="wave_simulation.gif?raw=true" width="400px">

## Quick Start

```bash
# Install MPI and NetCDF (Ubuntu/Debian)
sudo apt install libopenmpi-dev libnetcdf-dev

# Install MPI and NetCDF (Arch)
sudo pacman -S openmpi netcdf

# Build and run
mpicc -O3 -march=native waveEq.c -o wave.x -lnetcdf -lm
mpirun -np 4 ./wave.x

# View results
ncview grid.nc
```

## Output Formats

The solver supports three output formats, controlled by flags in `config.h`. Set to `1` to enable, `0` to disable.

| Flag | Default | Output file | Requires |
|---|---|---|---|
| `write_posix` | 0 | `t_NNNN.txt` | nothing |
| `write_netcdf` | 1 | `grid.nc` | `libnetcdf` |
| `write_tensogram` | 0 | `grid.tgm` | `libtensogram_ffi` |

### Building with Tensogram support

[Tensogram](https://github.com/ecmwf/tensogram) is an experimental N-Tensor binary format. To enable it:

1. Set `write_tensogram` to `1` in `config.h`
2. Build the Tensogram FFI library and copy its header:

```bash
cd /path/to/tensogram
cargo build --release -p tensogram-ffi
cp crates/tensogram-ffi/tensogram.h /path/to/2d_wave_mpi/
```

3. Compile with the Tensogram library:

```bash
TGM=/path/to/tensogram/target/release
mpicc -O3 -march=native -I. waveEq.c -o wave.x \
    -lnetcdf -L$TGM -ltensogram_ffi -lm -lpthread -ldl
```

4. Run with the library path:

```bash
LD_LIBRARY_PATH=$TGM mpirun -np 4 ./wave.x
```

### Building without any optional output

Set both `write_netcdf` and `write_tensogram` to `0` in `config.h`:

```bash
mpicc -O3 -march=native waveEq.c -o wave.x -lm
```

## Visualization

### NetCDF — use ncview or any NetCDF viewer

```bash
ncview grid.nc
```

### Tensogram — use the included Python visualizer

The visualizer uses the `tensogram` Python bindings. Install them first:

```bash
cd /path/to/tensogram
uv venv .venv && source .venv/bin/activate
uv pip install numpy matplotlib Pillow maturin
maturin develop --release -m crates/tensogram-python/Cargo.toml
```

Then visualize:

```bash
# Interactive window
python3 visualize_tgm.py grid.tgm

# Save as GIF (every 10th frame)
python3 visualize_tgm.py grid.tgm 10 --save

# Zoom into a sub-region using range decode (tensorjump)
python3 visualize_tgm.py grid.tgm 10 --zoom 100:400,100:400 --save
```

Usage: `visualize_tgm.py <file> [skip] [--save] [--zoom row_start:row_end,col_start:col_end]`

- `file` — path to the `.tgm` file (default: `grid.tgm`)
- `skip` — use every Nth frame (default: 10)
- `--save` — save as GIF instead of showing interactive window
- `--zoom` — decode only a rectangular sub-region via range decode (tensorjump). Only fetches the bytes needed for the specified rows and columns — no full-grid decode.

The `dt` value and grid size are read automatically from the Tensogram metadata.

## Running the Solver

The number of MPI ranks must be a perfect square (1, 4, 9, 16, ...) and the grid size (528 by default) must be divisible by `sqrt(ranks)`.

```bash
mpirun -np 1 ./wave.x     # single process
mpirun -np 4 ./wave.x     # 2x2 decomposition
mpirun -np 9 ./wave.x     # 3x3 decomposition
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
| `write_netcdf` | 1 | Enable NetCDF output |
| `write_tensogram` | 0 | Enable Tensogram output |

### Tensogram encoding

These are only used when `write_tensogram` is `1`:

| Parameter | Default | Description |
|---|---|---|
| `tgm_use_packing` | 1 | Use simple_packing lossy compression (0 = raw float64) |
| `tgm_bits_per_value` | 24 | Bits per value when packing (lower = smaller file, less precision) |
| `tgm_use_szip` | 1 | Apply szip (libaec) lossless compression on top of packing |
| `tgm_szip_rsi` | 128 | Records per segment for szip |
| `tgm_szip_block_size` | 16 | Samples per block for szip |
| `tgm_szip_flags` | 8 | libaec flags (8 = AEC_DATA_PREPROCESS) |

The encoding pipeline is: **raw float64 → simple_packing (lossy) → szip (lossless)**. Each stage is independently toggleable. With both enabled (default), a 528×528 float64 grid compresses from ~2.2 MB to ~0.6 MB per timestep.

Common `tgm_bits_per_value` choices:

| Bits | Compression vs raw | Precision loss | Use case |
|---|---|---|---|
| 24 | ~2.7x smaller | ~10⁻⁷ | Default, good balance |
| 16 | ~4x smaller | ~10⁻⁵ | Visualization, quick analysis |
| 12 | ~5.3x smaller | ~10⁻⁴ | Thumbnails, previews |
| 0 (packing off) | 1x (no compression) | None | Exact reproducibility |

When `tgm_use_szip` is enabled, the encoder automatically stores szip block offsets in the metadata. This enables **range decode (tensorjump)** — the visualizer's `--zoom` flag uses this to fetch only the rows/columns needed without decompressing the full grid.
