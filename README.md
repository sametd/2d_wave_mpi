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

## Project Structure

```
waveEq.c          Main solver (MPI topology, Gauss-Seidel, time loop)
config.h          Simulation parameters and output toggles
io_netcdf.h       NetCDF output (included when write_netcdf=1)
io_tensogram.h    Tensogram benchmark (included when write_tensogram=1)
tensogram.h       Tensogram C API header (from tensogram FFI build)
visualize_tgm.py  Animate a single .tgm file (PIL-optimized GIFs)
plot_bench.py     Benchmark report: comparison PNGs + per-codec GIFs
benchmark.md      Compression benchmark results and analysis
```

## Output Formats

Controlled by flags in `config.h`. Set to `1` to enable, `0` to disable.

| Flag | Default | Output file | Requires |
|---|---|---|---|
| `write_posix` | 0 | `t_NNNN.txt` | nothing |
| `write_netcdf` | 0 | `grid.nc` | `libnetcdf` |
| `write_tensogram` | 1 | `bench_*.tgm` | `libtensogram_ffi` |

### Building with Tensogram support

[Tensogram](https://github.com/ecmwf/tensogram) is a binary format for N-dimensional tensors with built-in compression. To enable it:

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
    -L$TGM -ltensogram_ffi -lm -lpthread -ldl
```

4. Run with the library path:

```bash
LD_LIBRARY_PATH=$TGM mpirun -np 16 ./wave.x
```

### Building with NetCDF

Set `write_netcdf` to `1` and `write_tensogram` to `0` in `config.h`:

```bash
mpicc -O3 -march=native waveEq.c -o wave.x -lnetcdf -lm
```

### Building without any optional output

Set both `write_netcdf` and `write_tensogram` to `0` in `config.h`:

```bash
mpicc -O3 -march=native waveEq.c -o wave.x -lm
```

## Compression Benchmark

When `write_tensogram` is enabled, the solver captures simulation frames and then benchmarks all 21 compression codecs automatically. Each codec writes a separate `bench_<name>.tgm` file and the results are printed as a table and saved to `results.md`.

The codecs tested include:

- **Lossless:** raw, lz4, zstd (levels 1/3/9), blosc2
- **Lossless + shuffle:** byte shuffle with zstd, lz4, blosc2
- **Lossy quantization:** simple_packing at 24/16/12 bits per value, combined with zstd, lz4, or szip
- **Floating-point lossy:** ZFP (fixed rate, fixed accuracy), SZ3 (absolute and relative error bounds)

See [benchmark.md](benchmark.md) for the full benchmark table with analysis and recommendations.

## Visualization

### NetCDF

```bash
ncview grid.nc
```

### Tensogram — single file animation

```bash
# Interactive window
python3 visualize_tgm.py grid.tgm

# Save as GIF (every 10th frame)
python3 visualize_tgm.py grid.tgm 10 --save

# Zoom into a sub-region using range decode
python3 visualize_tgm.py grid.tgm 10 --zoom 100:400,100:400 --save
```

### Tensogram — benchmark comparison

```bash
# Static comparison PNGs (all codecs side by side + error maps)
python3 plot_bench.py 150 --save

# Generate animated GIFs for every codec
python3 plot_bench.py 150 --save --gifs --skip 3
```

Both Python scripts require the `tensogram` Python bindings:

```bash
cd /path/to/tensogram
python3 -m venv .venv && source .venv/bin/activate
pip install numpy matplotlib Pillow maturin
cd crates/tensogram-python && maturin develop --release
```

## Running the Solver

The number of MPI ranks must be a perfect square (1, 4, 9, 16, ...) and the grid size (528 by default) must be divisible by `sqrt(ranks)`.

```bash
mpirun -np 1 ./wave.x     # single process
mpirun -np 4 ./wave.x     # 2x2 decomposition
mpirun -np 16 ./wave.x    # 4x4 decomposition (recommended)
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
