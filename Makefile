# 2D Wave Equation MPI Solver — Makefile
#
# Usage:
#   make                    Build the solver
#   make run                Build and run with default MPI ranks
#   make bench              Build, run, then plot benchmarks
#   make viz                Visualize the raw output
#   make gifs               Generate animated GIFs for all codecs
#   make clean              Remove build artifacts
#   make distclean          Remove build artifacts + simulation output

# ── Toolchain ────────────────────────────────────────────────────────
CC       := mpicc
CFLAGS   := -O3 -march=native -Wall -Wextra
LDLIBS   := -ltensogram_ffi -lm -lpthread -ldl

# ── Tensogram FFI path (override: make TGM_DIR=/your/path) ──────────
TGM_DIR  ?= $(HOME)/Desktop/RUST/tensogram-remote/target/release

CFLAGS   += -I.
LDFLAGS  += -L$(TGM_DIR)

# ── NetCDF (only linked when write_netcdf=1 in config.h) ────────────
NETCDF   ?= 0
ifeq ($(NETCDF),1)
  LDLIBS += -lnetcdf
endif

# ── MPI configuration ───────────────────────────────────────────────
NP       ?= 16
MPIFLAGS ?= --oversubscribe

# ── Python (for visualization scripts) ──────────────────────────────
PYTHON   ?= python3
SKIP     ?= 10
FRAME    ?= 150

# ── Targets ─────────────────────────────────────────────────────────
TARGET   := wave.x
SRC      := waveEq.c
HEADERS  := config.h io_tensogram.h io_netcdf.h tensogram.h

.PHONY: all run bench viz gifs plots clean distclean header

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ $(LDLIBS)

run: $(TARGET)
	LD_LIBRARY_PATH=$(TGM_DIR) mpirun $(MPIFLAGS) -np $(NP) ./$(TARGET)

bench: run plots gifs

# ── Visualization targets ───────────────────────────────────────────
plots: run
	LD_LIBRARY_PATH=$(TGM_DIR) MPLBACKEND=Agg $(PYTHON) plot_bench.py $(FRAME) --save

gifs:
	LD_LIBRARY_PATH=$(TGM_DIR) MPLBACKEND=Agg $(PYTHON) plot_bench.py $(FRAME) --gifs --skip $(SKIP)

viz:
	LD_LIBRARY_PATH=$(TGM_DIR) $(PYTHON) visualize_tgm.py grid.tgm $(SKIP)

viz-save:
	LD_LIBRARY_PATH=$(TGM_DIR) MPLBACKEND=Agg $(PYTHON) visualize_tgm.py grid.tgm $(SKIP) --save

# ── Copy tensogram.h from FFI build ────────────────────────────────
header:
	cp $(TGM_DIR)/../../crates/tensogram-ffi/tensogram.h .

# ── Cleanup ─────────────────────────────────────────────────────────
clean:
	rm -f $(TARGET)

distclean: clean
	rm -f bench_*.tgm bench_*.gif bench_*.png grid.tgm benchmark.md
