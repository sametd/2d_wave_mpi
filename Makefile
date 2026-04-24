#  2D Wave Equation MPI Solver — Makefile
#
#  The tensogram C FFI library is fetched from crates.io via a small
#  helper crate in `vendor/`. Run `make tensogram` once (or let `make`
#  run it as a dependency) to build `libtensogram_ffi.so` and drop
#  `tensogram.h` into the project root.
#
#  Common targets:
#    make                    Build wave.x (auto-builds tensogram)
#    make run                Build and run with default MPI ranks
#    make bench              Build, run, then plot benchmarks (sequential)
#    make viz FILE=bench_raw.tgm   Visualize a .tgm file
#    make gifs               Generate animated GIFs for all codecs
#    make tensogram          Build libtensogram_ffi.so + copy header
#    make clean              Remove build artifacts (keep vendor cache)
#    make distclean          Remove everything, including vendor/target
#    make print-tgm          Debug: show resolved tensogram paths

CC       := mpicc
CFLAGS   := -O3 -march=native -Wall -Wextra
LDLIBS   := -lm -lpthread -ldl

VENDOR_DIR := vendor
TGM_BUILD  := $(VENDOR_DIR)/target/release
TGM_LIB    := $(TGM_BUILD)/libtensogram_ffi.so
TGM_HEADER := tensogram.h

NETCDF   ?= 0
ifeq ($(NETCDF),1)
  LDLIBS += -lnetcdf
endif

NP       ?= 16
MPIFLAGS ?= --oversubscribe

VENV     ?= .venv
ifneq ($(wildcard $(VENV)/bin/python3),)
  PYTHON ?= $(VENV)/bin/python3
else
  PYTHON ?= python3
endif

SKIP     ?= 10
FRAME    ?= 150
FILE     ?= bench_raw.tgm

RUNENV   := LD_LIBRARY_PATH=$(abspath $(TGM_BUILD)):$$LD_LIBRARY_PATH

TARGET   := wave.x
SRC      := waveEq.c
HEADERS  := config.h io_tensogram.h io_netcdf.h $(TGM_HEADER)

.PHONY: all run bench viz viz-save gifs plots clean distclean tensogram py-setup print-tgm

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS) $(TGM_LIB)
	$(CC) $(CFLAGS) -I. $(SRC) -o $@ \
	    -L$(abspath $(TGM_BUILD)) -ltensogram_ffi \
	    -Wl,-rpath,$(abspath $(TGM_BUILD)) $(LDLIBS)

#  ── tensogram from crates.io via the vendor helper crate ──────────
#  `cargo build` downloads `tensogram-ffi` into the shared Cargo cache
#  the first time around, then builds its cdylib. The hashed artefact
#  under `target/release/deps/` is then stabilised via a symlink so
#  `mpicc -ltensogram_ffi` can find it.

tensogram: $(TGM_LIB) $(TGM_HEADER)

$(TGM_LIB): $(VENDOR_DIR)/Cargo.toml $(VENDOR_DIR)/src/lib.rs
	cd $(VENDOR_DIR) && cargo build --release
	@ln -sf "$$(ls $(abspath $(VENDOR_DIR))/target/release/deps/libtensogram_ffi-*.so | head -1)" $(TGM_LIB)
	@echo "  -> $(TGM_LIB)"

#  `cargo metadata` locates the tensogram-ffi crate root in the
#  registry cache; the shipped `tensogram.h` sits next to Cargo.toml.

$(TGM_HEADER): $(TGM_LIB)
	@hdr="$$(cargo metadata --manifest-path $(VENDOR_DIR)/Cargo.toml --format-version=1 \
	    | $(PYTHON) -c "import sys,json,os;m=json.load(sys.stdin);print(next(os.path.join(os.path.dirname(p['manifest_path']),'tensogram.h') for p in m['packages'] if p['name']=='tensogram-ffi'))")"; \
	cp "$$hdr" $(TGM_HEADER); \
	echo "  -> $(TGM_HEADER) (from $$hdr)"

print-tgm:
	@echo "VENDOR_DIR : $(VENDOR_DIR)"
	@echo "TGM_BUILD  : $(abspath $(TGM_BUILD))"
	@echo "TGM_LIB    : $(TGM_LIB)"
	@echo "TGM_HEADER : $(TGM_HEADER)"
	@echo "RUNENV     : $(RUNENV)"

run: $(TARGET)
	$(RUNENV) mpirun $(MPIFLAGS) -np $(NP) ./$(TARGET)

.NOTPARALLEL: bench
bench: run plots gifs

plots:
	$(RUNENV) MPLBACKEND=Agg $(PYTHON) plot_bench.py $(FRAME) --save

gifs:
	$(RUNENV) MPLBACKEND=Agg $(PYTHON) plot_bench.py $(FRAME) --gifs --skip $(SKIP)

viz:
	$(RUNENV) $(PYTHON) visualize_tgm.py $(FILE) $(SKIP)

viz-save:
	$(RUNENV) MPLBACKEND=Agg $(PYTHON) visualize_tgm.py $(FILE) $(SKIP) --save

clean:
	rm -f $(TARGET)

distclean: clean
	rm -f bench_*.tgm bench_*.gif bench_*.png grid.tgm benchmark.md $(TGM_HEADER)
	rm -rf $(VENDOR_DIR)/target $(VENDOR_DIR)/Cargo.lock $(VENV)

py-setup:
	python3 -m venv $(VENV)
	$(VENV)/bin/pip install --upgrade pip
	$(VENV)/bin/pip install tensogram numpy matplotlib Pillow
