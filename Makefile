#  2D Wave Equation MPI Solver — Makefile
#
#  The tensogram C FFI library is downloaded as a prebuilt release
#  asset from GitHub (https://github.com/ecmwf/tensogram/releases)
#  into the local `.tensogram/` prefix. Run `make tensogram` once
#  (or let `make` run it as a dependency).
#
#  Common targets:
#    make                    Build wave.x (auto-downloads tensogram)
#    make run                Build and run with default MPI ranks
#    make bench              Build, run, then plot benchmarks (sequential)
#    make viz FILE=bench_raw.tgm   Visualize a .tgm file
#    make gifs               Generate animated GIFs for all codecs
#    make tensogram          Download + unpack the prebuilt FFI tarball
#    make clean              Remove build artifacts (keep .tensogram)
#    make distclean          Remove everything, including .tensogram
#    make print-tgm          Debug: show resolved tensogram paths

CC       := mpicc
CFLAGS   := -O3 -march=native -Wall -Wextra
LDLIBS   := -lm -lpthread -ldl

#  ── tensogram prebuilt release asset ──────────────────────────────
TGM_VERSION ?= 0.21.0

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_S),Darwin)
  TGM_PLATFORM := macos-aarch64
  TGM_SOEXT    := dylib
else ifeq ($(UNAME_M),aarch64)
  TGM_PLATFORM := linux-aarch64
  TGM_SOEXT    := so
else
  TGM_PLATFORM := linux-x86_64
  TGM_SOEXT    := so
endif

TGM_PREFIX  := .tensogram
TGM_LIBDIR  := $(TGM_PREFIX)/lib
TGM_LIB     := $(TGM_LIBDIR)/libtensogram.$(TGM_SOEXT)
TGM_INCLUDE := $(TGM_PREFIX)/include
TGM_ASSET   := tensogram-ffi-$(TGM_VERSION)-$(TGM_PLATFORM).tar.gz
TGM_URL     := https://github.com/ecmwf/tensogram/releases/download/$(TGM_VERSION)/$(TGM_ASSET)

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

RUNENV   := LD_LIBRARY_PATH=$(abspath $(TGM_LIBDIR)):$$LD_LIBRARY_PATH

TARGET   := wave.x
SRC      := waveEq.c
HEADERS  := config.h io_tensogram.h io_netcdf.h

.PHONY: all run bench viz viz-save gifs plots clean distclean tensogram py-setup print-tgm

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS) $(TGM_LIB)
	$(CC) $(CFLAGS) -I. -I$(TGM_INCLUDE) $(SRC) -o $@ \
	    -L$(abspath $(TGM_LIBDIR)) -ltensogram \
	    -Wl,-rpath,$(abspath $(TGM_LIBDIR)) $(LDLIBS)

#  ── tensogram prebuilt C-FFI tarball from GitHub releases ─────────
#  Each release ships `tensogram-ffi-<version>-<platform>.tar.gz`
#  containing lib/libtensogram.{so,dylib}, lib/libtensogram.a, a
#  pkg-config descriptor, and include/tensogram/tensogram.h. We
#  unpack it into the local `.tensogram/` prefix — no Rust
#  toolchain, no sudo, no /usr/local pollution.

tensogram: $(TGM_LIB)

$(TGM_LIB):
	@echo "Fetching $(TGM_ASSET) ..."
	@mkdir -p $(TGM_PREFIX)
	curl -fsSL "$(TGM_URL)" | tar -xz -C $(TGM_PREFIX)
	@touch $(TGM_LIB)
	@echo "  -> $(TGM_LIB)"
	@echo "  -> $(TGM_INCLUDE)/tensogram/tensogram.h"

print-tgm:
	@echo "TGM_VERSION  : $(TGM_VERSION)"
	@echo "TGM_PLATFORM : $(TGM_PLATFORM)"
	@echo "TGM_URL      : $(TGM_URL)"
	@echo "TGM_PREFIX   : $(abspath $(TGM_PREFIX))"
	@echo "TGM_LIB      : $(TGM_LIB)"
	@echo "TGM_INCLUDE  : $(TGM_INCLUDE)"
	@echo "RUNENV       : $(RUNENV)"

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
	rm -f bench_*.tgm bench_*.gif bench_*.png grid.tgm benchmark.md
	rm -rf $(TGM_PREFIX) $(VENV)

py-setup:
	python3 -m venv $(VENV)
	$(VENV)/bin/pip install --upgrade pip
	$(VENV)/bin/pip install 'tensogram==$(TGM_VERSION)' numpy matplotlib Pillow
