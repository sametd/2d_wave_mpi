# Tensogram Compression Benchmark

## Reproducing these results

```bash
# 1. Build tensogram FFI library (with all codecs)
cd /path/to/tensogram
cargo build --release -p tensogram-ffi

# 2. Copy the header and build the benchmark
cd /path/to/2d_wave_mpi
cp /path/to/tensogram/crates/tensogram-ffi/tensogram.h .
TGM=/path/to/tensogram/target/release
mpicc -O3 -march=native -I. waveEq.c -o wave_bench.x \
    -L$TGM -ltensogram_ffi -lm -lpthread -ldl

# 3. Run (adjust -np to a perfect square that fits your cores)
LD_LIBRARY_PATH=$TGM mpirun -np 16 ./wave_bench.x

# 4. Visualize
source /path/to/tensogram/.venv/bin/activate
python3 plot_bench.py 150 --save --gifs --skip 3
```

The benchmark runs the full simulation, captures frames, then encodes each frame with all 21 codec configurations. Results are printed to stdout and saved to this file. The codec table and parameters are defined in `io_tensogram.h`.

## Dataset

- **Simulation:** 2D wave equation, two Gaussian bumps on a square domain
- **Grid:** 528 x 528 points
- **Data type:** float64 (64-bit double precision, 8 bytes per value)
- **Simulation time:** t = 0 to 3.0 seconds (3001 solver steps, every 10th captured)
- **Frames captured:** 301
- **Raw size per frame:** 528 x 528 x 8 bytes = **2.23 MB**
- **Total raw size:** 301 frames x 2.23 MB = **671 MB**
- **MPI ranks:** 16 (4x4 decomposition)

## How to read this table

- **Size** and **Ratio** are compared to the raw uncompressed baseline (671 MB)
- **Encode MB/s** is raw data throughput (how fast the encoder processes raw input)
- **Max Error** is the worst-case absolute difference between any decoded value and the original
- **RMSE** is the root-mean-square error across all values
- **Lossless** means Max Error = 0 (bit-for-bit identical after decompression)

---

## 1. Lossless Compression

These preserve the original float64 values exactly. No precision loss.

Float64 data has high-entropy mantissa bits that look like random noise to byte-level compressors, so lossless ratios on scientific floats are inherently low.

| Config | What it does | Size | Ratio | Speed | Lossy? |
|--------|-------------|------|-------|-------|--------|
| **raw** | No compression, just tensogram framing | 671 MB | 1.0x | 451 MB/s | No |
| **lz4** | Fast lossless compressor (no effort spent searching for patterns) | 674 MB | 1.0x | 1078 MB/s | No |
| **zstd level 1** | Zstandard at minimum effort | 639 MB | 1.05x | 260 MB/s | No |
| **zstd level 3** | Zstandard at default effort | 639 MB | 1.05x | 253 MB/s | No |
| **zstd level 9** | Zstandard at high effort (diminishing returns on float64) | 640 MB | 1.05x | 350 MB/s | No |
| **blosc2 (zstd)** | Blosc2 meta-compressor using zstd internally, chunk-based | 526 MB | 1.28x | 223 MB/s | No |

**Key takeaway:** Lossless compressors barely help on float64 data. LZ4 is fastest but provides zero compression. Zstd levels 1-9 all give ~5% reduction because the bottleneck is float64 entropy, not search effort. Blosc2 does slightly better (1.28x) because its chunked format finds patterns within smaller blocks.

---

## 2. Lossless with Byte Shuffle

Byte shuffling rearranges the 8 bytes of each float64 value so that all first bytes are together, all second bytes together, etc. This groups the exponent bytes (which are similar across values of similar magnitude) and makes them compressible.

| Config | What it does | Size | Ratio | Speed | Lossy? |
|--------|-------------|------|-------|-------|--------|
| **shuffle + zstd 3** | Byte shuffle then zstd | 528 MB | 1.27x | 497 MB/s | No |
| **shuffle + lz4** | Byte shuffle then lz4 | 545 MB | 1.23x | 627 MB/s | No |
| **shuffle + blosc2 (lz4)** | Byte shuffle then blosc2 with lz4 | 544 MB | 1.23x | 491 MB/s | No |
| **shuffle + blosc2 (zstd)** | Byte shuffle then blosc2 with zstd | 526 MB | 1.28x | 188 MB/s | No |

**Key takeaway:** Shuffle helps, but only gets to ~1.3x. Best lossless option is **shuffle + zstd** (1.27x at 497 MB/s). If you need exact float64 preservation, this is the ceiling.

---

## 3. Lossy: GRIB-Style Quantization (simple_packing)

Reduces each float64 to a fixed number of bits by computing `(value - min) / scale` and storing the integer result. The number after "pack" is bits per value (out of the original 64).

| Config | Bits | What it does | Size | Ratio | Speed | Max Error | RMSE |
|--------|------|-------------|------|-------|-------|-----------|------|
| **pack24** | 24 | Quantize to 24 bits, no further compression | 252 MB | 2.7x | 595 MB/s | 3.0e-08 | 1.7e-08 |
| **pack24 + lz4** | 24 | Quantize then lz4 | 243 MB | 2.8x | 572 MB/s | 3.0e-08 | 1.7e-08 |
| **pack24 + zstd** | 24 | Quantize then zstd | 236 MB | 2.8x | 426 MB/s | 3.0e-08 | 1.7e-08 |
| **pack16 + szip** | 16 | Quantize then szip (CCSDS standard) | 157 MB | 4.3x | 522 MB/s | 7.6e-06 | 3.6e-06 |
| **pack16 + zstd** | 16 | Quantize then zstd | 147 MB | 4.6x | 495 MB/s | 7.6e-06 | 3.6e-06 |
| **pack12 + zstd** | 12 | Quantize then zstd | 88 MB | 7.6x | 423 MB/s | 2.4e-04 | 1.1e-04 |

**What the bit depths mean for a wave amplitude range of ~[-0.3, 1.0]:**

| Bits | Precision | Significant digits | Good for |
|------|-----------|-------------------|----------|
| 24 | ~30 nanoseconds | ~7 digits | Publication-quality, indistinguishable from raw |
| 16 | ~8 microseconds | ~5 digits | Routine analysis, visualization |
| 12 | ~0.2 milliseconds | ~3-4 digits | Quick previews, screening |

**Key takeaway:** **pack16 + zstd** is the sweet spot for most scientific workflows: 4.6x compression with ~5 digits of precision and 495 MB/s throughput.

---

## 4. Lossy: Floating-Point Compressors (ZFP, SZ3)

These are designed specifically for scientific floating-point arrays. They exploit spatial correlations between neighboring values, not just bit patterns.

### ZFP

| Config | Mode | What it does | Size | Ratio | Speed | Max Error | RMSE |
|--------|------|-------------|------|-------|-------|-----------|------|
| **zfp rate 16** | Fixed rate | Each value gets exactly 16 bits (guaranteed 4x ratio) | 168 MB | 4.0x | 292 MB/s | 5.4e-05 | 3.1e-06 |
| **zfp rate 8** | Fixed rate | Each value gets exactly 8 bits (guaranteed 8x ratio) | 84 MB | 8.0x | 515 MB/s | 1.3e-02 | 9.9e-04 |
| **zfp tolerance 1e-6** | Fixed accuracy | Error guaranteed below 1e-6 per value | 188 MB | 3.6x | 263 MB/s | 5.7e-07 | 1.3e-07 |

- **Fixed rate** = you choose the size, error varies
- **Fixed accuracy** = you choose the max error, size varies

### SZ3

| Config | Mode | What it does | Size | Ratio | Speed | Max Error | RMSE |
|--------|------|-------------|------|-------|-------|-----------|------|
| **sz3 absolute 1e-6** | Absolute error bound | Every value within 1e-6 of original | 74 MB | 9.1x | 217 MB/s | 1.0e-06 | 5.1e-07 |
| **sz3 relative 1e-4** | Relative error bound | Every value within 0.01% of original | 21 MB | 31.8x | 339 MB/s | 1.0e-04 | 4.9e-05 |

**Key takeaway:** SZ3 with relative error bound achieves **31.8x compression** while keeping every value within 0.01% of the original. For error-bounded science, **sz3 absolute 1e-6** gives 9.1x with a hard guarantee that no value deviates by more than 0.000001.

---

## Summary: Which codec should I use?

| Use case | Recommended | Ratio | Error | Speed |
|----------|------------|-------|-------|-------|
| Exact reproduction required | shuffle + zstd | 1.3x | 0 | 497 MB/s |
| Near-lossless archival | pack24 + zstd | 2.8x | < 30 ns | 426 MB/s |
| Day-to-day analysis | pack16 + zstd | 4.6x | < 8 us | 495 MB/s |
| Error-bounded science | sz3 absolute 1e-6 | 9.1x | = 1e-6 | 217 MB/s |
| Visualization / screening | sz3 relative 1e-4 | 31.8x | 0.01% | 339 MB/s |
| Maximum throughput | lz4 | 1.0x | 0 | 1078 MB/s |

---

## Notes

- **Why is lz4 slightly larger than raw?** LZ4 cannot compress the random-looking float64 mantissa bits, so the output is the original data plus LZ4's framing overhead (~8.7 KB per frame).
- **Why do zstd levels 1, 3, and 9 all give the same ratio?** The bottleneck is float64 entropy, not search depth. Higher levels spend more CPU but find no additional patterns in random-looking mantissa bytes.
- **szip** (CCSDS 121.0-B-3) is the space agency standard used in satellite data. It works well with byte-aligned quantized data (16-bit, 32-bit) but does not support raw 64-bit floats.
