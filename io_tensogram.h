#ifndef IO_TENSOGRAM_H
#define IO_TENSOGRAM_H

#include "config.h"

#if write_tensogram

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <mpi.h>
#include "tensogram.h"

#define TGM_CHECK(e)                                                            \
    {                                                                           \
        tgm_error tgm_err = (e);                                               \
        if (tgm_err != TGM_ERROR_OK) {                                         \
            const char *tgm_detail = tgm_last_error();                          \
            const char *tgm_code   = tgm_error_string(tgm_err);                \
            fprintf(stderr, "Tensogram error [%s]: %s\n",                       \
                    tgm_code, tgm_detail ? tgm_detail : "(no details)");        \
            MPI_Abort(MPI_COMM_WORLD, (int)tgm_err);                            \
        }                                                                       \
    }

/* ── Codec configuration ── */

typedef struct {
    const char *name;
    const char *encoding;       /* "none" or "simple_packing" */
    int bits_per_value;
    const char *filter;         /* "none" or "shuffle" */
    const char *compression;    /* "none","szip","zstd","lz4","blosc2","zfp","sz3" */
    const char *comp_params;    /* JSON fragment, may contain %zu for num_values */
} tgm_codec_config_t;

static const tgm_codec_config_t TGM_CODECS[] = {
    /* lossless */
    { "raw",              "none", 0, "none", "none", "" },
    { "zstd-1",           "none", 0, "none", "zstd", "\"zstd_level\":1" },
    { "zstd-3",           "none", 0, "none", "zstd", "\"zstd_level\":3" },
    { "zstd-9",           "none", 0, "none", "zstd", "\"zstd_level\":9" },
    { "lz4",              "none", 0, "none", "lz4",  "" },
    { "blosc2-zstd",      "none", 0, "none", "blosc2",
      "\"blosc2_codec\":\"zstd\",\"blosc2_clevel\":5" },
    /* lossless + shuffle */
    { "shuffle+zstd-3",   "none", 0, "shuffle", "zstd", "\"zstd_level\":3" },
    { "shuffle+lz4",      "none", 0, "shuffle", "lz4",  "" },
    { "shuffle+blosc2-lz4",  "none", 0, "shuffle", "blosc2",
      "\"blosc2_codec\":\"lz4\",\"blosc2_clevel\":5" },
    { "shuffle+blosc2-zstd", "none", 0, "shuffle", "blosc2",
      "\"blosc2_codec\":\"zstd\",\"blosc2_clevel\":5" },
    /* lossy: simple_packing */
    { "pack24",           "simple_packing", 24, "none", "none", "" },
    { "pack24+lz4",       "simple_packing", 24, "none", "lz4",  "" },
    { "pack24+zstd-3",    "simple_packing", 24, "none", "zstd", "\"zstd_level\":3" },
    { "pack16+szip",      "simple_packing", 16, "none", "szip",
      "\"szip_rsi\":128,\"szip_block_size\":16,\"szip_flags\":8" },
    { "pack16+zstd-3",    "simple_packing", 16, "none", "zstd", "\"zstd_level\":3" },
    { "pack12+zstd-3",    "simple_packing", 12, "none", "zstd", "\"zstd_level\":3" },
    /* lossy: floating-point compressors */
    { "zfp-rate16",       "none", 0, "none", "zfp",
      "\"zfp_mode\":\"fixed_rate\",\"zfp_rate\":16.0,\"zfp_num_values\":%zu" },
    { "zfp-rate8",        "none", 0, "none", "zfp",
      "\"zfp_mode\":\"fixed_rate\",\"zfp_rate\":8.0,\"zfp_num_values\":%zu" },
    { "zfp-tol1e-6",      "none", 0, "none", "zfp",
      "\"zfp_mode\":\"fixed_accuracy\",\"zfp_tolerance\":1e-6,\"zfp_num_values\":%zu" },
    { "sz3-abs1e-6",      "none", 0, "none", "sz3",
      "\"sz3_error_bound_mode\":\"abs\",\"sz3_error_bound\":1e-6,\"sz3_num_values\":%zu" },
    { "sz3-rel1e-4",      "none", 0, "none", "sz3",
      "\"sz3_error_bound_mode\":\"rel\",\"sz3_error_bound\":1e-4,\"sz3_num_values\":%zu" },
};

#define TGM_NUM_CODECS (sizeof(TGM_CODECS) / sizeof(TGM_CODECS[0]))

typedef struct {
    const char *name;
    size_t file_bytes;
    double encode_time_sec;
    double max_abs_error;
    double rmse;
    int ok;
} tgm_bench_result_t;

/* ── JSON builder ── */

static int tgm_build_json(char *buf, int buflen,
                          const tgm_codec_config_t *c,
                          double T, const double *X, size_t num_values)
{
    int pos = 0, remaining = buflen, n;

#define JAPPEND(...)                                             \
    do {                                                         \
        n = snprintf(buf + pos, remaining, __VA_ARGS__);         \
        if (n < 0 || n >= remaining) return -1;                  \
        pos += n; remaining -= n;                                \
    } while (0)

    JAPPEND("{\"version\":1,\"descriptors\":[{\"type\":\"ndarray\","
            "\"ndim\":2,\"shape\":[%d,%d],\"strides\":[%d,1],"
            "\"dtype\":\"float64\",\"byte_order\":\"little\",",
            grid_size, grid_size, grid_size);

    if (strcmp(c->encoding, "simple_packing") == 0) {
        double ref_val;
        int32_t bin_scale;
        tgm_error err = tgm_simple_packing_compute_params(
            X, num_values, (uint32_t)c->bits_per_value, 0, &ref_val, &bin_scale);
        if (err != TGM_ERROR_OK) return -1;
        JAPPEND("\"encoding\":\"simple_packing\","
                "\"reference_value\":%.17g,\"binary_scale_factor\":%d,"
                "\"decimal_scale_factor\":0,\"bits_per_value\":%d,",
                ref_val, (int)bin_scale, c->bits_per_value);
    } else {
        JAPPEND("\"encoding\":\"none\",");
    }

    JAPPEND("\"filter\":\"%s\"", c->filter);
    if (strcmp(c->filter, "shuffle") == 0)
        JAPPEND(",\"shuffle_element_size\":8");
    JAPPEND(",");

    JAPPEND("\"compression\":\"%s\"", c->compression);
    if (c->comp_params[0] != '\0') {
        char expanded[512];
        snprintf(expanded, sizeof(expanded), c->comp_params, num_values);
        JAPPEND(",%s", expanded);
    }

    JAPPEND("}],\"wave\":{\"time\":%.10g,\"dx\":%.10g,\"dt\":%.10g}}",
            T, (double)dx, (double)dt);

#undef JAPPEND
    return pos;
}

/* ── Single codec benchmark ── */

static tgm_bench_result_t tgm_run_codec_bench(
    const tgm_codec_config_t *codec,
    double **snapshots, double *snap_times, int num_snaps)
{
    tgm_bench_result_t res = {0};
    res.name = codec->name;

    size_t num_values = (size_t)grid_size * grid_size;
    size_t raw_bytes  = num_values * sizeof(double);

    char fname[256];
    snprintf(fname, sizeof(fname), "bench_%s.tgm", codec->name);

    tgm_file_t *file = NULL;
    tgm_error err = tgm_file_create(fname, &file);
    if (err != TGM_ERROR_OK) {
        fprintf(stderr, "  [%s] create failed: %s\n", codec->name, tgm_last_error());
        return res;
    }

    double t0 = MPI_Wtime();

    for (int s = 0; s < num_snaps; s++) {
        char json[4096];
        if (tgm_build_json(json, sizeof(json), codec,
                           snap_times[s], snapshots[s], num_values) < 0) {
            fprintf(stderr, "  [%s] JSON overflow at step %d\n", codec->name, s);
            tgm_file_close(file);
            return res;
        }

        const uint8_t *ptrs[1] = { (const uint8_t *)snapshots[s] };
        size_t lens[1] = { raw_bytes };

        err = tgm_file_append(file, json, ptrs, lens, 1, "xxh3", 0);
        if (err != TGM_ERROR_OK) {
            fprintf(stderr, "  [%s] encode+append failed at step %d: %s\n",
                    codec->name, s, tgm_last_error());
            tgm_file_close(file);
            return res;
        }
    }

    res.encode_time_sec = MPI_Wtime() - t0;
    tgm_file_close(file);

    struct stat st;
    if (stat(fname, &st) == 0)
        res.file_bytes = (size_t)st.st_size;

    /* accuracy: decode first frame, compare with original */
    tgm_file_t *rf = NULL;
    err = tgm_file_open(fname, &rf);
    if (err == TGM_ERROR_OK) {
        tgm_message_t *msg = NULL;
        err = tgm_file_decode_message(rf, 0, 1, 1, 0, &msg);
        if (err == TGM_ERROR_OK) {
            size_t dlen = 0;
            const uint8_t *dptr = tgm_object_data(msg, 0, &dlen);
            if (dptr && dlen == raw_bytes) {
                const double *dec = (const double *)dptr;
                double max_err = 0, sum_sq = 0;
                for (size_t i = 0; i < num_values; i++) {
                    double e = fabs(dec[i] - snapshots[0][i]);
                    if (e > max_err) max_err = e;
                    sum_sq += e * e;
                }
                res.max_abs_error = max_err;
                res.rmse = sqrt(sum_sq / num_values);
            }
            tgm_message_free(msg);
        } else {
            fprintf(stderr, "  [%s] decode failed: %s\n", codec->name, tgm_last_error());
        }
        tgm_file_close(rf);
    }

    res.ok = 1;
    return res;
}

/* ── Results printer ── */

static void tgm_print_results(tgm_bench_result_t *results, int n, int num_snaps) {
    size_t num_values = (size_t)grid_size * grid_size;
    size_t raw_total  = num_values * sizeof(double) * num_snaps;

    const char *hdr =
        "| # | Config | Encoding | Filter | Compression | "
        "Size (MB) | Ratio | Encode MB/s | Max Error | RMSE |\n"
        "|---|--------|----------|--------|-------------|"
        "-----------|-------|-------------|-----------|------|\n";

    printf("\n%s", hdr);

    FILE *fp = fopen("benchmark.md", "w");
    if (fp) {
        fprintf(fp, "# Tensogram Compression Benchmark\n\n");
        fprintf(fp, "Grid: %d x %d float64, %d timesteps\n",
                grid_size, grid_size, num_snaps);
        fprintf(fp, "Raw size per step: %.3f MB, total: %.3f MB\n\n",
                (double)(num_values * 8) / 1e6, (double)raw_total / 1e6);
        fprintf(fp, "%s", hdr);
    }

    for (int i = 0; i < n; i++) {
        const tgm_codec_config_t *c = &TGM_CODECS[i];
        if (!results[i].ok) {
            const char *skip = "| %2d | %-20s | - | - | - | SKIP | - | - | - | - |\n";
            printf(skip, i + 1, results[i].name);
            if (fp) fprintf(fp, skip, i + 1, results[i].name);
            continue;
        }

        double size_mb = (double)results[i].file_bytes / 1e6;
        double ratio   = (raw_total > 0) ? (double)raw_total / results[i].file_bytes : 0;
        double enc_mbs = (results[i].encode_time_sec > 0)
            ? ((double)raw_total / 1e6) / results[i].encode_time_sec : 0;

        const char *fmt =
            "| %2d | %-20s | %-15s | %-7s | %-11s | %9.3f | %5.1fx | %11.1f | %.2e | %.2e |\n";
        printf(fmt, i + 1, c->name, c->encoding, c->filter, c->compression,
               size_mb, ratio, enc_mbs, results[i].max_abs_error, results[i].rmse);
        if (fp)
            fprintf(fp, fmt, i + 1, c->name, c->encoding, c->filter, c->compression,
                    size_mb, ratio, enc_mbs, results[i].max_abs_error, results[i].rmse);
    }

    if (fp) {
        fprintf(fp, "\n*Generated by 2d_wave_mpi benchmark*\n");
        fclose(fp);
        printf("\nResults written to benchmark.md\n");
    }
}

/* ── Top-level benchmark driver ── */

static void tgm_run_benchmark(double **snapshots, double *snap_times, int snap_count) {
    printf("\n======================================================\n");
    printf("  Tensogram Compression Benchmark: %d codecs x %d frames\n",
           (int)TGM_NUM_CODECS, snap_count);
    printf("  Grid: %d x %d float64 (%.3f MB/frame)\n",
           grid_size, grid_size, (double)(grid_size * grid_size * 8) / 1e6);
    printf("======================================================\n\n");

    tgm_bench_result_t results[TGM_NUM_CODECS];
    for (int c = 0; c < (int)TGM_NUM_CODECS; c++) {
        printf("  [%2d/%d] %-20s ... ", c + 1, (int)TGM_NUM_CODECS, TGM_CODECS[c].name);
        fflush(stdout);
        results[c] = tgm_run_codec_bench(&TGM_CODECS[c], snapshots, snap_times, snap_count);
        if (results[c].ok) {
            double raw_total = (double)grid_size * grid_size * 8.0 * snap_count;
            double ratio = raw_total /
                           (results[c].file_bytes > 0 ? results[c].file_bytes : 1);
            double mbs = results[c].encode_time_sec > 0
                ? (raw_total / 1e6) / results[c].encode_time_sec : 0;
            printf("%.3f MB, %.1fx, %.1f MB/s\n",
                   (double)results[c].file_bytes / 1e6, ratio, mbs);
        } else {
            printf("FAILED\n");
        }
    }

    tgm_print_results(results, (int)TGM_NUM_CODECS, snap_count);
}

#endif /* write_tensogram */
#endif /* IO_TENSOGRAM_H */
