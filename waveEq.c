#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "config.h"
#if write_netcdf
#include <netcdf.h>
#endif
#if write_tensogram
#include "tensogram.h"
#endif

typedef struct {
    // u->up ..
    int U, D, L, R, ME, MAX, TOP;
    MPI_Comm top_world;
} ranks_t;

ranks_t get_ranks(void) {
    ranks_t r;  // r-> return value
    int dims[2], periodic[2], mycoords[2];
    periodic[0] = 0;
    periodic[1] = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &(r.MAX));
    r.TOP = (int)sqrt(r.MAX);
    if (r.MAX != r.TOP * r.TOP) {
        fprintf(stderr, "Error: number of ranks (%d) must be a perfect square\n", r.MAX);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (grid_size % r.TOP != 0) {
        fprintf(stderr, "Error: grid_size (%d) must be divisible by sqrt(ranks) (%d)\n", grid_size, r.TOP);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    dims[0] = r.TOP;
    dims[1] = r.TOP;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &(r.top_world));
    MPI_Comm_rank(r.top_world, &(r.ME));              // get my rank in the new topology
    MPI_Cart_coords(r.top_world, r.ME, 2, mycoords);  // get my coordinates
    MPI_Cart_shift(r.top_world, 0, 1, &(r.U), &(r.D));
    MPI_Cart_shift(r.top_world, 1, 1, &(r.L), &(r.R));
    return r;
}

#if write_netcdf
typedef struct ncFileProps {
    NETCDFID ncDimIDs[3];
    NETCDFID ncFileID;
    NETCDFID ncXDimID, ncYDimID, ncTDimID;
    NETCDFID ncXVarID, ncYVarID, ncTVarID;
    NETCDFID ncWaveLenVarID;
    size_t start[3];
    size_t count[3];
} ncFileProps;


ncFileProps create_ncFile() {
    ncFileProps nc_file;
    nc_file.start[0] = 0;
    nc_file.start[1] = 0;
    nc_file.start[2] = 0;
    nc_file.count[0] = 1;
    nc_file.count[1] = 1;
    nc_file.count[2] = grid_size;
    CHECK(nc_create(ncFileName, NC_CLOBBER, &(nc_file.ncFileID)));
    CHECK(nc_def_dim(nc_file.ncFileID, "time", NC_UNLIMITED, &(nc_file.ncTDimID)));
    CHECK(nc_def_dim(nc_file.ncFileID, "X Dimension", grid_size, &(nc_file.ncXDimID)));
    CHECK(nc_def_dim(nc_file.ncFileID, "Y Dimension", grid_size, &(nc_file.ncYDimID)));

    CHECK(nc_def_var(nc_file.ncFileID, "time", NC_FLOAT, 1, &(nc_file.ncTDimID), &(nc_file.ncTVarID)));
    CHECK(nc_def_var(nc_file.ncFileID, "X Variable", NC_FLOAT, 1, &(nc_file.ncXDimID), &(nc_file.ncXVarID)));
    CHECK(nc_def_var(nc_file.ncFileID, "Y Variable", NC_FLOAT, 1, &(nc_file.ncYDimID), &(nc_file.ncYVarID)));

    nc_file.ncDimIDs[0] = nc_file.ncTDimID;
    nc_file.ncDimIDs[1] = nc_file.ncYDimID;
    nc_file.ncDimIDs[2] = nc_file.ncXDimID;

    CHECK(nc_def_var(nc_file.ncFileID, "WAVE LENGTH", NC_FLOAT, 3, nc_file.ncDimIDs, &(nc_file.ncWaveLenVarID)));

    CHECK(nc_enddef(nc_file.ncFileID));
    // assigning dimension variables
    float grid_x_points[grid_size], grid_y_points[grid_size];

    for (int i = 0; i < grid_size; ++i) {
        grid_x_points[i] = (-(grid_size / 2) + i) * dx;
        grid_y_points[i] = (-(grid_size / 2) + i) * dx;
    }

    CHECK(nc_put_var_float(nc_file.ncFileID, nc_file.ncXVarID, grid_x_points));
    CHECK(nc_put_var_float(nc_file.ncFileID, nc_file.ncYVarID, grid_y_points));
    return nc_file;
}

static ncFileProps g_nc_file;

void nc_write_step(double T, double* X) {
    float time_val = (float)T;
    size_t time_idx = g_nc_file.start[0];
    CHECK(nc_put_var1_float(g_nc_file.ncFileID, g_nc_file.ncTVarID, &time_idx, &time_val));

    for (int i = 0; i < grid_size; ++i) {
        g_nc_file.start[1] = i;
        CHECK(nc_put_vara_double(g_nc_file.ncFileID, g_nc_file.ncWaveLenVarID,
                                 g_nc_file.start, g_nc_file.count, &X[i * grid_size]));
    }
    g_nc_file.start[0] += 1;
    CHECK(nc_sync(g_nc_file.ncFileID));
}

void nc_close_file(void) {
    CHECK(nc_close(g_nc_file.ncFileID));
}

#endif

#if write_tensogram
static struct tgm_TgmFile *g_tgm_file = NULL;

void tgm_create_output(void) {
    TGM_CHECK(tgm_file_create(tgmFileName, &g_tgm_file));
}

void tgm_write_step(double T, double* X) {
    size_t num_values = (size_t)grid_size * grid_size;
    size_t num_bytes = num_values * sizeof(double);
    char json[2048];
    int remaining = (int)sizeof(json);
    int pos = 0;
    int n;

#define JSON_APPEND(...)                                        \
    do {                                                        \
        n = snprintf(json + pos, remaining, __VA_ARGS__);       \
        if (n < 0 || n >= remaining) {                          \
            fprintf(stderr, "Tensogram: JSON buffer overflow\n");\
            MPI_Abort(MPI_COMM_WORLD, 1);                       \
        }                                                       \
        pos += n;                                               \
        remaining -= n;                                         \
    } while (0)

    JSON_APPEND("{\"version\":1,"
        "\"objects\":[{\"type\":\"ntensor\",\"ndim\":2,"
            "\"shape\":[%d,%d],\"strides\":[%d,1],\"dtype\":\"float64\"}],"
        "\"payload\":[{\"byte_order\":\"little\",",
        grid_size, grid_size, grid_size);

#if tgm_use_packing
    double ref_val;
    int32_t bin_scale;
    TGM_CHECK(tgm_simple_packing_compute_params(
        X, num_values, tgm_bits_per_value, 0, &ref_val, &bin_scale));

    JSON_APPEND("\"encoding\":\"simple_packing\","
        "\"reference_value\":%.17g,\"binary_scale_factor\":%d,"
        "\"decimal_scale_factor\":0,\"bits_per_value\":%d,",
        ref_val, (int)bin_scale, tgm_bits_per_value);
#else
    JSON_APPEND("\"encoding\":\"none\",");
#endif

    JSON_APPEND("\"filter\":\"none\",");

#if tgm_use_szip
    JSON_APPEND("\"compression\":\"szip\","
        "\"szip_rsi\":%d,\"szip_block_size\":%d,\"szip_flags\":%d",
        tgm_szip_rsi, tgm_szip_block_size, tgm_szip_flags);
#else
    JSON_APPEND("\"compression\":\"none\"");
#endif

    JSON_APPEND("}],\"wave\":{\"time\":%.10g,\"dx\":%.10g,\"dt\":%.10g}}",
        T, (double)dx, (double)dt);

#undef JSON_APPEND

    const uint8_t *data_ptrs[1] = { (const uint8_t*)X };
    uintptr_t data_lens[1] = { num_bytes };
    struct tgm_TgmBytes encoded;

    TGM_CHECK(tgm_encode(json, data_ptrs, data_lens, 1, NULL, &encoded));
    TGM_CHECK(tgm_file_append_raw(g_tgm_file, encoded.data, encoded.len));
    tgm_bytes_free(encoded);
}

void tgm_close_output(void) {
    if (g_tgm_file != NULL) {
        tgm_file_close(g_tgm_file);
        g_tgm_file = NULL;
    }
}
#endif


void apply_initial_cond(double* vec) {
    double cx1 = f_to_go * 0.3;
    double cy1 = f_to_go * 0.3;
    double cx2 = f_to_go * 0.7;
    double cy2 = f_to_go * 0.7;
    double sigma = f_to_go / 15.0;
    double sigma2 = sigma * sigma;
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            double x = i * dx;
            double y = j * dx;
            double r2a = (x - cx1) * (x - cx1) + (y - cy1) * (y - cy1);
            double r2b = (x - cx2) * (x - cx2) + (y - cy2) * (y - cy2);
            vec[i * grid_size + j] = exp(-r2a / (2.0 * sigma2)) + exp(-r2b / (2.0 * sigma2));
        }
    }
}

void vec_copy(double* vec1, double* vec2) {
    for (int i = 0; i < grid_size * grid_size; i++) {
        vec2[i] = vec1[i];
    }
}

void calcB(double* R1, double* R2, double* b) {
    for (int i = 0; i < grid_size * grid_size; i++) {
        b[i] = 2 * R1[i] - R2[i];
    }
}

void Fill_X_Init(double* mat, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            mat[i * size + j] = ((double)rand()) / RAND_MAX;
        }
    }
}

double GaussSeidelRB(double* x, double* b, int n) {
    double err;
    int ib;
    double x_nxt;
    int up, dn, lt, rt;
    up  = -(n + 2);
    dn  = n + 2;
    lt  = -1;
    rt  = 1;
    err = 0.;
    for (int j = 0; j < n; j++) {
        ib = (n + 2) * (j + 1) + 1;
        for (int i = 0; i < n; i++) {
            x_nxt = off_diag * (x[ib + i + up] + x[ib + i + lt] + x[ib + i + rt] + x[ib + i + dn]);
            x_nxt = (b[n * j + i] + x_nxt) / on_diag;
            err += fabs(x[ib + i] - x_nxt);
            x[ib + i] = x_nxt;
        }
    }
    return err / ((double)(n * n));
}


void export_to_file(double T, double* X, int size) {
    if (write_posix) {
        FILE* fp;
        char fn[128];
        char ff[] = "t_%04d.txt";  // filename format (id_me, T)
        char nf[] = "%6.6lf ";     // number format in file

        sprintf(fn, ff, (int)(T / dt));
        fp = fopen(fn, "w");
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                fprintf(fp, nf, X[(size)*j + i]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
        fclose(fp);
    }

#if write_netcdf
    nc_write_step(T, X);
#endif
#if write_tensogram
    tgm_write_step(T, X);
#endif
}

void apply_boundary(ranks_t rank, double* local_x, int local_size) {
    if (rank.L == MPI_PROC_NULL) {
        for (int i = 0; i < local_size; i++) {
            local_x[i * local_size] = 0.;
        }
    }
    if (rank.R == MPI_PROC_NULL) {
        for (int i = 0; i < local_size; i++) {
            local_x[i * local_size + local_size - 1] = 0.;
        }
    }
    if (rank.U == MPI_PROC_NULL) {
        for (int i = 0; i < local_size; i++) {
            local_x[i] = 0.;
        }
    }
    if (rank.D == MPI_PROC_NULL) {
        for (int i = 0; i < local_size; i++) {
            local_x[local_size * (local_size - 1) + i] = 0.;
        }
    }
}


int main(int argc, char* argv[]) {
    ranks_t rank;
    MPI_Status stat;
    MPI_Request req;
    MPI_Request *scatter_reqs = NULL;
    MPI_Init(&argc, &argv);
    setbuf(stdout, NULL);

    rank = get_ranks();
    scatter_reqs = (MPI_Request*)malloc(rank.MAX * sizeof(MPI_Request));
    double *r1, *r2, *b, *local_x, *local_b, t, gerr, err;
    int local_size;
    local_size = (grid_size / rank.TOP) + 2;
    MPI_Datatype subgrid, vt_vertical, vt_horizontal, nonghost;
    MPI_Type_vector(local_size - 2, local_size - 2, grid_size, MPI_DOUBLE,
                    &subgrid);
    MPI_Type_vector(local_size, 1, local_size, MPI_DOUBLE, &vt_vertical);
    MPI_Type_vector(1, local_size, local_size, MPI_DOUBLE, &vt_horizontal);
    MPI_Type_vector(local_size - 2, local_size - 2, local_size, MPI_DOUBLE,
                    &nonghost);
    MPI_Type_commit(&subgrid);
    MPI_Type_commit(&vt_vertical);
    MPI_Type_commit(&vt_horizontal);
    MPI_Type_commit(&nonghost);
    local_b = (double*)malloc((local_size - 2) * (local_size - 2) * sizeof(double));
    if (rank.ME == 0) {
        r1 = (double*)malloc(grid_size * grid_size * sizeof(double));
        apply_initial_cond(r1);
        r2 = (double*)malloc(grid_size * grid_size * sizeof(double));
        vec_copy(r1, r2);
        b = (double*)malloc(grid_size * grid_size * sizeof(double));
#if write_netcdf
        g_nc_file = create_ncFile();
#endif
#if write_tensogram
        tgm_create_output();
#endif
        calcB(r1, r2, b);
        for (int i = 0; i < rank.TOP; i++) {
            for (int j = 0; j < rank.TOP; j++) {
                int dest = i * rank.TOP + j;
                MPI_Isend(&b[i * (grid_size * (local_size - 2)) + j * (local_size - 2)],
                          1, subgrid, dest, dest,
                          rank.top_world, &scatter_reqs[dest]);
            }
        }
    }
    MPI_Recv(local_b, (local_size - 2) * (local_size - 2), MPI_DOUBLE, 0, rank.ME,
             rank.top_world, &stat);
    if (rank.ME == 0) {
        MPI_Waitall(rank.MAX, scatter_reqs, MPI_STATUSES_IGNORE);
    }
    local_x = (double*)malloc(local_size * local_size * sizeof(double));
    Fill_X_Init(local_x, local_size);
    apply_boundary(rank, local_x, local_size);
    for (t = 0; t < t_to_go; t += dt) {
        do {
            err = GaussSeidelRB(local_x, local_b, local_size - 2);
            MPI_Sendrecv(&local_x[local_size], 1, vt_horizontal, rank.U, 0,
                         &local_x[local_size * (local_size - 1)], 1, vt_horizontal,
                         rank.D, 0, rank.top_world, &stat);
            MPI_Sendrecv(&local_x[local_size * (local_size - 2)], 1, vt_horizontal, rank.D, 1,
                         local_x, 1, vt_horizontal,
                         rank.U, 1, rank.top_world, &stat);
            MPI_Sendrecv(&local_x[1], 1, vt_vertical, rank.L, 2,
                         &local_x[local_size - 1], 1, vt_vertical,
                         rank.R, 2, rank.top_world, &stat);
            MPI_Sendrecv(&local_x[local_size - 2], 1, vt_vertical, rank.R, 3,
                         local_x, 1, vt_vertical,
                         rank.L, 3, rank.top_world, &stat);
            MPI_Allreduce(&err, &gerr, 1, MPI_DOUBLE, MPI_MAX, rank.top_world);
            if (rank.ME == 0) {
                printf("%e\t%f\t%d\n", gerr, t, grid_size);
            }
        } while (gerr > gerr_max);
        apply_boundary(rank, local_x, local_size);
        MPI_Isend(&local_x[local_size + 1], 1, nonghost, 0, rank.ME, rank.top_world, &req);

        if (rank.ME == 0) {
            vec_copy(r1, r2);

            for (int i = 0; i < rank.TOP; i++) {
                for (int j = 0; j < rank.TOP; j++) {
                    MPI_Recv(&r1[(i * (grid_size) * (grid_size / rank.TOP)) + (j * (grid_size / rank.TOP))],
                             1, subgrid, i * rank.TOP + j, i * rank.TOP + j,
                             rank.top_world, &stat);
                }
            }
            export_to_file(t, r1, grid_size);
            calcB(r1, r2, b);
            for (int i = 0; i < rank.TOP; i++) {
                for (int j = 0; j < rank.TOP; j++) {
                    int dest = i * rank.TOP + j;
                    MPI_Isend(&b[(i * (grid_size) * (grid_size / rank.TOP)) + (j * (grid_size / rank.TOP))],
                              1, subgrid, dest, dest,
                              rank.top_world, &scatter_reqs[dest]);
                }
            }
        }
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        MPI_Recv(local_b, (local_size - 2) * (local_size - 2), MPI_DOUBLE, 0,
                 rank.ME, rank.top_world, &stat);
        if (rank.ME == 0) {
            MPI_Waitall(rank.MAX, scatter_reqs, MPI_STATUSES_IGNORE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(local_b);
    free(local_x);
    free(scatter_reqs);
    MPI_Type_free(&subgrid);
    MPI_Type_free(&vt_vertical);
    MPI_Type_free(&vt_horizontal);
    MPI_Type_free(&nonghost);
    MPI_Comm_free(&rank.top_world);
    if (rank.ME == 0) {
#if write_netcdf
        nc_close_file();
#endif
#if write_tensogram
        tgm_close_output();
#endif
        free(r1);
        free(r2);
        free(b);
    }
    MPI_Finalize();
    return 0;
}
