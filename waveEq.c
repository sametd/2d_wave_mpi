/*
 * 2D Wave Equation Solver
 *
 * Solves ddu/dtt = laplacian(u) on a square domain with zero boundary
 * conditions using implicit time stepping with Gauss-Seidel iteration
 * and 2D MPI domain decomposition.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "config.h"
#include "io_netcdf.h"
#include "io_tensogram.h"

/* ── MPI topology ── */

typedef struct {
    int U, D, L, R, ME, MAX, TOP;
    MPI_Comm top_world;
} ranks_t;

static ranks_t get_ranks(void) {
    ranks_t r;
    int dims[2] = {0, 0}, periodic[2] = {0, 0};
    MPI_Comm_size(MPI_COMM_WORLD, &r.MAX);
    r.TOP = (int)sqrt(r.MAX);
    if (r.MAX != r.TOP * r.TOP) {
        fprintf(stderr, "Error: number of ranks (%d) must be a perfect square\n", r.MAX);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (grid_size % r.TOP != 0) {
        fprintf(stderr, "Error: grid_size (%d) must be divisible by sqrt(ranks) (%d)\n",
                grid_size, r.TOP);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    dims[0] = r.TOP;
    dims[1] = r.TOP;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &r.top_world);
    MPI_Comm_rank(r.top_world, &r.ME);
    MPI_Cart_shift(r.top_world, 0, 1, &r.U, &r.D);
    MPI_Cart_shift(r.top_world, 1, 1, &r.L, &r.R);
    return r;
}

/* ── Solver kernels ── */

static void apply_initial_cond(double *vec) {
    double cx1 = f_to_go * 0.3, cy1 = f_to_go * 0.3;
    double cx2 = f_to_go * 0.7, cy2 = f_to_go * 0.7;
    double sigma = f_to_go / 15.0, sigma2 = sigma * sigma;
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            double x = i * dx, y = j * dx;
            double r2a = (x - cx1) * (x - cx1) + (y - cy1) * (y - cy1);
            double r2b = (x - cx2) * (x - cx2) + (y - cy2) * (y - cy2);
            vec[i * grid_size + j] = exp(-r2a / (2.0 * sigma2))
                                   + exp(-r2b / (2.0 * sigma2));
        }
    }
}

static void vec_copy(const double *src, double *dst) {
    memcpy(dst, src, (size_t)grid_size * grid_size * sizeof(double));
}

static void calcB(const double *R1, const double *R2, double *b) {
    for (int i = 0; i < grid_size * grid_size; i++)
        b[i] = 2 * R1[i] - R2[i];
}

static void fill_random(double *mat, int size) {
    for (int i = 0; i < size * size; i++)
        mat[i] = ((double)rand()) / RAND_MAX;
}

static double gauss_seidel(double *x, const double *b, int n) {
    int up = -(n + 2), dn = n + 2, lt = -1, rt = 1;
    double err = 0.0;
    for (int j = 0; j < n; j++) {
        int ib = (n + 2) * (j + 1) + 1;
        for (int i = 0; i < n; i++) {
            double x_nxt = off_diag * (x[ib+i+up] + x[ib+i+lt] + x[ib+i+rt] + x[ib+i+dn]);
            x_nxt = (b[n * j + i] + x_nxt) / on_diag;
            err += fabs(x[ib + i] - x_nxt);
            x[ib + i] = x_nxt;
        }
    }
    return err / ((double)(n * n));
}

static void apply_boundary(ranks_t rank, double *local_x, int ls) {
    if (rank.L == MPI_PROC_NULL)
        for (int i = 0; i < ls; i++) local_x[i * ls] = 0.0;
    if (rank.R == MPI_PROC_NULL)
        for (int i = 0; i < ls; i++) local_x[i * ls + ls - 1] = 0.0;
    if (rank.U == MPI_PROC_NULL)
        for (int i = 0; i < ls; i++) local_x[i] = 0.0;
    if (rank.D == MPI_PROC_NULL)
        for (int i = 0; i < ls; i++) local_x[ls * (ls - 1) + i] = 0.0;
}

/* ── File output ── */

static void export_to_file(double T, double *X) {
    if (write_posix) {
        char fn[128];
        sprintf(fn, "t_%04d.txt", (int)(T / dt));
        FILE *fp = fopen(fn, "w");
        for (int j = 0; j < grid_size; j++) {
            for (int i = 0; i < grid_size; i++)
                fprintf(fp, "%6.6lf ", X[grid_size * j + i]);
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
#if write_netcdf
    nc_write_step(T, X);
#endif
}

/* ── Main ── */

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    setbuf(stdout, NULL);

    ranks_t rank = get_ranks();
    MPI_Status stat;
    MPI_Request req;
    MPI_Request *scatter_reqs = malloc(rank.MAX * sizeof(MPI_Request));

    int local_size = (grid_size / rank.TOP) + 2;

    MPI_Datatype subgrid, vt_vertical, vt_horizontal, nonghost;
    MPI_Type_vector(local_size - 2, local_size - 2, grid_size, MPI_DOUBLE, &subgrid);
    MPI_Type_vector(local_size, 1, local_size, MPI_DOUBLE, &vt_vertical);
    MPI_Type_vector(1, local_size, local_size, MPI_DOUBLE, &vt_horizontal);
    MPI_Type_vector(local_size - 2, local_size - 2, local_size, MPI_DOUBLE, &nonghost);
    MPI_Type_commit(&subgrid);
    MPI_Type_commit(&vt_vertical);
    MPI_Type_commit(&vt_horizontal);
    MPI_Type_commit(&nonghost);

    double *local_b = malloc((local_size - 2) * (local_size - 2) * sizeof(double));
    double *r1 = NULL, *r2 = NULL, *b = NULL;

#if write_tensogram
    int total_output_steps = (int)(t_to_go / dt) / write_interval + 1;
    int max_snaps = tgm_bench_steps > 0 ? tgm_bench_steps : total_output_steps;
    double **snapshots = NULL;
    double *snap_times = NULL;
    int snap_count = 0;
#endif

    if (rank.ME == 0) {
        r1 = malloc(grid_size * grid_size * sizeof(double));
        r2 = malloc(grid_size * grid_size * sizeof(double));
        b  = malloc(grid_size * grid_size * sizeof(double));
        apply_initial_cond(r1);
        vec_copy(r1, r2);
#if write_netcdf
        g_nc_file = nc_create_file();
#endif
#if write_tensogram
        snapshots  = malloc(max_snaps * sizeof(double *));
        snap_times = malloc(max_snaps * sizeof(double));
#endif
        calcB(r1, r2, b);
        for (int i = 0; i < rank.TOP; i++)
            for (int j = 0; j < rank.TOP; j++) {
                int dest = i * rank.TOP + j;
                MPI_Isend(&b[i * grid_size * (local_size - 2) + j * (local_size - 2)],
                          1, subgrid, dest, dest, rank.top_world, &scatter_reqs[dest]);
            }
    }

    MPI_Recv(local_b, (local_size - 2) * (local_size - 2), MPI_DOUBLE,
             0, rank.ME, rank.top_world, &stat);
    if (rank.ME == 0)
        MPI_Waitall(rank.MAX, scatter_reqs, MPI_STATUSES_IGNORE);

    double *local_x = malloc(local_size * local_size * sizeof(double));
    fill_random(local_x, local_size);
    apply_boundary(rank, local_x, local_size);

    /* ── Time loop ── */
    int step = 0;
    for (double t = 0; t < t_to_go; t += dt, step++) {
        double gerr, err;
        do {
            err = gauss_seidel(local_x, local_b, local_size - 2);
            MPI_Sendrecv(&local_x[local_size], 1, vt_horizontal, rank.U, 0,
                         &local_x[local_size * (local_size - 1)], 1, vt_horizontal,
                         rank.D, 0, rank.top_world, &stat);
            MPI_Sendrecv(&local_x[local_size * (local_size - 2)], 1, vt_horizontal, rank.D, 1,
                         local_x, 1, vt_horizontal, rank.U, 1, rank.top_world, &stat);
            MPI_Sendrecv(&local_x[1], 1, vt_vertical, rank.L, 2,
                         &local_x[local_size - 1], 1, vt_vertical,
                         rank.R, 2, rank.top_world, &stat);
            MPI_Sendrecv(&local_x[local_size - 2], 1, vt_vertical, rank.R, 3,
                         local_x, 1, vt_vertical, rank.L, 3, rank.top_world, &stat);
            MPI_Allreduce(&err, &gerr, 1, MPI_DOUBLE, MPI_MAX, rank.top_world);
            if (rank.ME == 0)
                printf("%e\t%f\t%d\n", gerr, t, grid_size);
        } while (gerr > gerr_max);

        apply_boundary(rank, local_x, local_size);
        MPI_Isend(&local_x[local_size + 1], 1, nonghost, 0, rank.ME, rank.top_world, &req);

        if (rank.ME == 0) {
            vec_copy(r1, r2);
            for (int i = 0; i < rank.TOP; i++)
                for (int j = 0; j < rank.TOP; j++)
                    MPI_Recv(&r1[i * grid_size * (grid_size / rank.TOP) + j * (grid_size / rank.TOP)],
                             1, subgrid, i * rank.TOP + j, i * rank.TOP + j,
                             rank.top_world, &stat);

            if (step % write_interval == 0) {
                export_to_file(t, r1);
#if write_tensogram
                if (snap_count < max_snaps) {
                    size_t nbytes = (size_t)grid_size * grid_size * sizeof(double);
                    snapshots[snap_count] = malloc(nbytes);
                    memcpy(snapshots[snap_count], r1, nbytes);
                    snap_times[snap_count] = t;
                    snap_count++;
                }
#endif
            }
            calcB(r1, r2, b);
            for (int i = 0; i < rank.TOP; i++)
                for (int j = 0; j < rank.TOP; j++) {
                    int dest = i * rank.TOP + j;
                    MPI_Isend(&b[i * grid_size * (grid_size / rank.TOP) + j * (grid_size / rank.TOP)],
                              1, subgrid, dest, dest, rank.top_world, &scatter_reqs[dest]);
                }
        }

        MPI_Wait(&req, MPI_STATUS_IGNORE);
        MPI_Recv(local_b, (local_size - 2) * (local_size - 2), MPI_DOUBLE,
                 0, rank.ME, rank.top_world, &stat);
        if (rank.ME == 0)
            MPI_Waitall(rank.MAX, scatter_reqs, MPI_STATUSES_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);

#if write_tensogram
    if (rank.ME == 0 && snap_count > 0)
        tgm_run_benchmark(snapshots, snap_times, snap_count);

    if (rank.ME == 0) {
        for (int i = 0; i < snap_count; i++) free(snapshots[i]);
        free(snapshots);
        free(snap_times);
    }
#endif

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
        free(r1);
        free(r2);
        free(b);
    }

    MPI_Finalize();
    return 0;
}
