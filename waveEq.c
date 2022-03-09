#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "math.h"
#include "mpi.h"
#include "netcdf.h"
#include "unistd.h"

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
    MPI_Comm_size(MPI_COMM_WORLD, &(r.MAX));  // get my rank in the new topology
    r.TOP = (int)sqrt(r.MAX);                 // r.TOP-> topology number
    if (r.MAX != r.TOP * r.TOP)
        exit(-1);
    dims[0] = r.TOP;
    dims[1] = r.TOP;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 1, &(r.top_world));
    MPI_Comm_rank(r.top_world, &(r.ME));              // get my rank in the new topology
    MPI_Cart_coords(r.top_world, r.ME, 2, mycoords);  // get my coordinates
    MPI_Cart_shift(r.top_world, 0, 1, &(r.U), &(r.D));
    MPI_Cart_shift(r.top_world, 1, 1, &(r.L), &(r.R));
    return r;
}

#ifdef write_netcdf
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

#endif


void apply_initial_cond(double* vec) {
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            vec[i * grid_size + j] = i * dx * ((f_to_go)-i * (dx)) * j * (dx) * ((f_to_go)-j * (dx));
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
            x_nxt = (b[n * j + i] - x_nxt) / on_diag;
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

    if (write_netcdf) {
        static ncFileProps nc_file;
        if (access(ncFileName, F_OK)) {
            nc_file = create_ncFile();
        }


        for (int i = 0; i < grid_size; ++i) {
            nc_file.start[1] = i;
            CHECK(nc_put_vara_double(nc_file.ncFileID, nc_file.ncWaveLenVarID, nc_file.start, nc_file.count, &X[i * grid_size]));
        }
        nc_file.start[0] += 1;
        if (fabs(T - (double)t_to_go) < 2e-5) {
            printf("NetCDF Closed!\n");
            nc_close(nc_file.ncFileID);
        }
    }
}

void apply_boundary(ranks_t rank, double* local_x, int local_size) {
    if (rank.L == -2) {
        for (int i = 0; i < local_size; i++) {
            local_x[i * local_size] = 0.;
        }
    }
    if (rank.R == -2) {
        for (int i = 0; i < local_size; i++) {
            local_x[i * local_size + local_size - 1] = 0.;
        }
    }
    if (rank.U == -2) {
        for (int i = 0; i < local_size; i++) {
            local_x[i] = 0.;
        }
    }
    if (rank.D == -2) {
        for (int i = 0; i < local_size; i++) {
            local_x[local_size * (local_size - 1) + i] = 0.;
        }
    }
}


int main(int argc, char* argv[]) {
    ranks_t rank;
    MPI_Status stat;
    MPI_Request req;
    MPI_Init(&argc, &argv);

    rank = get_ranks();
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
        calcB(r1, r2, b);
        for (int i = 0; i < rank.TOP; i++) {
            for (int j = 0; j < rank.TOP; j++) {
                MPI_Isend(&b[i * (grid_size * (local_size - 2)) + j * (local_size - 2)],
                          1, subgrid, i * rank.TOP + j, i * rank.TOP + j,
                          rank.top_world, &req);
            }
        }
    }
    MPI_Recv(local_b, (local_size - 2) * (local_size - 2), MPI_DOUBLE, 0, rank.ME,
             rank.top_world, &stat);
    local_x = (double*)malloc(local_size * local_size * sizeof(double));
    Fill_X_Init(local_x, local_size);
    apply_boundary(rank, local_x, local_size);
    for (t = 0; t < t_to_go; t += dt) {
        do {
            err = GaussSeidelRB(local_x, local_b, local_size - 2);
            MPI_Isend(&local_x[local_size], 1, vt_horizontal, rank.U, rank.U + 10000,
                      rank.top_world, &req);
            MPI_Recv(&local_x[local_size * (local_size - 1)], 1, vt_horizontal,
                     rank.D, rank.ME + 10000, rank.top_world, &stat);
            MPI_Isend(&local_x[local_size * (local_size - 2)], 1, vt_horizontal,
                      rank.D, rank.D + 10000, rank.top_world, &req);
            MPI_Recv(local_x, 1, vt_horizontal, rank.U, rank.ME + 10000,
                     rank.top_world, &stat);
            MPI_Isend(&local_x[1], 1, vt_vertical, rank.L, rank.L + 10000,
                      rank.top_world, &req);
            MPI_Recv(&local_x[local_size - 1], 1, vt_vertical, rank.R,
                     rank.ME + 10000, rank.top_world, &stat);
            MPI_Isend(&local_x[local_size - 2], 1, vt_vertical, rank.R,
                      rank.R + 10000, rank.top_world, &req);
            MPI_Recv(local_x, 1, vt_vertical, rank.L, rank.ME + 10000, rank.top_world,
                     &stat);
            MPI_Allreduce(&err, &gerr, 1, MPI_DOUBLE, MPI_MAX, rank.top_world);
            if (rank.ME == 0) {
                printf("%e\t%f\t%d\n", gerr, t, grid_size);
            }
        } while (gerr > gerr_max);
        apply_boundary(rank, local_x, local_size);
        MPI_Isend(&local_x[local_size + 1], 1, nonghost, 0, rank.ME, rank.top_world,
                  &req);


        if (rank.ME == 0) {
            vec_copy(r1, r2);

            export_to_file(t, r1, grid_size);

            for (int i = 0; i < rank.TOP; i++) {
                for (int j = 0; j < rank.TOP; j++) {
                    MPI_Recv(&r1[(i * (grid_size) * (grid_size / rank.TOP)) + (j * (grid_size / rank.TOP))],
                             1, subgrid, i * rank.TOP + j, i * rank.TOP + j,
                             rank.top_world, &stat);
                }
            }
            calcB(r1, r2, b);
            for (int i = 0; i < rank.TOP; i++) {
                for (int j = 0; j < rank.TOP; j++) {
                    MPI_Isend(&b[(i * (grid_size) * (grid_size / rank.TOP)) + (j * (grid_size / rank.TOP))],
                              1, subgrid, i * rank.TOP + j, i * rank.TOP + j,
                              rank.top_world, &req);
                }
            }
        }
        MPI_Recv(local_b, (local_size - 2) * (local_size - 2), MPI_DOUBLE, 0,
                 rank.ME, rank.top_world, &stat);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(local_b);
    free(local_x);
    if (rank.ME == 0) {
        free(r1);
        free(r2);
        free(b);
    }
    MPI_Finalize();
    return 0;
}
