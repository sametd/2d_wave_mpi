#ifndef IO_NETCDF_H
#define IO_NETCDF_H

#include "config.h"

#if write_netcdf

#include <stdio.h>
#include <netcdf.h>
#include <mpi.h>

#define ncFileName "grid.nc"

typedef int NETCDFID;

#define NC_CHECK(e)                                                 \
    {                                                               \
        int err = e;                                                \
        if (err) {                                                  \
            fprintf(stderr, "NetCDF error: %s\n", nc_strerror(err));\
            MPI_Abort(MPI_COMM_WORLD, err);                         \
        }                                                           \
    }

typedef struct {
    NETCDFID ncDimIDs[3];
    NETCDFID ncFileID;
    NETCDFID ncXDimID, ncYDimID, ncTDimID;
    NETCDFID ncXVarID, ncYVarID, ncTVarID;
    NETCDFID ncWaveLenVarID;
    size_t start[3];
    size_t count[3];
} ncFileProps;

static ncFileProps g_nc_file;

static ncFileProps nc_create_file(void) {
    ncFileProps nc;
    nc.start[0] = 0; nc.start[1] = 0; nc.start[2] = 0;
    nc.count[0] = 1; nc.count[1] = 1; nc.count[2] = grid_size;

    NC_CHECK(nc_create(ncFileName, NC_CLOBBER, &nc.ncFileID));
    NC_CHECK(nc_def_dim(nc.ncFileID, "time", NC_UNLIMITED, &nc.ncTDimID));
    NC_CHECK(nc_def_dim(nc.ncFileID, "X Dimension", grid_size, &nc.ncXDimID));
    NC_CHECK(nc_def_dim(nc.ncFileID, "Y Dimension", grid_size, &nc.ncYDimID));

    NC_CHECK(nc_def_var(nc.ncFileID, "time", NC_FLOAT, 1, &nc.ncTDimID, &nc.ncTVarID));
    NC_CHECK(nc_def_var(nc.ncFileID, "X Variable", NC_FLOAT, 1, &nc.ncXDimID, &nc.ncXVarID));
    NC_CHECK(nc_def_var(nc.ncFileID, "Y Variable", NC_FLOAT, 1, &nc.ncYDimID, &nc.ncYVarID));

    nc.ncDimIDs[0] = nc.ncTDimID;
    nc.ncDimIDs[1] = nc.ncYDimID;
    nc.ncDimIDs[2] = nc.ncXDimID;
    NC_CHECK(nc_def_var(nc.ncFileID, "WAVE LENGTH", NC_FLOAT, 3, nc.ncDimIDs, &nc.ncWaveLenVarID));
    NC_CHECK(nc_enddef(nc.ncFileID));

    float grid_x[grid_size], grid_y[grid_size];
    for (int i = 0; i < grid_size; ++i) {
        grid_x[i] = (-(grid_size / 2) + i) * dx;
        grid_y[i] = (-(grid_size / 2) + i) * dx;
    }
    NC_CHECK(nc_put_var_float(nc.ncFileID, nc.ncXVarID, grid_x));
    NC_CHECK(nc_put_var_float(nc.ncFileID, nc.ncYVarID, grid_y));
    return nc;
}

static void nc_write_step(double T, double *X) {
    float time_val = (float)T;
    size_t time_idx = g_nc_file.start[0];
    NC_CHECK(nc_put_var1_float(g_nc_file.ncFileID, g_nc_file.ncTVarID, &time_idx, &time_val));
    for (int i = 0; i < grid_size; ++i) {
        g_nc_file.start[1] = i;
        NC_CHECK(nc_put_vara_double(g_nc_file.ncFileID, g_nc_file.ncWaveLenVarID,
                                    g_nc_file.start, g_nc_file.count, &X[i * grid_size]));
    }
    g_nc_file.start[0] += 1;
    NC_CHECK(nc_sync(g_nc_file.ncFileID));
}

static void nc_close_file(void) {
    NC_CHECK(nc_close(g_nc_file.ncFileID));
}

#endif /* write_netcdf */
#endif /* IO_NETCDF_H */
