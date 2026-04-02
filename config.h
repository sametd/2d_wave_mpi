#define dx 0.00594998
#define dt 0.0001
// r² = dt²/dx², the CFL-like ratio for the implicit wave equation
#define off_diag ((dt * dt) / (dx * dx))
// diagonal: 1 + 4r²
#define on_diag (1.0 + 4.0 * off_diag)
#define f_to_go (3.14159)
#define grid_size ((int)(f_to_go / dx))
#define t_to_go 0.2
#define gerr_max 1e-10

#define write_posix 0
#define write_netcdf 1

#if write_netcdf
#define ncFileName "grid.nc"
#define CHECK(e)                                                    \
    {                                                               \
        int err = e;                                                \
        if (err) {                                                  \
            fprintf(stderr, "NetCDF error: %s\n", nc_strerror(err));\
            MPI_Abort(MPI_COMM_WORLD, err);                         \
        }                                                           \
    }

typedef int NETCDFID;
#endif
