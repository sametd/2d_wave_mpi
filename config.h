#define dx 0.02
#define dt 0.001
// Gauss Seidel off diag
#define off_diag ((0.25) * (dt * dt) / (dx * dx))
// Gauss Seidel on diag
#define on_diag (1 + (4 * off_diag))
#define f_to_go (3.14159) / 2
#define grid_size ((int)(f_to_go / dx))
#define t_to_go 1.5
#define gerr_max 1e-10

#define write_posix 0
#define write_netcdf 1

#ifdef write_netcdf
#define ncFileName "grid.nc"
#define CHECK(e)                                   \
    {                                              \
        int err = e;                               \
        if (err) {                                 \
            printf("Error: %s\n", nc_strerror(e)); \
            exit(EXIT_FAILURE);                    \
        }                                          \
    }

typedef int NETCDFID;
#endif
