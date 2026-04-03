#define dx 0.00594998
#define dt 0.001
// r² = dt²/dx², the CFL-like ratio for the implicit wave equation
#define off_diag ((dt * dt) / (dx * dx))
// diagonal: 1 + 4r²
#define on_diag (1.0 + 4.0 * off_diag)
#define f_to_go (3.14159)
#define grid_size ((int)(f_to_go / dx))
#define t_to_go 3.0
#define gerr_max 1e-6

#define write_posix 0
#define write_netcdf 1
#define write_tensogram 0

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

#if write_tensogram
#define tgmFileName "grid.tgm"
#define tgm_use_packing 1
#define tgm_bits_per_value 24
#define tgm_use_szip 1
#define tgm_szip_rsi 128
#define tgm_szip_block_size 16
#define tgm_szip_flags 8
#define TGM_CHECK(e)                                                            \
    {                                                                           \
        enum tgm_TgmError tgm_err = (e);                                       \
        if (tgm_err != OK) {                                                    \
            const char *tgm_msg = tgm_last_error();                             \
            fprintf(stderr, "Tensogram error: %s\n", tgm_msg ? tgm_msg : "?"); \
            MPI_Abort(MPI_COMM_WORLD, (int)tgm_err);                            \
        }                                                                       \
    }
#endif
