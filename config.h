#ifndef CONFIG_H
#define CONFIG_H

/* ── Physics ── */
#define dx 0.00594998
#define dt 0.001
#define off_diag ((dt * dt) / (dx * dx))       /* r² = dt²/dx² */
#define on_diag (1.0 + 4.0 * off_diag)         /* 1 + 4r² */
#define f_to_go (3.14159)
#define grid_size ((int)(f_to_go / dx))
#define t_to_go 3.0
#define gerr_max 1e-6

/* ── Output toggles ── */
#define write_posix 0
#define write_netcdf 0
#define write_tensogram 1

/* ── Output frequency: write every Nth solver step ── */
/* dt=0.001 means 3001 steps for t=3.0                */
/*   write_interval=1  → 3001 frames, 6.69 GB raw     */
/*   write_interval=10 →  301 frames,  670 MB raw     */
/*   write_interval=20 →  151 frames,  337 MB raw     */
#define write_interval 10

/* ── Tensogram benchmark: max frames to capture (0 = all) ── */
#define tgm_bench_steps 0

#endif
