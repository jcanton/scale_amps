# AMPS Reference-Data Runs — Cluster Instructions (M0)

Branch: `cloudlab_port` of jcanton/scale_amps (instrumented
`scale_atmos_phy_mp_amps.F90`). Build exactly as usual (same SCALE_SYS,
same Makefile flags); the instrumentation is inert unless enabled by
namelist.

## 1. Enable dumps

Add to the `&PARAM_ATMOS_PHY_MP_AMPS_bin` group of the run config
(`scale-rm/test/case/cloudlab/scripts/run.conf` for the warm spin-up;
`restart_run.conf` for the ice seeding run):

    l_amps_dump          = .true.,
    amps_dump_dir        = "./amps_dump",
    amps_dump_step_stride = 300,
    amps_dump_is         = 3,
    amps_dump_ie         = 4,
    amps_dump_js         = 3,
    amps_dump_je         = 4,

Create the directory next to the run: `mkdir -p amps_dump`.

Notes:
- Indices are LOCAL (per-rank, halo-inclusive; first interior point is 3
  with IHALO=JHALO=2). Every rank dumps its own 2x2 box -> 32 ranks x 4
  columns.
- `amps_dump_step_stride=300` with TIME_DT_ATMOS_PHY_MP=1.0 s dumps every
  300 s: 12 dump times over the 1 h spin-up.
- Volume estimate, warm run (npr=6,nbr=40, npi=18,nbi=20, npa=5,nba=1,nca=4):
  <= ~0.35 MB per record pair per cloudy column-call; worst case
  12 x 4 x 2 x 0.35 MB ~ 34 MB per rank, ~1 GB total. Usually much less
  (only cloudy levels are in the packed vectors; micro dumps only fire
  when the column has cloudy points).

## 2. Runs wanted (in order of usefulness)

1. **Warm spin-up** (`run.conf`, ice off in AMPSTASK sense — spin-up config
   as-is): primary M2 validation target.
2. **Seeding run** (`restart_run.conf`, ice on, 40 ice bins): M3/M4 target.
3. **Seeding run with deterministic habit**: copy `AMPSTASK.F`, set
   `ihabit_gm_random=0`, rerun 2. Gives the tight ice comparison target
   (counter-based RNG in the port cannot match the Fortran LCG stream).

## 3. Collect results

    tar czf amps_dump_<runname>.tar.gz amps_dump/

Copy the tarballs back. Locally they are parsed with:

    python3 scripts/amps_dump_reader.py amps_dump/ -o amps_ref_<runname>.npz

## 4. Sanity checks after the run

- `amps_dump/` contains `amps_dump_r<rank>_t<thread>.bin` files with
  nonzero size for ranks whose local 2x2 box had cloud.
- The reader completes without "bad magic" errors — if it fails, the
  record layout in the Fortran patch and `scripts/amps_dump_reader.py`
  have diverged; do not proceed, report back.
