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
  with IHALO=JHALO=2). Every rank dumps its own 2x2 box -> 4 columns per
  rank; total rank count depends on which launcher config is actually used
  (see the ranks-count caveat below) — 34 MB per rank either way.
- `amps_dump_step_stride=300` with TIME_DT_ATMOS_PHY_MP=1.0 s dumps every
  300 s: 12 dump times over the 1 h spin-up. `TIME_AMPS` starts at 0 and the
  dump condition is `mod(TIME_AMPS, amps_dump_step_stride) == 0`, so the
  FIRST MP step (TIME_AMPS=0) is always dumped whenever dumping is enabled,
  regardless of the stride value. `amps_dump_step_stride` must be >= 1 — a
  value of 0 makes the `mod()` a divide-by-zero and crashes the run; there
  is no "disable via stride=0", use `l_amps_dump = .false.` instead.
- Volume estimate, warm run (npr=6,nbr=40, npi=18,nbi=20, npa=5,nba=1,nca=4):
  <= ~0.35 MB per record pair per cloudy column-call; worst case
  12 x 4 x 2 x 0.35 MB ~ 34 MB per rank. Usually much less (only cloudy
  levels are in the packed vectors; micro dumps only fire when the column
  has cloudy points).
- **Ranks-count inconsistency**: the committed `job_run.sh`/`job_restart_run.sh`
  launch `mpirun -n 16` with `--ntasks=16`, `--cpus-per-task=3` (16 MPI
  ranks x 3 OpenMP threads each), but the committed `run.conf`/
  `restart_run.conf` specify `PRC_NUM_X=4, PRC_NUM_Y=8` (32 MPI domains).
  These do not match as committed; whichever one is actually used at
  submission time (job script's `-n` count must equal `PRC_NUM_X *
  PRC_NUM_Y`, so one of the two needs fixing before the run, not the dump
  code) determines the true rank count, and total dump volume scales with
  it (rank count x 34 MB, not a fixed "32 ranks" figure).

## 2. Runs wanted (in order of usefulness)

1. **Warm spin-up** (`run.conf`, ice off in AMPSTASK sense — spin-up config
   as-is): primary M2 validation target.
2. **Seeding run** (`restart_run.conf`, ice on, 40 ice bins): M3/M4 target.
3. **Seeding run with deterministic habit**: copy `AMPSTASK.F`, set
   `ihabit_gm_random=0`, rerun 2. Gives the tight ice comparison target
   (counter-based RNG in the port cannot match the Fortran LCG stream).

## 3. Collect results

Also copy the exact `run.conf` (or `restart_run.conf` for the seeding runs)
and `AMPSTASK.F` that were actually used for the run into the tarball —
per-call replay/validation against the reference data depends on knowing the
exact namelist and habit settings that produced it (e.g. run 3 in section 2
above only differs from run 2 by an edited `AMPSTASK.F`; without a copy of
that file alongside the dump, the two runs' tarballs are indistinguishable):

    cp run.conf AMPSTASK.F amps_dump/          # or restart_run.conf
    tar czf amps_dump_<runname>.tar.gz amps_dump/

Copy the tarballs back. Locally they are parsed with:

    python3 scripts/amps_dump_reader.py amps_dump/ -o amps_ref_<runname>.npz

## 4. Sanity checks after the run

- `amps_dump/` contains `amps_dump_r<rank>_t<thread>.bin` files with
  nonzero size for ranks whose local 2x2 box had cloud.
- The reader completes without "bad magic" errors — if it fails, the
  record layout in the Fortran patch and `scripts/amps_dump_reader.py`
  have diverged; do not proceed, report back.
