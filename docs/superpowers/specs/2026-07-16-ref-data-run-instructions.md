# AMPS Reference-Data Runs — Cluster Instructions (M0)

Branch: `cloudlab_port` of jcanton/scale_amps (instrumented
`scale_atmos_phy_mp_amps.F90`). Build exactly as usual (same SCALE_SYS,
same Makefile flags); the instrumentation is inert unless enabled by
namelist.

## 1. Pre-flight

Do these before launching either run:

- **Fix the AMPS_DATA paths.** `AMPSTASK.F` hard-codes `DRCETB`/`DRAPTB`/
  `DRSTTB` at `/cluster/scratch/congchia/scale_amps/...`
  (`AMPSTASK.F:203,209,212`). Edit these three paths to point at the
  AMPS_DATA location available on whatever machine/account is actually
  running the job — unresolved paths abort at AMPS setup.
- **Reconcile the 16-vs-32 rank mismatch.** The committed `job_run.sh`/
  `job_restart_run.sh` launch `mpirun -n 16` with `--ntasks=16,
  --cpus-per-task=3` (16 MPI ranks x 3 OpenMP threads each), but the
  committed `run.conf`/`restart_run.conf` specify `PRC_NUM_X=4,
  PRC_NUM_Y=8` (32 MPI domains). These do not match as committed — fix
  one of the two before submitting (the job script's `-n` count must
  equal `PRC_NUM_X * PRC_NUM_Y`). Whichever is actually used determines
  the true rank count, and total dump volume scales with it (rank count
  x per-rank volume, not a fixed "32 ranks" figure — see §2 for per-rank
  estimates).
- `mkdir -p amps_dump` next to the run before enabling dumps (§2) — the
  directory is not created for you; a missing one aborts cleanly at AMPS
  setup with a clear error.
- **Dump files are clobbered between runs.** `AMPS_DUMP_open` uses
  `status='replace'` with fixed file names
  (`amps_dump_r<rank>_t<thread>.bin`). Running the spin-up and then the
  seeding run with the same `amps_dump_dir` overwrites the first run's
  files. Before launching the next run, either set a distinct
  `amps_dump_dir` per run (in each run.conf), or tar off the previous
  run's directory (§4) and clear it —
  `rm -rf amps_dump && mkdir -p amps_dump` — before the next launch.

## 2. Enable dumps

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

Notes:
- Indices are LOCAL (per-rank, halo-inclusive; first interior point is 3
  with IHALO=JHALO=2). Every rank dumps its own 2x2 box -> 4 columns per
  rank.
- `amps_dump_step_stride=300` with TIME_DT_ATMOS_PHY_MP=1.0 s dumps every
  300 s. `TIME_AMPS` starts at 0 and the dump condition is
  `mod(TIME_AMPS, amps_dump_step_stride) == 0`, so the FIRST MP step
  (TIME_AMPS=0) is always dumped whenever dumping is enabled, regardless
  of the stride value. `amps_dump_step_stride` must be >= 1 — a value of
  0 makes the `mod()` a divide-by-zero and crashes the run; there is no
  "disable via stride=0", use `l_amps_dump = .false.` instead.
- **Dump counts (corrected):** the spin-up runs `TIME_AMPS=0..3600`
  inclusive of both ends (a restart-time force-tendency call at t=0, plus
  the final step at t=3600 executing before the loop exits), so it dumps
  at `t=0,300,600,...,3600` — **13** dump times, not 12. Each seeding run
  resets `TIME_AMPS` to 0 and runs `0..1200`, dumping at
  `t=0,300,600,900,1200` — **5** dump times.
- Volume estimate, warm run (npr=6,nbr=40, npi=18,nbi=20, npa=5,nba=1,nca=4):
  <= ~0.35 MB per record pair per cloudy column-call; worst case
  13 x 4 x 2 x 0.35 MB ~ 37 MB per rank. Usually much less (only cloudy
  levels are in the packed vectors; micro dumps only fire when the column
  has cloudy points). Seeding-run records are larger: `num_h_bins=40,40`
  (vs `40,20` in the spin-up) roughly doubles the ice-array size, so each
  seeding record pair is ~1.5-2x the warm-run estimate (~0.5-0.7 MB);
  worst case 5 x 4 x 2 x ~0.6 MB ~ 12 MB per rank, still well under the
  spin-up total.

## 3. Runs wanted

1. **Warm spin-up** (`run.conf`): primary M2 validation target. It is
   warm only because the IN aerosol category is empty (`N_ap_ini(2)=0`,
   `AMPSTASK.F:198`) and the drop-freezing modes (contact/immersion/
   homogeneous, `micexfg(14/16/17)`) are off — **not** because ice
   microphysics is disabled. Deposition nucleation and every other
   ice-process code path still execute on every cloudy micro call and
   simply no-op (no dust/IN aerosol present to nucleate on). Do not
   describe this run as "ice off in AMPSTASK sense"; the port must match
   a run where ice routines execute and no-op, not one where they are
   absent.
2. **Seeding run** (`restart_run.conf`, 40 ice bins): M3/M4 target.
   `RESTART_SKIP_READING_AMPS_ICE=.true.` makes this run start ice-free by
   design — it skips reading *all* ice tracers (`*_imass`, `*_axis`,
   `ex_cry`) from the spin-up restart file, regardless of what the
   spin-up produced (and regardless of the 20-vs-40-bin change).
   Caveat: the skipped `QTRC` slices are never explicitly zeroed by the
   driver — `QTRC` is allocated without initialization and there is no
   `QTRC=0` assignment before the restart read, so this is formally
   undefined even though it is benign in practice (fresh allocations are
   zero pages on Linux). Recommendation, **not applied**: a one-line
   `QTRC(:,:,:,iq) = 0.0_RP` patch inside the skip branch of
   `ATMOS_vars_restart_read` (`mod_atmos_vars.F90`, near line 1152) would
   make this safe by construction. Leave it as a documented option for
   whoever hardens this path later — do not apply it for these reference
   runs, since it would silently change what the reference dumps
   represent.

Other verified config facts that apply to both runs:

- `MP_do_precipitation=.false.` (`run.conf:137`, `restart_run.conf:138`)
  is **required** by the AMPS driver (it aborts if `.true.`) and does
  **not** disable AMPS-internal sedimentation. AMPS does its own
  sedimentation under `l_sediment=.true.`, and the sed dumps (record
  phases 3/4) fire normally in both runs whenever the dumped column has
  liquid or ice. Flipping `l_sediment=.false.` would silently remove all
  phase-3/4 records — don't.
- Aerosols are fully prognostic in both runs (`l_fix_aerosols=.false.`):
  the per-step refill to the initial profile is skipped, so aerosol state
  (`qapvm`/`qapv`) evolves via activation/nucleation scavenging,
  advection, and seeding only. Port validation must take aerosol state
  from the dump records, never regenerate it from the `AMPSTASK.F`
  initial profiles — it will diverge after the first MP step.
- `l_no_ice_heat=.true.` would silently produce **empty** dumps (all six
  `AMPS_DUMP_active` call sites live only in the main code path used when
  `l_no_ice_heat=.false.`). Both configs currently set `.false.`
  (`run.conf:152`, `restart_run.conf:153`) — fine as committed, but treat
  this as a hard requirement, not an incidental default, if either config
  is ever touched.
- Cosmetic: the spin-up's cloud-seeding window (`DO_CLOUD_SEEDING=.true.`
  in `run.conf`) opens exactly at spin-up end, so the last iteration
  computes a seeding tendency that is then discarded (the restart is
  written first). Harmless as committed. For cleaner provenance you may
  optionally set `DO_CLOUD_SEEDING=.false.` in the spin-up `run.conf` —
  no behavioral change for the committed durations.

**`ihabit_gm_random` — verified to have zero effect on cloudlab** (see
`docs/superpowers/facts/verification/rngv_verdict.md`). The flag is `1`
in `AMPSTASK.F:85`. Every consumer of the RNG-drawn habit is guarded at
`T <= 253.16 K` (−20 °C), and the multi-habit `frq` lookup tables that
feed the draw are themselves only ever exercised at T ≤ −20 °C (that's
their warm edge; warmer temperatures clamp to that same row). No other
code anywhere in the repo consumes the RNG stream (`rand2_ty`/`rdsd`).
The cloudlab domain's coldest level starts at 266.09 K (−7.06 °C, ~13 K
of margin above −20 °C); reaching −20 °C within the ~1.3-2 h combined
run is physically implausible under the configured weak forcings. Net
effect: `ihabit_gm_random=1` vs `=0` is expected to be **bit-for-bit
identical** on cloudlab, so a "deterministic-habit rerun" adds no
validation value and has been dropped from the run list above. The
port's counter-based RNG (which cannot reproduce the Fortran LCG stream
bit-for-bit) is therefore irrelevant to cloudlab validation regardless of
which flag setting is used.

Residual check (the one thing not settled by code inspection alone):
after each run, confirm `min(T)` from that run's history output stays
**> −20.0 °C**. If a future config or forcing change ever pushes any cell
to ≤ −20 °C, this verdict no longer holds and the RNG-stream divergence
becomes real.

## 4. Collect results

Also copy the exact `run.conf` (or `restart_run.conf` for the seeding
run) and `AMPSTASK.F` that were actually used for the run into the
tarball — per-call replay/validation against the reference data depends
on knowing the exact namelist and habit settings that produced it (the
spin-up and seeding run share the same `AMPSTASK.F` but differ in
`run.conf`/`restart_run.conf` — ice-bin count, seeding window,
`RESTART_SKIP_READING_AMPS_ICE`; without a copy of the config alongside
the dump, the two runs' tarballs are indistinguishable):

    cp run.conf AMPSTASK.F amps_dump/          # or restart_run.conf
    tar czf amps_dump_<runname>.tar.gz amps_dump/

Copy the tarballs back. Locally they are parsed with:

    python3 scripts/amps_dump_reader.py amps_dump/ -o amps_ref_<runname>.npz

## 5. Sanity checks after the run

- `amps_dump/` contains `amps_dump_r<rank>_t<thread>.bin` files with
  nonzero size for ranks whose local 2x2 box had cloud.
- The reader completes without "bad magic" errors — if it fails, the
  record layout in the Fortran patch and `scripts/amps_dump_reader.py`
  have diverged; do not proceed, report back.
- `min(T)` from the run's history output is `> −20.0 °C` (the residual
  check from §3 on `ihabit_gm_random` irrelevance).
