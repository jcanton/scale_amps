# Audit: ref-data run instructions vs. cloudlab configs (branch `cloudlab_port`)

Spec audited: `/Users/jcanton/projects/scale_amps/docs/superpowers/specs/2026-07-16-ref-data-run-instructions.md`. The 16-vs-32 rank note and `ihabit_gm_random` question are excluded as instructed.

---

## Quirk 1 — `MP_do_precipitation=.false.` vs AMPS sedimentation (sed dumps DO fire)

- **Instructions say**: nothing about `MP_do_precipitation`; they implicitly assume phase-3/4 sed records will be produced.
- **Config/code does**: `run.conf:137` and `restart_run.conf:138` set `MP_do_precipitation=.false.`. This flag controls the SCALE-level sedimentation path only (`scale-rm/src/atmos/mod_atmos_phy_mp_driver.F90:71,451`). For AMPS it is *mandatory*: the driver **aborts** if it is `.true.` (`mod_atmos_phy_mp_driver.F90:479-483`, "Precipitation should be off if AMPS is used"). AMPS does its own sedimentation internally, gated by the separate AMPS-namelist flag `l_sediment=.true.` (`run.conf:151`, `restart_run.conf:152`; `Sedimention_IF` at `scalelib/.../scale_atmos_phy_mp_amps.F90:2359`). The phase-3/4 dump calls sit inside that `Sedimention_IF` block, before/after `sclsedprz_original`, guarded only by `AMPS_DUMP_active` and `k2r>0` (liquid, lines 2402/2438) / `k2i>0` (ice, lines 2457/2493). Since `calc_precipitation=.false.`, the driver zeroes `SFLX_rain/snow` once at setup (`mod_atmos_phy_mp_driver.F90:508-516`) and AMPS itself fills them.
- **Consequence**: sed-record dumps fire normally in both runs whenever the dumped column has liquid (phase isn=0) or ice (isn=1) at a dump step. No conflict.
- **Recommended fix**: add one clarifying sentence to the spec: "`MP_do_precipitation=.false.` is required by the AMPS driver and does NOT disable AMPS-internal sedimentation (`l_sediment=.true.`); sed records (phases 3/4) are produced normally." Also warn that flipping `l_sediment=.false.` would silently remove all phase-3/4 records.

---

## Quirk 2 — `l_fix_aerosols=.false.` (module default is `.true.`)

- **Instructions say**: nothing about aerosol treatment.
- **Config/code does**: both configs set `l_fix_aerosols=.false.` (`run.conf:150`, `restart_run.conf:151`), overriding the module default `.true.` (`scale_atmos_phy_mp_amps.F90:135`). With `.false.`, the per-step refill of `qapv` back to the initial profile is skipped (`scale_atmos_phy_mp_amps.F90:1982-1994`), so aerosols are fully prognostic: depleted by activation/nucleation scavenging, advected, and replenished only by the seeding tendency. (Note also `fix_aerosol_type=.false.,...` — even with `l_fix_aerosols=.true.` nothing would refill, since the refill is additionally gated per-category at line 1985.) `l_fill_aerosols=.false.` too (line 153/154).
- **Consequence**: the aerosol arrays (`qapvm`, `qapv`) captured in the phase-1/3 dumps evolve over the run and are NOT reconstructible from `AMPSTASK.F` initial profiles. Per-call replay is unaffected (inputs are in the record), but any validation logic in the port that regenerates aerosols from `N_ap_ini`/`ini_aerosol_prf` (the default-`.true.` behaviour) will diverge from the reference after the first MP step.
- **Recommended fix**: state in the spec that aerosols are prognostic in these runs (`l_fix_aerosols=.false.`) and that the port must take aerosol state from the dump record, never from the initial profile.

---

## Quirk 3 — spin-up "warm phase" is emergent, not switched off; exact ice-formation conditions

- **Instructions say**: run 1 is "Warm spin-up (`run.conf`, ice off in AMPSTASK sense — spin-up config as-is)".
- **Config/code does**: `AMPSTASK.F` (shared by BOTH runs — read from the run CWD, `contrib/AMPS/mod_amps_utility.F90:1398`) does **not** turn ice off. Active `micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,1,0,0,0,0,1,0,1` (`AMPSTASK.F:286`): ice-ice/liq-ice collision (3,4), vapor deposition on ice (7), melting (8), master ice nucleation (10)=1 and **depositional nucleation (13)=1** are all ON; contact (14), splinter (15), immersion (16), homogeneous (17) are OFF. Ice bins exist in the spin-up (`num_h_bins=40,20`, `run.conf:145`). `Ice_Nucleation2` runs every cloudy micro call (`contrib/AMPS/class_Cloud_Micro.F90:1323-1330`) and calls `deposition_mode_vec` because `micexfg(13)/=0` (`contrib/AMPS/mod_amps_core.F90:3240-3243`).
  The code conditions for ANY ice to nucleate (`mod_amps_core.F90:17153-17155`) are:
  1. `T < 273.16 K`,
  2. `Si >= 0` (supersaturation over ice at t+dt, `s_v_n(2)`),
  3. **dust/IN aerosol category `ga(2)` has `con >= 1e-22` AND `mass(1) >= 1e-22`**.
  Condition 3 is the only real gate: `N_ap_ini = 317.0, 0.0, 0.0, 0.0` (`AMPSTASK.F:198`) leaves category 2 empty, `l_fix_aerosols=.false.` never refills it, and immersion/contact/homogeneous freezing (the drop-freezing paths, `class_Cloud_Micro.F90:1131-1135`, `mod_amps_core.F90:13513` requires supercooling ≥ 5 K anyway) are all disabled. The nucleation rate additionally uses `nucleation_halflife=0.0001 > 0` → the `ICE_FF` sigmoid (`mod_amps_core.F90:17320-17328`), which is ~0 above ≈ −9 °C.
  Edge case: `run.conf` has `DO_CLOUD_SEEDING=.true.` with a window opening exactly at spin-up end — LOWER = 18:00:00, spin-up runs 17:00:00→18:00:00 (`run.conf:38,40,209-218`; window test `mod_user.F90:490-492` is on NOWDATE, which is the end-of-step time). `USER_calc_tendency` at the final iteration therefore *does* compute a seeding tendency, but the main loop writes the restart (`mod_rm_driver.F90:561`) **before** `USER_calc_tendency` (line 570) and exits (line 583) without another update, so the tendency is discarded and the restart is uncontaminated.
- **Consequence**: spin-up produces no ice **only because the IN aerosol category is empty**, not because ice microphysics is off. The reference micro records still exercise `Ice_Nucleation2`/`deposition_mode_vec` code paths (with zero dust) and, with `ihabit_gm_random=1`, still consume RNG draws in `cal_growth_mode_hex_inl_vec` (`mod_amps_core.F90:17162`, draws gated at `mod_amps_utility.F90:1012`). The port must match a run where ice routines execute and no-op, not a run where they are absent.
- **Recommended fix**: reword run 1's description to "warm because the IN category (`N_ap_ini(2)=0`) is empty and immersion/contact/homogeneous freezing are off in `micexfg`; deposition nucleation and all ice process code are compiled in and executed". Warn that changing `N_ap_ini`, `micexfg(14/16/17)`, or extending `TIME_DURATION` past 18:00:00 (seeding window opens then) makes the spin-up icy.

---

## Quirk 4 — `RESTART_SKIP_READING_AMPS_ICE=.true.` and the 20→40 ice-bin change

- **Instructions say**: run 2 is "Seeding run (`restart_run.conf`, ice on, 40 ice bins)" — no mention of how the restart tolerates the tracer-count change.
- **Config/code does**: `restart_run.conf:62` sets `RESTART_SKIP_READING_AMPS_ICE=.true.` (declared/read at `scale-rm/src/admin/mod_admin_restart.F90:52,215`, passed at line 500). In `ATMOS_vars_restart_read` (`scale-rm/src/atmos/mod_atmos_vars.F90:1145-1160`) any tracer whose name contains `imass`, `axis`, or `ex_cry` is **not read** from the restart file. All 16 ice-PPV name templates (`scale_atmos_phy_mp_amps.F90:182-197`: `*_imass`, `*_axis`, `ex_cry`) match; no liquid (`*_lmass`) or aerosol (`*_amass`) name does. So the 40,20→40,40 bin change (`run.conf:145` vs `restart_run.conf:146`) is handled by simply never reading ANY ice tracer: the run never asks the file for the nonexistent `ice_imass21..40`, and the seeding run starts ice-free (spin-up ice, had there been any, would be dropped too).
- **Consequence / caveat**: the skipped `QTRC` slices are never explicitly zeroed — `QTRC` is allocated without initialization (`mod_atmos_vars.F90:546`) and no `QTRC=0` assignment exists before the restart read. In practice fresh large allocations are zero pages on Linux, but this is formally undefined; a compiler/allocator change (or `-finit-real=snan` debugging builds) would poison every ice tracer and crash `ATMOS_vars_check` (`ATMOS_VARS_CHECKRANGE=.true.`, `restart_run.conf:164`).
- **Recommended fix**: add a note that the seeding run intentionally starts with zero ice regardless of spin-up contents; optionally patch `mod_atmos_vars.F90` to set `QTRC(:,:,:,iq)=0.0_RP` when a tracer is skipped (one line inside the `if (skip...)` branch at line 1152) before these runs are treated as canonical.

---

## Quirk 5 — dump-time counts: 13 (not 12) spin-up dumps, 5 seeding dumps

- **Instructions say**: "12 dump times over the 1 h spin-up"; volume 12 × 4 × 2 × 0.35 MB ≈ 34 MB/rank; nothing about the seeding run's dump count or record size.
- **Config/code does**: on resume, `restart_read` calls `ATMOS_driver_calc_tendency(force=.true.)` (`mod_rm_driver.F90` restart_read, line ~98 of the subroutine) → MP runs at t=0 with `TIME_AMPS=0`. The main loop then runs MP at every whole second including the **final** step: `calc_tendency` executes at the last iteration before `TIME_DOend` exits (`mod_rm_driver.F90:565,583`), and the MP flag fires at t=3600 (`mod_admin_time.F90:1017-1021`). `TIME_AMPS` increments once per MP call (`scale_atmos_phy_mp_amps.F90:3798`), so spin-up covers `TIME_AMPS=0..3600` → dumps at 0,300,…,3600 = **13 dump times** (~37 MB/rank worst case, not 34). Seeding run: `TIME_AMPS` is reset to 0 at setup of the new executable (`scale_atmos_phy_mp_amps.F90:834`), `TIME_AMPS=0..1200` → dumps at 0,300,600,900,1200 = **5 dump times**. Seeding-run records are also bigger: `num_h_bins=40,40` doubles `nbi` (ice arrays dominate the record at npi=18), so a record pair is roughly 1.5–2× the warm 0.35 MB estimate — still well under the warm run's total (≈5×4×~0.6 MB ≈ 12 MB/rank worst case).
- **Consequence**: tarball-size expectations and the reader's expected time axis are slightly off; a validator asserting exactly 12 dump times will fail.
- **Recommended fix**: update the spec to "13 dump times (t = 0, 300, …, 3600 s) for the spin-up; 5 (t = 0, 300, 600, 900, 1200 s) for each seeding run", and note the seeding run's larger per-record ice block (nbi=40).

---

## Quirk 6 — other mismatches

**(a) Hard-coded user paths in `AMPSTASK.F`.** `DRCETB`/`DRAPTB`/`DRSTTB` point at `/cluster/scratch/congchia/scale_amps/...` (`AMPSTASK.F:203,209,212`); `AMPSTASK.F` itself is opened by relative path from the run CWD (`mod_amps_utility.F90:1398`). Instructions say "build exactly as usual" and never mention these. Anyone other than that user gets an abort at AMPS setup. **Fix**: add a pre-flight step "edit the three AMPS_DATA paths in AMPSTASK.F to your scratch location".

**(b) Dump files are silently clobbered between runs.** `AMPS_DUMP_open` uses `status='replace'` with fixed names `amps_dump_r<rank>_t<thread>.bin` (`scale_atmos_phy_mp_amps.F90:5205-5212`). Spin-up, seeding run, and the deterministic-habit rerun all executed in the same directory with the same `amps_dump_dir` overwrite each other. The spec's section 3 tars after each run but never says to move/empty the directory first. **Fix**: instruct "tar and then delete/rename `amps_dump/` (or set a distinct `amps_dump_dir` per run) before launching the next run".

**(c) Dumps exist only in the `l_no_ice_heat=.false.` code path.** All six `AMPS_DUMP_active` call sites (lines 2091-2493) are in the main OpenMP region; the alternate water-only region used when `l_no_ice_heat=.true.` (second parallel region starting ~line 3076, see comment at 3047) has none. Both configs set `.false.` (`run.conf:152`, `restart_run.conf:153`) so this is fine as committed — but flipping that flag silently produces empty dump files. **Fix**: one warning line in the spec.

**(d) Spin-up seeding window overlap (cosmetic).** As detailed in Quirk 3: `run.conf`'s seeding window [18:00:00, 19:01:30) touches the spin-up's final NOWDATE; the tendency is computed but discarded, and history at the last output time can show nonzero `RHOQ_t_SEED_*`. Harmless as committed; becomes real seeding if `TIME_DURATION` is ever extended. **Fix**: note it, or set `DO_CLOUD_SEEDING=.false.` in the spin-up `run.conf` for the reference runs (cleaner provenance, zero behavioural change for the committed durations).

**(e) Namelist placement verified correct.** The dump keys and group `&PARAM_ATMOS_PHY_MP_AMPS_bin` in the instructions match the instrumented code exactly (`scale_atmos_phy_mp_amps.F90:357-381`); stride-0 divide-by-zero and "first step always dumped" claims are confirmed (`AMPS_DUMP_active`, lines 5229-5235); missing `amps_dump` directory aborts cleanly with a clear error (lines 5213-5216).