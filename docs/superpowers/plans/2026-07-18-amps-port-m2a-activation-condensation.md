# AMPS→icon4py Port — M2a (Warm-Phase: host loop + liquid diagnosis + activation + condensation) Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** First validated warm-phase physics: the operator-split host loop (warm subset), liquid `diag_pq` (per-bin derived quantities incl. terminal velocity), CCN activation + vapor advancement (`cal_aptact_var8_kc04dep`, the Brent supersaturation solve), liquid vapor deposition (Chen–Lamb growth + mass-space bin remap), and liquid budget repair — validated per-call against the cloudlab spin-up reference dumps.

**Architecture:** Numpy reference implementations are the M2a deliverable (they become `StencilTest`/box-replay references and are what we validate against Fortran dumps); GT4Py DSL kernels follow where pointwise/spike-proven, but M2a's correctness gate is the numpy path matching dumped Fortran state to tight tolerance. Host loop is Python orchestration (substep loops), each process a numpy function now, DSL-migratable later. Per the M0 gate + M1 foundations: CGS internally, Fortran-exact table semantics, counter-RNG (idle in warm phase — activation/condensation don't roll habit dice above −20°C).

**Ground truth (committed on `cloudlab_port`), referenced as G1–G6:**
- G1 = `docs/superpowers/facts/m2/micro-tendency-orchestration.md` (cal_micro_tendency substep skeleton, ifc_cloud_micro sequence, diag_t, update_mesrc/update_group_all, refresh preamble)
- G2 = `docs/superpowers/facts/m2/activation.md` (cal_aptact_var8_kc04dep, zbrent_act_vec, func_liqvap/icevap_vec, cal_coef_svsteady_init, Köhler radius)
- G3 = `docs/superpowers/facts/m2/vapor-deposition.md` (vapor_deposition liquid path, Chen–Lamb, cal_transbin/shift_bin remap, diag_pq liquid, capacitance/ventilation coeffs)
- G4 = `docs/superpowers/facts/m2/coalescence.md` (M2b — not this plan)
- G5 = `docs/superpowers/facts/m2/sedimentation-terminalvel.md` (w_terminal_vel — needed by diag_pq; repair driver)
- G6 = `docs/superpowers/facts/m2/breakup-tables-and-icon4py-m2.md` (M2b; PART B icon4py composition template applies here)
- M1 fact files F1–F5 (`docs/superpowers/facts/m1/`) still apply for constants/thermo/index/state.

**Reference data:** cloudlab spin-up dumps (per `docs/superpowers/specs/2026-07-16-ref-data-run-instructions.md`), parsed by `driver/ref_data.py`. The warm spin-up exercises activation + condensation + liquid diagnosis (no ice: IN category empty) — the exact M2a validation target. Per-call replay via `driver/box.py::case_from_micro_record`.

**Tech Stack:** Python 3.12, gt4py==1.1.11, numpy, pytest; icon4py worktree `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics` (branch `amps_microphysics`, `$WT`).

## Global Constraints

- Transcribe physics VERBATIM from G1–G3/G5 + F1–F5; if a fact section is truncated/ambiguous, the dispatch authorizes reading the named Fortran routine directly (quote what you read into the report) — do NOT guess a formula.
- Units CGS internally (AMPS-native), float64. Reuse M1 modules: `core/constants.AmpsConst`, `core/thermo` (esat, coefficients, diag_t/theta_il), `core/bin_grid`, `config.AmpsConfig`, `core/index_maps`, `core/state`, `core/lookup_tables`, `core/rng`, `core/packing`, `driver/ref_data`/`box`.
- Iterative solvers: fixed-max-iteration masked loops (Brent → the exact ITMAX + bracket-expansion from G2; do NOT substitute scipy — reproduce the Fortran iteration so results match). Repair: host-side loop around the rescale with the Fortran's NITER/NITER_TOTAL bounds and global any() convergence (G5).
- Bin remap (cal_transbin/shift_bin): gather-formulated over fixed bin count (G3) — numpy now; the Fortran is already dense i×ibx.
- License headers; ruff/mypy/tach/pre-commit clean; TDD; per-task report to `$WT/.superpowers/sdd/m2a-task-N-report.md`; branch `amps_microphysics`; commit per task, push after each REVIEWED task.
- Validation tolerance: warm-phase per-call replay rtol ~1e-8 initially (float64, same ops, RNG idle), tightening toward 1e-10 as each process is isolated; where a process has a documented cancellation/iteration-order sensitivity (like activation's Brent), loosen with a recorded justification (spike-B precedent).
- Test evidence: unit tests assert against values recomputed from the fact formulas in-test AND, where a dump exists, against the dumped Fortran arrays. If dumps are not yet available locally, structure the replay tests to be skipped-with-marker (datatest-style) pending the cluster run, but the numpy-vs-formula unit tests must be unconditional.

---

### Task 1: `implementations/warm_loop.py` — operator-split host-loop skeleton

**Files:**
- Create: `$WT/.../amps/implementations/__init__.py`, `$WT/.../amps/implementations/warm_loop.py`
- Test: `tests/amps/unit_tests/test_warm_loop.py`

**Interfaces:**
- Produces: `WarmLoopState` (bundles the M1 LiquidState/AerosolState/ThermoState + air group scalars for one packed column-vector) and `run_warm_micro_tendency(state, config, dt, luts) -> WarmLoopState` — the warm subset of `cal_micro_tendency` (G1 §1): the col_loop×n_step_cl / vap_loop×n_step_vp structure, the refresh preamble (update_mesrc→diag_t→update_airgroup→diag_pq) as ONE function `_refresh_state` (collapsing G1's 4× duplication), and STUBBED process hooks (`_activation`, `_vapor_deposition_liquid`, `_repair` raise NotImplementedError with the task that fills them). This task delivers the skeleton + refresh + substep bookkeeping + dt_cl/dt_vp setup ONLY; processes are Tasks 3–6.
- Also `ifc_warm(...)` mirroring G1 §2's ifc_cloud_micro sequence for the warm path (ini→reality_check→refresh→run_warm_micro_tendency→return_output), reusing packing/state from M1.

**Steps:** TDD. Test: substep counts (n_step_cl=1, n_step_vp=10 for cloudlab) drive the right number of `_refresh_state` calls; dt_cl=dt/n_step_cl, dt_vp=dt_cl/n_step_vp exact; the process stubs raise NotImplementedError (loud); `_refresh_state` on a known state reproduces diag_t's T from theta_il (reuse M1 thermo, assert vs F1 §5). Commit `feat(amps): warm-phase operator-split host-loop skeleton`.

---

### Task 2: `core/liquid_diag.py` — liquid `diag_pq` + terminal velocity

**Files:** Create `$WT/.../amps/core/liquid_diag.py`; Test `tests/amps/unit_tests/test_liquid_diag.py`

**Interfaces:** `diag_pq_liquid(liquid_state, thermo, config, luts) -> LiquidDiag` computing per-bin: mean mass, length (diameter, geometric from bin mass), density, terminal velocity (`w_terminal_vel`/`cal_terminal_vel_vec` from G5/G3 — Best-number/Reynolds fall speed), capacitance + ventilation coefficients (G3, needed by vapor deposition). VERBATIM from G3 (diag_pq liquid branch) + G5 (terminal velocity). This is a dependency of Tasks 3–6, so it comes first after the skeleton.

**Steps:** TDD. Tests: mean-mass = bin-mass-weighted per G3; terminal velocity spot-values recomputed from the G5 fall-speed formula for 3 drop sizes; density = liquid (den_w) for pure drops; ventilation/capacitance vs G3 formulas. Commit `feat(amps): liquid diag_pq (per-bin mean mass, length, density, terminal velocity, ventilation)`.

---

### Task 3: `core/activation.py` part 1 — the supersaturation Brent solver + vapor/temperature objective functions

**Files:** Create `$WT/.../amps/core/activation.py`; Test `tests/amps/unit_tests/test_activation_solver.py`

**Interfaces:** the contained solver primitives from G2: `func_liqvap` (liquid vapor objective), `func_icevap` (ice vapor objective incl. the KC04 deposition branch — kept for completeness though warm-phase Si path dominates), `cal_coef_svsteady_init` (steady-state coefficients), `cal_air_temp`, `zbrent_act` (Brent root-finder, EXACT ITMAX + bracket-expansion from G2), `func` (the top-level objective). Köhler critical/haze radius via M1... (M1 didn't port Köhler — this task adds `kohler_radius` from G2's get_critrad_anal/itr; prefer the analytic variant per M0 spike, reproduce the iterative only if analytic absent). All numpy, masked fixed-iteration.

**Steps:** TDD. Tests: zbrent converges to a known root of a synthetic monotone function within ITMAX; func_liqvap/icevap recomputed vs G2 for a known (T, S, coefficients) triple; Köhler critical radius vs the analytic formula for a known aerosol; bracket-expansion path exercised (a root outside the initial bracket). Commit `feat(amps): activation supersaturation Brent solver + vapor/temp objectives`.

---

### Task 4: `core/activation.py` part 2 — `cal_aptact_var8_kc04dep` driver (activation + vapor advancement)

**Files:** Extend `$WT/.../amps/core/activation.py`; Test `tests/amps/unit_tests/test_activation.py`

**Interfaces:** `activate_and_advance_vapor(liquid, aerosol, thermo, config, dt_vp, luts, diag) -> (liquid', aerosol', thermo')` — the full cal_aptact_var8_kc04dep production path (G2): activate dry aerosol → droplets, advance vapor over the substep, place activated droplets into liquid bins (add_simple/add_samebin), the iflg_dhf (micexfg 19) branches. Uses Task 3 solver + Task 2 diag. Wire into warm_loop `_activation` hook.

**Steps:** TDD. Tests: mass conservation (vapor + condensate) across an activation step; number of activated droplets vs a hand-computed Köhler-activation case; droplet placement into the correct bin; the DHF branch toggles with config.ice_nucleation_dhf. Then a per-call replay test (skipped-pending-dumps marker if no local dumps): feed a spin-up pre-record's aerosol+thermo through activation, compare the activated liquid + advanced vapor to the post-record, rtol per the tolerance ladder. Commit `feat(amps): CCN activation + vapor advancement (cal_aptact_var8_kc04dep)`.

---

### Task 5: `core/vapor_deposition.py` — liquid condensation/evaporation growth + bin remap

**Files:** Create `$WT/.../amps/core/vapor_deposition.py`; Test `tests/amps/unit_tests/test_vapor_deposition.py`

**Interfaces:** `vapor_deposition_liquid(liquid, thermo, config, dt_vp, diag) -> liquid'` — Chen–Lamb 1994 semidiscrete growth for liquid (G3): excess vapor density (cal_ex_vapor_density), per-bin mass growth via capacitance+ventilation coefficients, mass-space bin remap (cal_transbin/shift_bin, gather-formulated over fixed bins), mass-component distribution (dep_mass3). Wire into warm_loop `_vapor_deposition_liquid`.

**Steps:** TDD. Tests: a supersaturated single-bin state grows and shifts to the adjacent bin conserving mass+number per G3; subsaturated evaporates; the remap is mass-conserving (sum before==after to 1e-12); excess-vapor-density vs G3 formula. Replay test (marker-gated) vs spin-up dumps for the vapor substep. Commit `feat(amps): liquid vapor deposition (Chen-Lamb growth + mass-space bin remap)`.

---

### Task 6: `core/repair.py` — liquid budget repair + wire the full warm loop

**Files:** Create `$WT/.../amps/core/repair.py`; extend `implementations/warm_loop.py`; Test `tests/amps/unit_tests/test_repair.py`, extend `test_warm_loop.py`

**Interfaces:** `repair_liquid(liquid, config) -> liquid'` — the mass/concentration/volume non-negativity closure (G5 repair): the host-side NITER/NITER_TOTAL loop with global any() convergence, one liquid rescale block per G5. Replace the `_repair` stub. Then wire all three process hooks (`_activation`, `_vapor_deposition_liquid`, `_repair`) into `run_warm_micro_tendency` so a full warm substep runs end-to-end.

**Steps:** TDD. Tests: a state with a small negative bin mass is repaired to non-negative conserving total; the convergence loop terminates within NITER; end-to-end `run_warm_micro_tendency` on a spin-up pre-record runs without NotImplementedError and produces finite state. Commit `feat(amps): liquid budget repair; wire full warm micro-tendency loop`.

---

### Task 7: box-driver warm-phase wiring + per-call replay harness

**Files:** Extend `$WT/.../amps/driver/box.py`; Test `tests/amps/integration_tests/test_warm_replay.py` (new integration_tests dir + __init__)

**Interfaces:** implement `run_box(case)` for the warm path (case_from_micro_record → ifc_warm → BoxResult) — replacing the M1 NotImplementedError skeleton for warm cases (ice cases still raise). The replay harness: load spin-up reference (`driver/ref_data`), for each micro pre/post pair run the warm loop on the pre-state and compare to post, aggregate pass/fail with per-field rtol, report worst offenders.

**Steps:** TDD. Test (datatest-marker, skipped without dumps + a clear message pointing at the run instructions): end-to-end replay over available spin-up dump pairs, assert per-field rtol within the ladder; a synthetic-input smoke test runs unconditionally (no dump needed) proving run_box executes a full step. Commit `feat(amps): box-driver warm path + per-call replay harness`.

---

### Task 8: M2a wrap — full suite, mapping-doc update, README, push

**Steps:** full amps suite green + pristine; pre-commit clean; update `docs/superpowers/facts/m1/fortran-mapping.md` (in scale_amps) with the M2a modules→Fortran rows and the warm-phase not-yet list (coalescence/breakup=M2b, sedimentation=M2c); update package README (M2a status); commit + push both repos. Ledger + memory: M2a closed, note whether replay validation ran (dumps available?) or is pending the cluster run.

---

## Plan Self-Review (write-time)

- Spec coverage: M2a scope from the milestone def (activation+condensation+liquid diag_pq, validate vs dumps) — host loop (T1), liquid diag_pq+terminal vel (T2), activation solver (T3) + driver (T4), vapor deposition (T5), repair + wiring (T6), box replay (T7), wrap (T8). All covered.
- Dependency order: skeleton→diag(needed by all)→solver→activation→vapor dep→repair+wire→replay. Correct.
- Placeholders: physics referenced to G1–G5/F1–F5 fact sections (committed ground truth), not inlined — same pattern as M1, each task names its sections; dispatch authorizes direct Fortran read on truncation.
- Known risk: warm spin-up dumps may not be available at implementation time (cluster run is the user's). Mitigation: replay tests are marker-gated + skipped-with-pointer; unit tests (numpy-vs-formula) are unconditional and are the primary correctness gate until dumps arrive; M2a can complete and be reviewed on unit tests, with replay validation a follow-on once dumps land.
- Type consistency: LiquidDiag (T2) consumed by T3–T5; WarmLoopState (T1) threads through; reuses M1 state/thermo/packing signatures.
- M2b/M2c carry-forward: coalescence+breakup (port cal_breakfragment Cluster-1-only + validate vs dumped bu_fd/bu_tmass — add dump instrumentation), sedimentation Euler+PPM w/ momentum/energy. Not this plan.
