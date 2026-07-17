# AMPS→icon4py Port — M1 (Foundations) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver M1 of the AMPS port (spec §8): the `amps` package's foundations — constants, thermodynamics (Fortran-table semantics), bin grid, configuration, LUT conversion+loading, state bundles + codegen layer skeleton, SI↔CGS packing, in-package reference-data reader, box-driver skeleton, RNG/util graduation from spikes — plus the nine M0 carry-forward items.

**Architecture:** Pure-Python/numpy reference implementations are the M1 deliverable for every physics formula (they become the `StencilTest` references in M2+); GT4Py DSL variants are added in M1 only where the idiom is already spike-proven and pointwise (esat table accessor, thermo coefficient functions). All Fortran formulas/constants are transcribed from the committed fact extracts — NOT from memory.

**Fact extracts (ground truth, committed on `cloudlab_port`):**
- F1 = `docs/superpowers/facts/m1/constants-thermo.md`
- F2 = `docs/superpowers/facts/m1/bins-config-indexmaps.md`
- F3 = `docs/superpowers/facts/m1/lut-files.md`
- F4 = `docs/superpowers/facts/m1/state-packing-si-cgs.md`
- F5 = `docs/superpowers/facts/m1/icon4py-m1-conventions.md`

**Tech Stack:** Python 3.12, gt4py==1.1.11, numpy, pytest; icon4py worktree `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics` (branch `amps_microphysics`), `$WT` below. Fortran side: `/Users/jcanton/projects/scale_amps` branch `cloudlab_port`.

## Global Constraints

- Transcribe every constant/formula VERBATIM from the fact files (F1–F5 above); if a fact file section is ambiguous or missing something, STOP and report NEEDS_CONTEXT naming the section — do not guess and do not read the 28k-line Fortran yourself unless the dispatch says to.
- Units: CGS internally, exactly as Fortran (F1 documents each). Float64 everywhere.
- Package layout per spec §4: modules under `model/atmosphere/subgrid_scale_physics/amps/src/icon4py/model/atmosphere/subgrid_scale_physics/amps/` — `core/` for physics, package root for `state.py`/`config.py`, `codegen/`, `data/`, `driver/`. Tests under `tests/amps/unit_tests/` (plain pytest, no grid fixtures needed for numpy tests; F5 §5 confirms).
- Every new .py file: 7-line ICON4Py license header (as in existing package files). Ruff line length 100. `uv run` from worktree root. Run `uv run --group dev --frozen --isolated pre-commit run --all-files` before each commit; hooks may auto-fix — re-add and re-run until clean.
- Established idioms (M0 gate report — binding): LUTs as tiled or explicitly-broadcast-safe fields, never K-only gathered; esat/thermo tables keep EXACT Fortran index arithmetic — truncate-then-clamp `I=MAX(1,MIN(int(T)-163,...))` and the `wt` clamp, per F1 §3c (M0 carry-forward #4); codegen layer must be designed for SPLIT operators (per-destination-bin/banded), never one monolithic 40-bin operator; deep generated expressions need recursion limit 1e5 at compile; RNG = quadratic round + 3 Lehmer rounds mod M31 with constants as int64 params (F5 §6, spike_e).
- DSL limits: no for/while; no bitwise/shifts on ints; `offset_provider={}` for pointwise; never name a parameter `shift`; numpy>=2.0 — never embed `repr(np.float64)` in generated source (use `float()`).
- Commits: conventional style, one per task. icon4py commits on `amps_microphysics`, push to origin after each REVIEWED task. scale_amps commits on `cloudlab_port`, push after each task.
- Test evidence: every physics module's tests must assert against values COMPUTED from the fact-file formulas inside the test (self-derived expected values), plus at least one hand-checkable literal (e.g. a constant's exact value from F1).

---

### Task 1: `core/constants.py` — physical constants

**Files:**
- Create: `$WT/.../amps/core/constants.py`, `$WT/.../amps/core/__init__.py` (header-only)
- Test: `$WT/.../amps/tests/amps/unit_tests/test_constants.py`

**Interfaces:**
- Produces: `class AmpsConst(float, enum.Enum)` (muphys `GraupelConsts` pattern, F5 §1) with EVERY constant from F1 §1 (`mod_amps_const.F90`: p00, T_0, gg, M_w, R_u, MR, den_w, den_i, R_d, R_v, C_pa, Racp, Rdvchiarui, M_a, c_w, L_e, L_f, L_s, a_cliq, k_w, undef) and F1 §2 (`acc_amps.F90` precomputed: PI, sq_three, coef3s, coef4pi3, coef2p, coefpi6, coefsq2p, coef3sq3, coef4p, coef3i4p1i3, coedpi6, coedsq2p, coed3sq3, coefpi180). Derived constants (MR, Racp, Rdvchiarui, M_a) computed from the same expressions, not hardcoded duplicates. Docstring: CGS units, per-constant unit comments transcribed from F1.

- [ ] Step 1: failing test — `test_constants.py` asserts e.g. `AmpsConst.p00 == 1.0e6`, `AmpsConst.T_0 == 273.16`, `AmpsConst.R_d == 287.04e4`, `AmpsConst.MR == 18.016/8.31436e7`, `AmpsConst.Racp == 287.04e4/1004.64e4`, `abs(AmpsConst.coef4pi3 - 4*math.pi/3) < 1e-12` (and the literal 4.18879020478639 exactly), one assert per constant (literal from F1). Run: `uv run --group test --frozen pytest -n0 model/atmosphere/subgrid_scale_physics/amps/tests/amps/unit_tests/test_constants.py -v` → ImportError.
- [ ] Step 2: implement from F1 §§1–2 verbatim. Test passes.
- [ ] Step 3: pre-commit; commit `feat(amps): core constants (CGS, from mod_amps_const/acc_amps)`.

---

### Task 2: `core/thermo.py` — saturation tables, thermo coefficients, θ_il closure

**Files:**
- Create: `$WT/.../amps/core/thermo.py`
- Test: `tests/amps/unit_tests/test_thermo.py`

**Interfaces (all numpy; array-in/array-out, float64):**
- `make_esat_tables() -> tuple[np.ndarray, np.ndarray]` — estbar(150), esitbar(111) EXACTLY per F1 §3a QSPARM2 (Murphy–Koop; T starts 163., increments BEFORE evaluation → first entry at T=164; tables in Pa).
- `esat_lk(phase: int, t, estbar, esitbar) -> np.ndarray` — F1 §3c verbatim semantics: `I=MAX(1,MIN(int(T)-163,149))` (int() = truncation toward zero as Fortran), `wt=MAX(MIN(T-(I+163),1.0),0.0)`, linear interp, ×10.0 (Pa→dyn/cm²). Ice: 110 cap.
- `esat_analytic(phase: int, t)` — F1 §3b (Lowe & Ficke) verbatim incl. ×1000.
- `t_from_esat_lk(phase, t_guess, esat, estbar, esitbar)` — F1 §3d reverse search verbatim (±5-entry window then full-scan fallback — transcribe the full goto-free equivalent from F1's quoted code; preserve exact window bounds and fallback behavior).
- Thermo coefficients from F1 §3e+ (whatever F1 quotes for D_v, k_a, d_vis, GTP/MGTP, sig_wa — one function each, formula verbatim).
- `theta_il` closure: `cal_til`, `cal_thetail`, and `diag_t` equivalents per F1 §5 (T from θ_il + condensate loadings), numpy.
- Also: `_esat_lk_dsl` gt4py field_operator using the TILED table idiom (table passed as CellK field; spike_a/d precedent per F5 §6) + a thin `@gtx.program` wrapper — the ONE DSL deliverable of this task.

**Tests (all self-derived from F1 formulas):**
- Table construction: `estbar[0] == exp(54.842763 - 6763.22/164 - 4.210*log(164) + 0.000367*164 + tanh(0.0415*(164-218.8))*(53.878 - 1331.22/164 - 9.44523*log(164) + 0.014025*164))` (recompute in the test); same pattern for `esitbar[0]` and the last entries.
- Accessor: at T exactly 200.0, esat_lk == estbar[int(200)-163-1]*10 within fp (wt==0 case — mind the 1-based→0-based shift; document it in code); at T=200.5 the interpolated value; boundary bands T=163.2 and T=350.0 exercise the clamps — assert the exact Fortran clamped results (I=1 / I=149).
- Truncate-vs-clamp order regression: T=163.7 → Fortran gives `int(163.7)-163 = 0 → I=max(1,0)=1`, wt = max(min(163.7-164,1),0) = 0 → value == estbar[0]*10. Assert exactly this (this is carry-forward #4's named requirement).
- Round-trip: `t_from_esat_lk(1, 250.0, esat_lk(1, 250.3, ...), ...)` ≈ 250.3 within interpolation error.
- DSL accessor: embedded-backend comparison vs numpy `esat_lk` at rtol 1e-12 over a (32, 61) random T field in [180, 310]; gtfn_cpu same (small compile, fine in CI-local run).
- diag_t: construct θ_il from a known (T, q) state via `cal_thetail`, invert with `diag_t`, assert T recovered (tolerance per the Fortran iteration's own convergence, from F1 §5 — if F1 shows a fixed-point loop, replicate its iteration count/exit).

- [ ] Step 1: failing tests (write ALL above). Step 2: implement from F1. Step 3: run full unit_tests dir. Step 4: pre-commit; commit `feat(amps): thermo — Murphy-Koop tables, Fortran-exact accessors, theta_il closure`.

---

### Task 3: `core/bin_grid.py` — bin boundaries and masses

**Files:**
- Create: `$WT/.../amps/core/bin_grid.py`
- Test: `tests/amps/unit_tests/test_bin_grid.py`

**Interfaces:**
- `make_bin_grid(token: str, nbins: int, ...) -> BinGrid` — frozen dataclass with `binb` (boundaries, len nbins+1), `a_b`, `b_b`, mean masses, haze-split index; parameters and math transcribed VERBATIM from F2 §1 (`binmicrosetup_scale`): the liquid 40/80-bin and ice 10/20/40-bin constructions, the `binb(i+1)=a_b*binb(i)+b_b` recurrence with the exact seed values/exponents F2 quotes, `nbin_h` haze handling, and any aerosol-bin construction present.
- Validation in constructor: nbins ∈ allowed sets per F2 (`40|80` liquid, `10|20|40` ice) — ValueError otherwise.

**Tests:** boundaries strictly increasing; recurrence self-check `binb[i+1] == a_b*binb[i] + b_b` for all i at 1e-12; first/last boundary literals for the cloudlab config (40 liquid / 20 ice / 40 ice) computed inside the test from F2's quoted seed constants; mass-doubling-ish ratio checks if F2 shows them.

- [ ] TDD steps as Task 1/2; commit `feat(amps): bin grid construction (binmicrosetup_scale port)`.

---

### Task 4: `config.py` + `core/index_maps.py` — AmpsConfig and PPV index maps

**Files:**
- Create: `$WT/.../amps/config.py`, `$WT/.../amps/core/index_maps.py`
- Test: `tests/amps/unit_tests/test_config.py`

**Interfaces:**
- `AmpsConfig` frozen dataclass: every field of the `AMPSTASK.F` namelist (F2 §4/§5: level_comp, act_type, micexfg→named booleans/enums per the spec's naming, n_step_cl, n_step_vp, coll_level, ihabit_gm_random, CCNMAX, frac_dust, CRIC_RN_IMM, nucleation_halflife, aerosol chemistry params, debug/out settings) PLUS every `PARAM_ATMOS_PHY_MP_AMPS_bin` field (F2 §6: num_h_bins, nbin_h, iadvv, ini_aerosol_prf, l_* switches, l_gaxis_version, ...). Defaults = the Fortran defaults from F2 §6; `classmethod cloudlab()` returns the exact cloudlab AMPSTASK.F + run.conf values (F2 §5 + the run.conf namelist recorded in the M0 plan §Task 4 instructions doc). Derived params in `__post_init__` (dt_cl/dt_vp given dt at run time stay OUT — runtime, not config).
- `micexfg` mapping: keep BOTH representations — the named booleans AND a `micexfg_array() -> tuple[int, ...]` reconstructing the Fortran 20-slot array (validation against dumps needs the array).
- `core/index_maps.py`: the `par_amps` index constants (F2 §2) as IntEnums per group (LiquidPPV: rmt_q=1, rcon_q=2, rmat_q=3, rmas_q=4; IcePPV: imt_q..; AerosolPPV: ...) with 1-based Fortran values preserved and helpers `.py_idx` (0-based); `maxdims` true parameters (mxnmasscomp=8, mxnvol=2, mxntend=12, mxnaxis=5, mxnnonmc=7, mxnmasscomp_r=3) per F2 §3.

**Tests:** cloudlab() field-by-field equality with the F2-quoted values (write them as literals in the test, citing F2 line numbers in comments); micexfg_array() reproduces the AMPSTASK array exactly; index enums match par_amps literals.

- [ ] TDD steps; commit `feat(amps): AmpsConfig (AMPSTASK+namelist) and PPV index maps`.

---

### Task 5: LUT conversion + `core/lookup_tables.py`

**Files:**
- Create: `$WT/.../amps/codegen/convert_luts.py` (offline converter script), `$WT/.../amps/core/lookup_tables.py` (loaders), `$WT/.../amps/data/` (converted `.npz` artifacts, committed), `$WT/.../amps/data/README.md` (provenance)
- Test: `tests/amps/unit_tests/test_lookup_tables.py`

**Interfaces:**
- Converter: parses the AMPS_DATA files per the reader formats in F3 (RDCETB: drop_drop_Rey4.dat → drpdrp(201,201)+aux headers; hex/bbc/col_drop → (64,71)/(64,71)/(62,71); grp01/04/08 → (37,125)/(27,125)/(21,125); habit frequency files; stdnorm.dat → znorm(4,451)) into ONE `amps_luts.npz` with keys `<name>` and `<name>_aux` (xs, dx, ys, dy, nr, nc). Source dir CLI arg; defaults documented to `/Users/jcanton/projects/scale_amps/scale-rm/test/case/cloudlab/AMPS_DATA`. Run it ONCE and commit the npz (sizes are sub-MB per F3 — verify total < 5 MB before committing; if larger, STOP and report).
- Loaders: `load_luts() -> AmpsLuts` (frozen dataclass of arrays) via `importlib.resources` on the package `data/` (F5 §3 pattern); computed tables generated at load: osmotic (init_osmo_par), normal/inverse-normal, IGP spline knots (init_inherent_growth_par), Low–List breakup fragment tables (cal_breakfragment) — each transcribed VERBATIM from F3's quoted generator code.

**Tests:** shapes + aux headers match F3 literals; spot values: first data line of drop_drop_Rey4.dat (quoted verbatim in F3) appears at the expected indices; computed tables: IGP knot count == 23 (class_Group vap_igp_aux), breakup tables sized per bin count formula from F3; all arrays finite.

- [ ] TDD steps (converter tested against the real AMPS_DATA dir — it exists locally). Commit `feat(amps): LUT conversion pipeline and loaders (AMPS_DATA -> packaged npz)`.

---

### Task 6: `state.py` + `codegen/` skeleton — bin-field bundles and the split-operator generator

**Files:**
- Create: `$WT/.../amps/state.py`, `$WT/.../amps/codegen/generate.py`, `$WT/.../amps/codegen/templates.py`, `$WT/.../amps/core/generated/` (emitted, committed)
- Test: `tests/amps/unit_tests/test_state.py`, `tests/amps/unit_tests/test_codegen.py`

**Interfaces:**
- `state.py`: numpy-first bundles — `LiquidState`, `IceState`, `AerosolState`, `ThermoState` dataclasses holding arrays shaped `(nprops, nbins, ncat, npoints)` matching the Fortran qrpv/qipv/qapv layout EXACTLY (F4 §1 packing), plus `.to_fields()` / `.from_fields()` converting to/from per-(prop,bin) gt4py `(Cell, K)` field dicts named `f"{group}_{prop_name}_{bin:02d}"` (prop names from Task 4 index maps). Round-trip must be lossless.
- `codegen/generate.py`: `emit_module(name: str, build: Callable[..., str], **params) -> Path` — writes generated DSL source into `core/generated/{name}.py` with license header + `# GENERATED by codegen/generate.py — do not edit` banner, imports it (spike `load_generated_operator` pattern, F5 §6), returns module. `check_regenerated()` helper for a test that re-emits and diffs against the committed file (drift guard). Design constraint stated in module docstring: builders must emit SPLIT operators — the gate report's amendment — and the first real builder (M2's collection kernel) must chunk by destination bin; M1 ships only a trivial demonstration builder (`build_axpy_per_bin(nbins)` emitting one field_operator per bin group of ≤8 bins) to prove the pipeline + the drift test.
- `amps/utils.py`: `recursion_limit` context manager (own copy — no muphys import; carry-forward #6) + the spike-E RNG graduated as `counter_hash01` field_operator + numpy replica in `core/rng.py` with docstrings covering the int64 domain bounds (axes < 2^31, constants < 2^20) and the quadratic-round 2-to-1 collision caveat (carry-forwards #5/#8).

**Tests:** state round-trip lossless (random arrays, exact equality); generated demo module: emitted file identical on re-emit; embedded execution of the demo operator matches numpy; `core/rng.py` bit-exact vs its numpy replica on embedded (copy the spike-E assertion pattern at package scale).

- [ ] TDD steps; commit `feat(amps): state bundles, split-operator codegen skeleton, rng/utils graduation`.
- [ ] Follow-up in same task: remove the now-satisfiable entries from root `tach.toml [external] exclude` (numpy/gt4py definitely imported by src now; icon4py_common too via constants/type aliases) — run `uv run --group dev --frozen --isolated pre-commit run --all-files` to prove `tach check-external` passes without them (carry-forward #3). Commit `chore: drop temporary tach external excludes for amps`.

---

### Task 7: `core/packing.py` — SI↔CGS packing/unpacking (SCALE⇄AMPS state)

**Files:**
- Create: `$WT/.../amps/core/packing.py`
- Test: `tests/amps/unit_tests/test_packing.py`

**Interfaces (numpy, mirrors mp_amps driver blocks in F4):**
- `pack_scale_to_amps(...)` — the Z_LOOP_01 equivalences (F4 §1): factor_mxr1 = QDRY+QV mixing-ratio conversion, aerosol-mass-into-liquid-mass addition, ice total = QI + melt + aerosol, the ×0.001 non-mass-PPV conversions applied to EXACTLY the property indices F4 lists, Emoist bookkeeping (both l_no_ice_heat branches).
- `unpack_amps_to_scale(...)` — the RHOQ_t recipes (F4 §2): liquid aerosol-mass subtraction, ice melt+aerosol subtraction, ×0.001 reverse, CPtot_t/CVtot_t per Δq (l_no_ice_heat: CP_WATER-for-ice variant), l_gaxis_version=1 path only (2/3 documented as M5 stubs raising NotImplementedError).
- `moistthermo_mask(...)` — micptrv cloudy-mask criterion + thil/qtp diagnosis per F4 §3.

**Tests:** pack→unpack round-trip on synthetic states conserves each tracer identically when no physics ran (tendencies == 0); the ×0.001 factors hit only the F4-listed indices (construct a state with sentinel values per index and check); thil/qtp values recomputed independently in the test from F4's quoted formulas.

- [ ] TDD steps; commit `feat(amps): SI<->CGS packing (Z_LOOP_01/unpack recipes port)`.

---

### Task 8: `driver/ref_data.py` + `driver/box.py` — reference reader (pytest-ified) and box-driver skeleton

**Files:**
- Create: `$WT/.../amps/driver/ref_data.py`, `$WT/.../amps/driver/box.py`, `$WT/.../amps/driver/__init__.py`
- Test: `tests/amps/unit_tests/test_ref_data.py`

**Interfaces:**
- `ref_data.py`: port of `scale_amps/scripts/amps_dump_reader.py` (same record layouts v1, rank-aware keys, endian auto-detect) as a proper module: `read_dump_file(path) -> list[Record]` (typed dataclasses `MicroRecord`/`SedRecord`), `load_reference(npz_or_dir) -> RefDataset` with `.micro_pairs()` yielding (pre, post) per (rank, t, i, j) and `.sed_pairs()` per (rank, t, i, j, isn). Version asserts become `raise ValueError` (not assert — python -O safe; carry-forward #2). Keep `scripts/amps_dump_reader.py` in scale_amps authoritative for cluster-side conversion; this module reads BOTH raw dumps and the converted npz.
- `box.py`: skeleton only — `BoxCase` dataclass (column thermo profile + initial spectra + AmpsConfig + dt/steps) and `run_box(case) -> BoxResult` raising NotImplementedError with a docstring naming the M2 wiring; plus `case_from_micro_record(rec) -> BoxCase` implemented (turns one dumped pre-record into a runnable case — the M2 replay entry point).
- Sed-input derivation notes (carry-forward #7): document in `ref_data.py` docstring the derivations from the final review (waccv=-1, spdsfcv=0, momv=momz_col, thskinv from thetav/qvv, pgnd via Exner, dz1v from fz_col) with the F90 line references.

**Tests:** pytest-converted round-trip suite: reuse the byte-builders from scale_amps `scripts/test_amps_dump_reader.py` (copy into the test file; keep in sync note) covering little+big endian, two ranks, direct asserts on dt/kmicvm/k1b/k2b, pair-matching by (rank,t,i,j), and `case_from_micro_record` field mapping.

- [ ] TDD steps; commit `feat(amps): reference-data reader (typed, pytest) and box-driver skeleton`.

---

### Task 9: scale_amps housekeeping (Fortran dt kind + mapping doc start)

**Files (repo: /Users/jcanton/projects/scale_amps, branch cloudlab_port):**
- Modify: `scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90` — `AMPS_DUMP_micro`/`AMPS_DUMP_sed` dummy `dt`: `real(RP)` → `real(DP)` (carry-forward #1; two lines).
- Create: `docs/superpowers/facts/m1/fortran-mapping.md` — the mapping-doc seed: table of M1 Python modules → Fortran sources (from this plan's task interfaces), the not-ported-yet list, the known divergences so far (counter RNG vs LCG; l_gaxis_version 2/3 stubs; npz vs NetCDF).

- [ ] Edit + self-review diff (no local compile; visual check that only the two declarations changed); write doc; commit `fix(amps-dump): dt dummy kind DP; docs: start fortran mapping table`; push.

---

### Task 10: M1 wrap — full-suite run, docs, push

- [ ] Run the whole package: `uv run --group test --frozen pytest -n0 model/atmosphere/subgrid_scale_physics/amps/ -v` → all green, output pristine.
- [ ] `uv run --group dev --frozen --isolated pre-commit run --all-files` → clean.
- [ ] Update `$WT/.../amps/README.md`: module map + M1 status.
- [ ] Commit `docs(amps): M1 foundations status`; push `amps_microphysics`.

---

## Plan Self-Review (done at write time)

- Spec coverage: M1 items from spec §8 all mapped (constants T1, thermo T2, bin grid T3, AmpsConfig T4, LUT conversion/loading T5, state bundles + codegen T6, SI↔CGS adapters T7, ref-data reader + box skeleton T8). All nine M0 carry-forwards placed: #1→T9, #2→T8, #3→T6, #4→T2, #5→T6, #6→T6, #7→T8, #8→T6 (design note in codegen docstring) + M2 plan, #9→noted in T6 docstring for M2.
- Placeholder scan: none; formulas intentionally referenced to fact files (F1–F5) rather than duplicated — the fact files are committed ground truth and each task names its exact sections.
- Type consistency: index maps (T4) feed state bundles (T6) and packing (T7); `esat_lk` (T2) signature reused nowhere else yet; `AmpsConfig.cloudlab()` used by T8's `case_from_micro_record`.
- Known risk: fact-file gaps (an extractor may have truncated a formula) — mitigated by the NEEDS_CONTEXT rule in Global Constraints; controller resolves against the Fortran directly.
