# AMPS Microphysics → icon4py (GT4Py) Port — Design Spec

**Date:** 2026-07-16
**Status:** approved design, pre-implementation
**Repos involved:**
- Source: `scale_amps` (this repo), branch `cloudlab` — AMPS Fortran in `contrib/AMPS/` (~72k lines), SCALE driver `scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90`.
- Target: `~/projects/icon4py`, new branch `amps_microphysics` off `main`, developed in worktree `.worktrees/amps_microphysics` (matches existing repo convention).

## 1. Goal

Port the AMPS bin (spectral, habit-predicting) microphysics scheme from Fortran to Python/GT4Py as a new icon4py package, with an aggressive clean redesign (the Fortran contains heavy duplication and ~28% dead commented-out code), validated against serialized reference data from the Fortran implementation, eventually upstreamable to C2SM/icon4py.

Entry point being ported: `ifc_cloud_micro` (`contrib/AMPS/mod_amps_lib.F90:39–326`) plus the SCALE-driver-side pieces that belong to the scheme (packing/unpacking, sedimentation, bin shift, negative adjustment, thermodynamic diagnosis).

## 2. Decision log (user-approved)

| Decision | Choice |
|---|---|
| Fortran baseline | **cloudlab branch** (includes INP-cap change in `mod_amps_core.F90:~17326`, widened `dbintendl(7,…)` diagnostics, `l_no_ice_heat` driver machinery). vs master the AMPS-core diff is tiny; the WBF/no-ice-heat work lives in the SCALE driver layer and defaults off. |
| Scope | **Full port: all switches, all processes used in `mod_amps_core.F90`** — including processes off in every repo config (hydrodynamic breakup, homogeneous freezing, contact nucleation, Hallett–Mossop, `act_type=2` activation). Staged so warm phase validates first. |
| Staging | Warm-phase first (matches cloudlab spin-up config with ice off), then ice. |
| Validation | Serialized per-call reference data around `ifc_cloud_micro`, tight rtol (~1e-10 warm phase). Instrumentation patch written here; runs executed by user on cluster. |
| Refactoring latitude | Aggressive clean redesign; Fortran↔Python mapping doc for traceability. |
| Backends | embedded + gtfn_cpu now; stay in DSL-compilable subset so gpu/dace come later for free. |
| Upstream | Eventually C2SM/icon4py: strict coding guidelines, CI, incremental PRs. |
| RNG (stochastic habit) | **Counter-based RNG** (hash of cell,k,bin,step) — parallel-safe, reproducible, statistically equivalent. Bit-for-bit with the Fortran Park–Miller LCG is impossible either way. |
| Sedimentation | **Full replication** including MOMZ/RHOU/RHOV/RHOE (momentum/energy) tendencies and surface precip fluxes; both Euler (`iadvv=1`, default) and PPM (`iadvv=2`) variants. |
| Bin architecture | **Approach A (codegen unrolling, pure DSL) as target, entered via spike (M0); per-kernel numpy fallback (Approach B) as sanctioned escape valve.** Sparse LOCAL BinDim (Approach C) rejected as highest framework risk. |
| Units | **CGS internally** (as Fortran AMPS), SI↔CGS adapters at the boundary where Fortran has its `0.001` factors. Keeps reference comparison tight. |
| Precision | Uniform float64 (SCALE default build makes `real(8)` accumulators same-kind as `PS` anyway; `scale_precision.F90:37-41`). |

## 3. Key recon facts the design rests on

- **True scope:** ~72k lines → ~45% comments (28 pts of that dead commented-out code), heavy duplication. Genuine unique compute ~14–18k lines; kernels to reproduce ~10–12k.
- **Processing shape:** strictly per-column; the scheme gathers cloudy levels of one column into a packed vector (L ≤ nz+1) per call (`scale_atmos_phy_mp_amps.F90:2004-2007`). Port replaces gather/scatter with dense masked (Cell, K) computation.
- **Time structure** (`cal_micro_tendency`, `class_Cloud_Micro.F90:788–1423`): outer collision substeps `n_step_cl` × {collection rain-rain / ice-ice / ice-rain, melting-shedding, hydrodynamic breakup ×2, Ice_Nucleation1, repair}, inner vapor substeps `n_step_vp` × {activation+vapor advance, vapor deposition rain/ice, Ice_Nucleation2, repair}. Production: `n_step_cl=1`, `n_step_vp=10`.
- **Real config source is `AMPSTASK.F`** (read by `read_AMPSTASK`, `mod_amps_utility.F90:1323-1421`), not the SCALE namelist. Three instances in repo (cloudlab / sheba / data/mp). All use `level_comp=7`, `act_type=1` (kc04dep), `ihabit_gm_random=1`. The port's `AmpsConfig` absorbs both AMPSTASK.F and `PARAM_ATMOS_PHY_MP_AMPS_bin`.
- **State layout:** liquid 4 props × 40 bins; ice 16 props × 20–40 bins (mass, rime/agg/crystal/melt/frozen/aerosol components, concentration, circumscribing volume, a/c/d axes³, gravity centers, extra crystals); aerosol 3 props × 1 bin × up to 4 categories. 493 SCALE tracers at defaults. `ifc_cloud_micro` updates state in place; SCALE forms tendencies by differencing.
- **Duplication to collapse:** `cal_aptact_var8_vec` vs `cal_aptact_var8_kc04dep` (~6k lines, ~90% identical incl. six contained solver routines — unify with flags); `*_vec`/`*_vec2` pairs (spheroid, bulk density, dep_mass, assign_Qp); immersion/KC04/homfreez shared skeleton; the 4×-repeated diagnosis preamble; six `cal_dmtend*` variants (only `cal_dmtend_scale` live).
- **GT4Py (1.1.11) feasibility:** no for/while in DSL; fixed-count unrolled iteration (muphys 6× Newton precedent) or host-side loop + mask (icon4py satad precedent); `scan_operator` fits sedimentation exactly (muphys precedent); `as_offset` gives runtime gather (no scatter — remaps must be gather-formulated); `gamma` builtin exists, `erf` missing (rational approximation); LOCAL dims are input-only; DaCe rejects a second horizontal dim — hence codegen unrolling.
- **Known Fortran quirks:** `dbintendl` declared `(3,2,mxnbin,L)` at `mod_amps_lib.F90:122` but `(7,…)` at callers (works via sequence association only) — port treats 7 as truth, flag upstream. Dead second `ifc_cloud_micro` call site behind `if (.false.)` (`mp_amps.F90:2975`) — ignored. History-diagnostic loop bug `k = KS, IE` at `mp_amps.F90:2940/2948`. `acc_amps.F90` is precision kinds, not OpenACC; the Fortran has no GPU support at all.

## 4. Package architecture

New uv-workspace package `model/atmosphere/subgrid_scale_physics/amps/`, registered in root `pyproject.toml` (workspace + pythonpath), `tach.toml` (depends only on `icon4py.model.common`), `noxfile.py`. Mirrors muphys conventions.

```
src/icon4py/model/atmosphere/subgrid_scale_physics/amps/
  core/                 # pure physics field_operators, no orchestration
    constants.py        # physical constants as float Enums, CGS
    thermo.py           # e_sat (analytic Murphy–Koop replaces estbar/esitbar LUTs),
                        #   diffusivity, conductivity, GTP; _scalar variants for scans
    bin_grid.py         # bin boundaries (binb(i+1)=a·binb(i)+b), split indices
    kohler.py           # critical/haze radius: analytic + fixed-iteration variants
    activation.py       # UNIFIED var8/kc04dep; svsteady coefficients, func_liqvap/icevap,
                        #   Brent → fixed-count masked bisection/Newton
    vapor_deposition.py # Chen–Lamb semidiscrete growth, growth-mode selection, inherent
                        #   growth ratio, capacitance, ventilation
    collection.py       # generic two-group collection engine (coalescence/aggregation/riming)
    breakup.py          # Low–List collisional + hydrodynamic
    melting.py          # melting/shedding, mv_ice2liq
    nucleation.py       # immersion (std + KC04), contact, homogeneous (Koop),
                        #   Hallett–Mossop, deposition/sorption
    habit.py            # Ice_Shape geometry, spheroid, bulk density, axis growth,
                        #   habit selection with counter-based RNG
    bin_remap.py        # cal_transbin / shift_bin, gather-formulated
    repair.py           # reality_check + budget repair (mass/concentration/volume closure)
    tendencies.py       # assign_tendency bookkeeping (11 process slots; slot 3 = activation/DHF,
                        #   per class_Mass_Bin.F90:80 — the core:20916 comment is stale)
    lookup_tables.py    # collision-efficiency/habit-frequency/osmotic/normal tables;
                        #   erf rational approximation
  codegen/              # bin-unrolling generator; emits real .py sources, committed to git
  implementations/      # composition + host-side substep orchestration
    warm_phase.py
    mixed_phase.py
  sedimentation.py      # Euler + PPM k-scans incl. momentum/energy tendencies
  state.py              # generated NamedTuple bin-field bundles; AmpsConfig
  driver/               # standalone column/box driver + NetCDF I/O
  data/                 # LUTs converted from scale-rm/test/case/cloudlab/AMPS_DATA → .npz/.nc
docs/fortran_mapping.md # routine → module table; dropped/changed/diverged list
tests/amps/             # unit_tests / stencil_tests / integration_tests (muphys pattern)
```

Codegen emits committed, inspectable `.py` files — reviewable upstream, debuggable, diffable when bin-count config changes trigger regeneration.

## 5. Data model

- **Bin state = generated NamedTuple bundles of `(Cell, K)` fields** per (property, bin). Property names derived from `par_amps` index maps.
- **Persistent between calls:** bin bundles + thermo fields only. All `Mass_Bin`/`Ice_Shape` diagnostic members (terminal velocity, capacitance, ventilation, semi-axes, …) are intermediates recomputed per substep (matches Fortran `diag_pq` refresh).
- **Config:** frozen dataclasses. `micexfg(1:19)` → named booleans/enums (index 19 = DHF activation flag `iflg_dhf` into `cal_aptact_var8_kc04dep`, per `class_Cloud_Micro.F90:1285`; index 20 genuinely unreferenced). Derived parameters in `__post_init__` (diffusion-granule style).
- **Global state eliminated:** `par_amps`/`maxdims`/`com_amps` → constants + config + init-time tables. Sequential LCG → counter-based hash RNG in DSL integer ops.
- **Sparsity:** dense masked computation with `where`; no compressed index lists.

## 6. Orchestration & kernel strategy

Host-side Python `AmpsScheme.run()` reproduces the operator-split skeleton exactly (substep loops in Python, each process one/few precompiled gt4py programs via `model_options.setup_program`):

```
for it_cl in range(n_step_cl):                     # dt_cl = dt/n_step_cl
    refresh_state()          # the 4×-duplicated Fortran preamble becomes one function
    collection(rain, rain);  collection(ice, ice);  collection(ice, rain)
    melting_shedding()
    hydrodyn_breakup(rain);  hydrodyn_breakup(ice)
    ice_nucleation_1()       # contact / immersion(std|KC04) / homogeneous / Hallett–Mossop
    repair('af_col'); update_groups()
    for it_vp in range(n_step_vp):                 # dt_vp = dt_cl/n_step_vp
        refresh_state()
        activation_and_vapor_advance()             # unified var8/kc04dep
        vapor_deposition(rain); vapor_deposition(ice)
        ice_nucleation_2()   # deposition/sorption modes
        repair('af_vap'); update_groups()
final_diagnosis(); terminal_velocities()
```

- **Collection engine:** one generator parameterized by (group pair, kernel LUT set, breakup flag). Dense generated bin-pair expressions; destination-bin deposits gather-formulated (per destination bin: masked sum over source pairs). LUT bilinear gathers via `as_offset` on table-as-field, or compile-time-folded per-bin-pair constants where indices are bin-static.
- **Solvers:** fixed-count masked iteration in-DSL where bounded and small; analytic variants preferred where Fortran offers them (`get_critrad_anal`). Budget repair (`NITER=300`, global `any()`) → host-side loop around a repair program with device mask and `.ndarray` convergence check (icon4py satad precedent).
- **Bin remap** (`cal_transbin`, `shift_bin`): gather formulation over the fixed bin count (Fortran is already dense i×ibx).
- **Sedimentation:** two k-scan implementations (Euler default, PPM opt-in) computing bin mass flux, energy transport, momentum tendencies, surface precip. CFL substep count = host-side scalar from a field reduction. PPM departure-cell logic reformulated to CFL-bounded fixed offsets (substepping bounds the offset, as in Fortran).
- **Kernel contract:** every kernel has a defined field-in/field-out signature; numpy fallback implementations (`core/*_np.py`) are drop-in, config-selected, tracked in the mapping doc, target count zero by upstream time.
- **Error paths:** debug writes/aborts → error-flag fields reduced host-side per phase. Error branches that also repair state keep the repair as physics.

## 7. Validation & testing

**Fortran instrumentation** (patch on `cloudlab` in scale_amps; user runs on cluster):

1. Per-call NetCDF dumps around `ifc_cloud_micro` call site #1: full in/out state (XC/XR/XS/XA, RV/DEN/PT/T/W, qtp/thil, RNG seeds, k-maps, dt, config echo). Subsampled (every Nth step / column subset) — full dumps every call would be enormous at 493 tracers.
2. Driver-level dumps around sedimentation (it lives outside `ifc_cloud_micro`): spectra + MOMZ_t/RHOU_t/RHOV_t/RHOE_t/precip in/out.
3. Optional per-process dumps inside `cal_micro_tendency` behind a flag, added where mismatches need bisecting.
4. Reference runs: cloudlab spin-up (ice off) → warm-phase; seeding run (ice on) → ice-phase; plus one ice run with `ihabit_gm_random=0` for tight ice comparison (counter-RNG ≠ Fortran LCG: stochastic-habit-affected fields validate statistically, everything else tightly vs the RNG-off run).

**Test pyramid** (muphys conventions):
- `StencilTest` unit tests per kernel with numpy references — no external data, CI-safe.
- Per-process integration tests vs per-process dumps (as added).
- Per-call replay tests: dumped inputs through the Python scheme, `assert_dallclose`, rtol ~1e-10 warm phase.
- Box-driver trajectory test over the spin-up hour; drift tracked.
- Data via `ICON4PY_TEST_DATA_PATH` locally first, CSCS-hosted tarball when upstreaming; `datatest` marker keeps CI green without data.

## 8. Milestones

- **M0 — Spike + scaffolding (gate).** Worktree/branch/package skeleton registered. Spike on embedded + gtfn_cpu: (a) compile/runtime of generated 40-bin collection kernel; (b) wide scan carry (40-field NamedTuple); (c) `as_offset` gather remap; (d) LUT-as-field vs analytic fit; (e) counter-RNG hash in DSL. Output: spike report fixing codegen design; go/no-go per pattern; fallback list seeded. Also: instrumentation patch + run instructions delivered.
- **M1 — Foundations.** constants, thermo, bin grid, `AmpsConfig`, LUT conversion/loading, state bundles + codegen layer, box driver skeleton, ref-data reader, SI↔CGS adapters.
- **M2 — Warm phase.** Unified activation, vapor deposition (liquid), collection (rain-rain) + collisional breakup, liquid repair, `diag_pq`(liquid), terminal velocity, warm sedimentation Euler+PPM incl. momentum/energy. **Validated per-call vs spin-up dumps.**
- **M3 — Ice growth & habit.** Ice thermo, Ice_Shape/habit geometry, counter-RNG habit selection, vapor deposition (ice) with growth modes, deposition nucleation, `diag_pq`(ice).
- **M4 — Ice interactions.** Aggregation, riming, melting/shedding, `mv_ice2liq`, immersion (both), contact, homogeneous, Hallett–Mossop, hydrodynamic breakup, full repair, ice sedimentation. **Validated vs seeding-run dumps.**
- **M5 — Hardening.** `act_type=2` unified path, remaining switches exercised, numpy-fallback burn-down, performance pass gtfn_cpu, dace smoke check, mapping doc complete.
- **M6 — Integration & upstreaming.** `AmpsComponent` + `State` adapter against post-PR#1301 `physics_interface` (Component/PhysicsState/ProcessTimeControl triple, muphys wiring as template); standalone-driver wiring; upstream PR slicing (earlier slices possible: foundations → warm → ice).

Each milestone gets its own implementation plan and PR(s); this spec is the single source of design truth.

## 9. Risks

| Risk | Mitigation |
|---|---|
| gtfn compile-time blow-up on 40-bin generated kernels | M0 spike measures before commitment; per-kernel numpy fallback contract |
| Scatter→gather remap reformulation changes numerics subtly | Standalone remap validation vs Fortran `cal_transbin` on synthetic inputs; per-process dumps to bisect |
| Ice validation vs stochastic habit selection | RNG-off Fortran rerun for tight compare; statistical tests for stochastic fields |
| PPM departure-index reformulation | CFL-bounded offsets + host substeps; Euler is the validated default |
| PR #1301 interface churn | Integration is M6 (last); component boundary matches muphys shape, which survives either protocol outcome |
| Reference-data volume | Subsampled dumps (configurable stride, column subset) |
| erf/normal-CDF accuracy | Rational approximation unit-tested against Fortran SPECFUN values |
| Replicating upstream AMPS bugs | Port the intent; document divergences in mapping doc; flag upstream (`dbintendl` 3-vs-7, `k=KS,IE` loop bug) |

## 10. Out of scope

- SCALE-side `l_no_ice_heat` experiment overlay (driver-level Emoist/cp bookkeeping) — a host-model coupling concern, not part of the scheme; revisit at M6 if wanted as a config option.
- The dead second `ifc_cloud_micro` call site (`if (.false.)` block).
- x/y PPM advection routines (`cal_flux_x/y`) — legacy, unused by SCALE.
- Fortran-side refactoring of `contrib/AMPS` itself.
- GPU performance tuning (design keeps it reachable; not a deliverable here).
