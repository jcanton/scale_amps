# Fortran ↔ icon4py `amps` Package Mapping (M1)

Seed of the Fortran↔Python mapping doc for the `amps` port. Ground truth for
this table is the committed M1 task reports in the icon4py worktree
(`/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/.superpowers/sdd/m1-task-{3..8}-report.md`;
tasks 1–2's reports were not retained in that directory — sourced instead
from `core/constants.py`/`core/thermo.py`'s own module docstrings, which cite
the same fact files) and the fact extracts under
`docs/superpowers/facts/m1/` (F1–F5, see the M1 plan
`docs/superpowers/plans/2026-07-17-amps-port-m1-foundations.md`). Update this
file as later milestones port more of `mod_amps_core.F90`.

icon4py package root: `model/atmosphere/subgrid_scale_physics/amps/src/icon4py/model/atmosphere/subgrid_scale_physics/amps/`
(paths below are relative to this root unless stated otherwise).

## Module → Fortran source map

| Python module | Fortran source(s) | Fact file | Notes |
|---|---|---|---|
| `core/constants.py` (`AmpsConst`) | `contrib/AMPS/mod_amps_const.F90` (base constants); `contrib/AMPS/acc_amps.F90` (precomputed geometry constants) | F1 §1–§2 (`constants-thermo.md`) | Derived constants (`MR`, `Racp`, `Rdvchiarui`, `M_a`) computed from expressions, not hardcoded. CGS units throughout. |
| `core/thermo.py` — `make_esat_tables`/`esat_lk` | `mod_amps_utility.F90` `QSPARM2` (table build, Murphy–Koop 2005); table accessor semantics | F1 §3a, §3c | Fortran-exact truncate-then-clamp index arithmetic (`I=MAX(1,MIN(int(T)-163,149))`) preserved verbatim (M0 carry-forward #4). |
| `core/thermo.py` — `esat_analytic` | Lowe & Ficke analytic formula | F1 §3b | Alternate to the table accessor; NOT the default path used elsewhere in the port. |
| `core/thermo.py` — `t_from_esat_lk` | reverse-search-from-table routine | F1 §3d | ±5-entry window + full-scan fallback preserved. |
| `core/thermo.py` — thermo coefficient fns (D_v, k_a, d_vis, GTP/MGTP, sig_wa) | `class_Thermo_Var.F90` coefficient routines | F1 §3e+ | One function per Fortran routine. |
| `core/thermo.py` — `cal_til`/`cal_thetail`/`diag_t` equivalents | `class_Thermo_Var.F90` θ_il closure (`cal_til`, `cal_thetail`, `diag_t`) | F1 §5 | T-from-θ_il+condensate-loadings closure, numpy. |
| `core/bin_grid.py` | `binmicrosetup_scale` (liquid + ice bin-boundary construction) | F2 §1 (`bins-config-indexmaps.md`) | Liquid 40/80-bin, ice 10/20/40-bin; `binb[i+1]=a_b*binb[i]+b_b` recurrence, per-bin (not scalar) `a_b`/`b_b`. Fortran "last-loop-iteration overwrites `binbr(nbr+1)`" quirk deliberately preserved (rel=1e-12, not bit-exact, on that one boundary). |
| `config.py` (`AmpsConfig`) | `contrib/AMPS/AMPSTASK.F` (`/AMPS_param/` namelist, read via `read_AMPSTASK`); `PARAM_ATMOS_PHY_MP_AMPS_bin` namelist (scale-side) | F2 §4–§6 | `cloudlab()`/`cloudlab_seeding()` classmethods pin exact cloudlab/restart `run.conf` values. `micexfg` kept both as named booleans/ints and as `micexfg_array()`. |
| `core/index_maps.py` | `contrib/AMPS/par_amps.F90` (9 index-space groups: `LiquidPPV`, `IcePPV`, `AerosolPPV`, `*MassIndex`, `*MassRatioIndex`, `IceAxisIndex`); `contrib/AMPS/maxdims.F90` (6 true `parameter`s) | F2 §2–§3 | `FortranIndex.py_idx` gives 0-based access while preserving 1-based Fortran `.value`. |
| `codegen/convert_luts.py` (offline converter, not packaged) | `RDCETB`/`RDSTTB` readers (collision-efficiency, habit-frequency, diagnostic-map, `stdnorm.dat` files) | F3 (`lut-files.md`) | Parses `AMPS_DATA/*` into `amps_luts_collision.npz` + `amps_luts_misc.npz` (split to satisfy `check-added-large-files`, 500 KB/file default). |
| `core/lookup_tables.py` — `load_luts` | packaged npz loader (`importlib.resources`) | F3 §4 | `AmpsLuts` frozen dataclass mirrors Fortran `sequence` derived types. |
| `core/lookup_tables.py` — `init_osmo_par` equivalent (`_osm_ammsul`/`_osm_sodchl`) | `mod_amps_utility.F90` `init_osmo_par`/`osm_ammsul`/`osm_sodchl` | F3 §5.1 | x-grid for `osm_ammsul` transcribed verbatim from the real Fortran read (`mod_amps_utility.F90:12996`) after Addendum 1 corrected an earlier inferred-uniform-grid guess; `osm_sodchl`'s x-grid likewise verbatim. |
| `core/lookup_tables.py` — `init_normal_lut`/`init_inv_normal_lut` | DCDFLIB `cdfnor`-based normal/inverse-normal CDF LUT generators | F3 §6 | Ported via `scipy.stats.norm.sf`/`.ppf` (F3 explicitly sanctions this as the intended equivalent). |
| `core/lookup_tables.py` — `init_inherent_growth_par` | IGP spline-knot generator | F3 §5 | Verbatim, no gaps; 23 knots. |
| `core/lookup_tables.py` — `breakup_fragment_table_sizes`/`make_breakup_fragment_tables` | `cal_breakfragment` (`mod_amps_lib.F90:1831-2017`) — SIZING ONLY | F3 §5.5 | See "Not ported yet" — fragment-mass/-count physics fill is NOT implemented; `bu_fd`/`bu_tmass` are zero-filled and hard-guarded (`NotImplementedError` unless `allow_placeholder=True`, plus a runtime `is_placeholder` flag on the returned dataclass). |
| `core/state.py` (`LiquidState`/`IceState`/`AerosolState`) | `qrpv`/`qipv`/`qapv` array layouts (`(nprops, nbins, ncat, npoints)`), PPV property axis order per `par_amps.F90` | F2 §2, F4 §1.1 | Uses `core/index_maps.py`'s `LiquidPPV`/`IcePPV`/`AerosolPPV` for the property axis order (1-based Fortran `.value` ↔ axis index). |
| `core/state.py` (`ThermoState`) | `Z_LOOP_01` thermo block (`ptotv, tv, thv, piv, pbv, moist_denv, qvv, thetav, wbv, momv`) | F4 §1.1 | New, package-local `ThermoProp` enum — no Fortran `qXpv` array backs it. `ncat != 1` raises `NotImplementedError` (matches Fortran's own "assumed to be 1" scope). |
| `core/packing.py` — `pack_scale_to_amps` | `Z_LOOP_01`, `scale_atmos_phy_mp_amps.F90:1625-1828` | F4 §1 | `factor_mxr1` mixing-ratio conversion; aerosol-into-drop-mass addition (`rmt_q`, `imt_q`); the ×0.001 unit factor applied to concentration/axis slots only, never mass slots. |
| `core/packing.py` — `unpack_amps_to_scale` | tendency block `scale_atmos_phy_mp_amps.F90:2676-3010`; energy block `:2305-2352` | F4 §2 | Every `RHOQ_t`/`CPtot_t`/`CVtot_t`/`RHOE_t` recipe. `l_gaxis_version=1` fully implemented (incl. v1-only `l_axis_limit` clip-overwrite); `l_gaxis_version ∈ {2,3}` raises `NotImplementedError` in both pack and unpack. |
| `core/packing.py` — `moistthermo_mask` | `moistthermo2_scale`, `scale_atmos_phy_mp_amps.F90:1053-1117` | F4 §3.2 | qc/qr/qi bin partition against `RRLMTB`/`RILMTB=1e-22`; closed-form `qtp`/`thp`. SS3.3's iterative T-refinement is NOT ported (out of task scope, not needed for the closed-form path used). |
| `core/rng.py` | none — counter-based design, graduated from `spikes/spike_e_counter_rng.py` | F5 §6 (`icon4py-m1-conventions.md`) | See "Known divergences" — this is NOT a port of the Fortran LCG. |
| `driver/ref_data.py` | `scale-amps/scripts/amps_dump_reader.py` (stays authoritative for cluster-side `.bin`→npz conversion); record layouts produced by `AMPS_DUMP_micro`/`AMPS_DUMP_sed` (`scale_atmos_phy_mp_amps.F90`, ~lines 5283-5372) | — | Typed port (`MicroRecord`/`SedRecord` dataclasses) of the same binary record layout; endian auto-detect and version checks ported verbatim (version-mismatch now `raise ValueError`, not `assert` — M0 carry-forward #2). |
| `driver/box.py` | box-driver skeleton (no single Fortran routine — orchestration scaffold) | — | `run_box()` raises `NotImplementedError` (M2 wiring target); `case_from_micro_record()` is implemented and reconstructs a `BoxCase` from one `MicroRecord`. |

## Not ported yet (M2+)

All physics process routines in `contrib/AMPS/mod_amps_core.F90` remain unported:

- Activation (`cal_aptact_var8_kc04dep` and friends)
- Vapor deposition
- Coalescence (`cal_Coalescence_Efficiency` — only the liquid-liquid branch is *read*, for the breakup call-tree inventory below; not implemented)
- Breakup (`cal_breakup_dis_LL`, `cal_sig_sf`, `cal_Hmusig`, `zbrent`, `getznorm2` — full call-tree inventoried in `m1-task-5-report.md` Addendum 1, ~1285 lines across 6 files; BLOCKED as a physics-port task, not a LUT-conversion task)
- Melting
- Nucleation (`Ice_Nucleation1`, homogeneous/immersion freezing dispatch)
- Habit (habit-frequency table consumers)
- `diag_pq` (only the liquid/`phase==1` branch was read for the breakup call-tree, ~200-270 of 648 lines; not implemented)
- Repair
- Sedimentation (physics; `driver/ref_data.py` reads the *dumped* sed records but does not replay the sedimentation math)
- The real `cal_breakfragment` fill: `bu_fd`/`bu_tmass` are zero-initialized placeholders in `core/lookup_tables.py`, guarded against silent use (`NotImplementedError` unless `allow_placeholder=True`). Cloudlab has rain collisional breakup ON (`micexfg[18]`) — M2 needs the real values. Call-tree inventory (routine, file:lines, why) is in `m1-task-5-report.md` for the follow-up task to start from.

## Known divergences (intentional, documented)

| Area | Fortran | Python (M1) | Why / status |
|---|---|---|---|
| RNG (habit selection) | Park-Miller-style LCG | Counter-based hash: affine combine of (cell, level, bin, step) → 1 quadratic mixing round → 3 Park-Miller Lehmer rounds, mod `M31=2**31-1` | Bit-for-bit match with the Fortran LCG is impossible by construction (different algorithm family); statistically equivalent (mean~0.5, var~1/12, lag-1 corr <0.01, flat 16-bin histogram — validated in `spike_e_counter_rng.py`). Verification report (`docs/superpowers/facts/verification/rngv_verdict.md`) finds this irrelevant for cloudlab. |
| `l_gaxis_version` | 3 methods: v1 ("ORIGINAL"), v2, v3 ("METHOD 1"/"METHOD 2") | Only v1 implemented; v2/v3 raise `NotImplementedError` in both `pack_scale_to_amps` and `unpack_amps_to_scale` | M5 stub — F4 SS1.3/SS2.4 only summarize v2/v3 tersely; cloudlab config uses v1. |
| Reference data format | Raw binary streams (`AMPS_DUMP_micro`/`AMPS_DUMP_sed`, per-rank `.bin` files) | Binary streams read directly (`driver/ref_data.py`), or pre-aggregated into `.npz` via `scale-amps/scripts/amps_dump_reader.py` | No NetCDF; the Fortran-side dump format is intentionally kept as flat unformatted Fortran writes for cluster-side simplicity, converted off-cluster. |
| `esat` (saturation vapor pressure) | Table-lookup semantics (`QSPARM2` + truncate-then-clamp accessor) is what production Fortran actually calls | `esat_lk` (table, default used elsewhere) + `esat_analytic` (Lowe & Ficke, available but NOT the default path) | Both ported; table-lookup semantics preserved bit-exact (index arithmetic, M0 carry-forward #4) since that's what the Fortran call sites use. |
| Breakup fragment tables | `cal_breakfragment` computes real `bu_fd`/`bu_tmass` per bin pair | Zero-filled placeholders, hard-guarded (`NotImplementedError` unless `allow_placeholder=True`, plus runtime `is_placeholder` flag) | Physics fill deferred to M2 (see "Not ported yet"); guard prevents silent consumption of zeros as real data. |
| Osmotic-coefficient x-grid (`osm_ammsul`) | Non-uniform grid (0.1 step 0.0–1.0, 0.2 step 1.0–2.0, 0.5 step 2.0–5.5) | Initially inferred as uniform `np.linspace` (wrong); corrected to the verbatim non-uniform grid after reading `mod_amps_utility.F90:12996` directly | Resolved in Task 5 Addendum 1 — flagged here since the original submission's guess was wrong and worth remembering as a lesson (fact-file paraphrase ≠ verbatim quote). |
| KC04 deposition-freezing `flg1` selector (`cal_aptact_var8_kc04dep`'s `func_icevap_vec`, M2a) | `flg1=x(n)` assigned OUTSIDE the `do m=1,Lbx` box loop it later gates (`mod_amps_core.F90:8919-8925`), so `n` is whatever index was left over from a PRECEDING loop — every box's classical-nucleation (`iflg_dep==1`) gate ends up keyed off ONE arbitrary, stale box's supersaturation, not its own | `flg1 = x_arr` (per-box array) in `core/activation.py::func_icevap_vec` | Reads as an accidental Fortran scalar-aliasing artifact (leftover loop variable reused post-loop), not intended physics — no plausible reason box A's gate should depend on box B's `x`. Ported the evident per-box intent instead, per this project's divergence policy (port intent, document divergence, flag upstream). **Upstream-bug candidate — not verified against a scale_amps maintainer; worth confirming before relying on either behavior.** |

## Carry-forward items referenced above

- **#1** (this task): `AMPS_DUMP_micro`/`AMPS_DUMP_sed` `dt` dummy kind `real(RP)`→`real(DP)` — fixed in `scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90`.
- **#2**: dump-reader version-mismatch checks are `raise ValueError`, not `assert` (Python `-O` safe) — done in `driver/ref_data.py` (M1 Task 8).
- **#4**: esat table accessor's truncate-then-clamp index arithmetic preserved bit-exact — done in `core/thermo.py` (M1 Task 2).
