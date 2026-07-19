# AMPS→icon4py Port — M2b (Warm-Phase: collection + collisional breakup) Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Warm-phase collision-coalescence (rain-rain stochastic collection) + Low-List collisional breakup, wired into the warm loop, plus porting `cal_breakfragment` (real fragment tables replacing the M1 placeholder) — and the M2a whole-branch-review hardening batch as opening tasks. Validated by conservation + physical-behavior tests; per-call dump replay pending user cluster runs.

**Architecture:** numpy reference implementations first (the correctness gate + StencilTest references). The coalescence bin-pair kernel is the O(nbins²) split-operator codegen target (M0 gate: chunk ≤8 bins/operator; numpy fallback sanctioned). CRITICAL basis fact (M2b facts H1): **coalescence operates per-VOLUME** (`gr%MS = qrpv·den`), opposite the per-mass activation/vapor-deposition ports — collection must scale state ×den before the kernel and ÷den after. This is the exact per-mass/per-volume trap that caused the M2a activation bug; the plan handles it explicitly (Task H0/6).

**Ground truth (committed on `cloudlab_port`), H1–H4:**
- H1 = `docs/superpowers/facts/m2b/coalescence-engine.md` (+ `docs/superpowers/facts/m2/coalescence.md` = G4, the earlier extraction; H1 fills G4's gaps) — coalescence engine, kernel, collision-eff LUT bilinear gather, scatter-to-target `collector_loop1`, per-volume density convention.
- H2 = `docs/superpowers/facts/m2b/collisional-breakup.md` — ibreak path, cal_breakup_dis_LL (Low-List), add_fragments_col_vec, P_breakup/Q_breakup2, bu_fd/bu_tmass runtime indexing.
- H3 = `docs/superpowers/facts/m2b/breakfragment-full-chain.md` — cal_breakfragment Cluster-1 port spec (cal_breakup_dis_LL, cal_Coalescence_Efficiency token==1, cal_sig_sf/cal_Hmusig, zbrent, getznorm2/cdfnor chain; minimal len/vtm — use core/liquid_diag._terminal_velocity, NOT diag_pq).
- H4 = `docs/superpowers/facts/m2b/codegen-and-plan-context.md` — codegen skeleton API, warm-phase module signatures, the B1/B2/D1-D5 hardening code sites, muphys collection precedent (none — AMPS collection is the first O(n²) bin-interaction kernel in icon4py).
- M2/M1 facts + fortran-mapping.md still apply. M2a whole-branch review findings: in the ledger `.superpowers/sdd/progress.md`.

**Tech Stack:** Python 3.12, gt4py==1.1.11, numpy, pytest; icon4py worktree `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics` (branch `amps_microphysics`, `$WT`), tip 17dac3d56.

## Global Constraints

- Transcribe physics VERBATIM from H1–H4/G4; on truncation the dispatch authorizes reading the named Fortran directly (quote into report). Never guess a formula.
- **Per-volume/per-mass discipline (the M2a lesson):** collection is per-volume (`con=q·den`, `mass=q·den`). Convert state→per-volume (×den) at the collection boundary, run the kernel in per-volume, convert back (÷den). Add a conservation test at REAL magnitudes (den≈1.2e-3 CGS) that fails if the ×den/÷den is wrong — do NOT let a fixture cancel it.
- **Test fixtures MUST use real magnitudes** (den≈1.2e-3, ptotv≠p00, per-cm³ number ~10²–10³) — degenerate values hid two systemic M2a bugs. Every conservation test asserts absolute balance, never a ratio that cancels the error.
- Units CGS internally, float64. Reuse M2a: core/state, core/liquid_diag (diag_pq_liquid/LiquidDiag incl. terminal_velocity), core/lookup_tables (AmpsLuts 7 collision-eff tables + BreakupFragmentTables guard to replace), config, implementations/warm_loop, core/constants/index_maps.
- License headers; ruff/mypy/tach/pre-commit clean; TDD; per-task report `$WT/.superpowers/sdd/m2b-task-N-report.md`; branch `amps_microphysics`; commit per task, push after each REVIEWED task (icon4py) + push scale_amps after any scale_amps task.
- Split-operator codegen (M0 gate, binding): the bin-pair kernel emits chunked operators (≤8 bins/op) via codegen/generate.py; numpy fallback allowed per-kernel, documented. gtfn compile tests may error if the local C++/Ninja toolchain is absent — that's infra, not a physics regression; numpy tests are the gate.

---

### Task 1 (HARDENING — opening): M2a whole-branch-review carry-forwards

**Files:** activation.py, vapor_deposition.py, liquid_diag.py, warm_loop.py, test_activation.py, test_vapor_deposition.py, test_warm_replay.py, test_warm_loop.py (all in `$WT/.../amps/`).

Address the M2a review items (details + file:lines in H4 and the ledger):
- **B1** — activation droplet placement is applied inline; Fortran defers it (`update_group_all(update_vapor=1)` AFTER `vapor_deposition(rain)`), so vapor-dep grows a pre-activation bin snapshot. EITHER defer activation's `rmt/rcon/rmat/rmas` bin increments until after `vapor_deposition_liquid` (mirroring the Fortran apply seam) OR, if that's disruptive to the current state-threading, document the divergence explicitly with the bounded-impact analysis and a test characterizing the overlap case. Choose deferral if clean (it's the faithful option and collision/M2b consumes the same seam); else document + flag for M2c/M3. Thermo threading is already correct — don't touch it.
- **B2** — delete the two orphan replay-stub tests (bare `pytest.skip` in test_activation.py, test_vapor_deposition.py) — the centralized `test_warm_replay.py` covers the composed path; OR implement+gate them like test_warm_replay.py. Delete is fine (state so in the report).
- **D1** — `np.select` without `default` in liquid_diag.py (_terminal_velocity, _den_aclen) silently returns 0 on NaN. Add `default=np.nan` (or an assert) so a bad value fails loud.
- **D2** — add one activation unit test at a non-p00 pressure (the current tests pin ptotv==p00, the ratio-1 masking value).
- **D3** — fix `_supersaturated_box_case` (test_warm_replay.py) moist_denv to a real ~1.2e-3 CGS magnitude.
- **D4** — note (not necessarily fix): n_step_cl>1 is never exercised with real physics — this task doesn't wire collision yet, so just add a ledger note; real coverage lands with Task 6.
- **D5** — it_cl==1 diag currency: Task 6 (collision) will consume `diag` at it_cl==1; ensure the refresh makes it current when collision lands (note here, enforce in Task 6).
- Cleanup: stale `/den` comment (activation.py ~2157) + dead `ActivationBoxState.den` field.

**Steps:** TDD where a behavior changes (B1, D1); mechanical for deletes/comments. Full suite green. Commit `fix(amps): M2a whole-branch-review hardening (B1 deferred-apply, B2/D1-D5, cleanup)`.

---

### Task 2: `core/collision_kernel.py` — collision efficiency + kernel (rain-rain)

**Files:** Create `$WT/.../amps/core/collision_kernel.py`; Test `test_collision_kernel.py`

**Interfaces:** `collision_efficiency(diag_i, diag_j, luts) -> E_c` (the drpdrp bilinear gather, H1 §3: axis index math from Reynolds `log10(Nre)` and radius ratio, clamped 2×2 stencil, the rrat<0 guard, ec_min/15 clamps) + `coalescence_efficiency(...)` (H1 cal_Coalescence_Efficiency token==1: E_coal, CKE, D_L/D_S/S_T/S_C) + `collision_kernel(diag_i, diag_j, con_j, dt) -> KC` (H1 §kernel: `E_c·(vtm_i−vtm_j)·A_c·con_j·dt`, per-volume con_j). Numpy, per-volume, CGS. Consumes LiquidDiag (terminal velocity, length, Nre).

**Tests:** bilinear gather vs hand-computed drpdrp interpolation for a known (Nre, rrat); rrat>1 reciprocal path; efficiency clamps; E_coal/CKE vs H1 formulas for a known drop pair; kernel value vs hand computation. Commit `feat(amps): collision efficiency + stochastic-collection kernel (rain-rain)`.

---

### Task 3: `core/coalescence.py` — the bin-pair collection engine (numpy)

**Files:** Create `$WT/.../amps/core/coalescence.py`; Test `test_coalescence.py`

**Interfaces:** `coalesce_rain(liquid_pv, diag, config, dt, luts) -> liquid_pv'` — the full rain-rain collector_loop1 (H1 §1.6): for each collector bin i (loop bounds/order per H1), each collected bin j, compute N_col = con_i·KC(i,j), the coalesced mass m_i+m_j, the destination-bin search + scatter (add_simple_vec/add_samebin_vec, H1 §2), number/mass conservation, cal_ratio_mass_col_vec property transfer, cal_needgive inter-bin borrowing. Operates on PER-VOLUME state (con/mass = q·den). Self-collection i==j and symmetry handled per H1. Breakup (ibreak) is Task 5 — leave a hook (ibreak=0 path only here).

**Tests:** two-bin collection conserves total number-reduces / mass-conserves (per-volume, real magnitudes); coalesced drops land in the correct destination bin (mass m_i+m_j → bin search); self-collection; a Golovin-kernel analytic check if tractable (constant kernel → known moment evolution); mass conservation to 1e-12 at den≈1.2e-3. Commit `feat(amps): rain-rain coalescence bin-pair engine (numpy, per-volume)`.

---

### Task 4: split-operator codegen for the collection kernel + gtfn feasibility

**Files:** Create `$WT/.../amps/codegen/collection_gen.py`, generated module(s) under `core/generated/`; Test `test_collection_codegen.py`

**Interfaces:** generate the bin-pair collection as SPLIT gt4py operators (≤8 destination bins/operator, per M0 gate + spike B) using the codegen skeleton (H4). The gather-formulated scatter (per destination bin, sum over source pairs) — reuse spike-A's tiled-table idiom for the collision-eff LUT if the kernel goes DSL. This task PROVES the codegen path for collection at nbins=40 (the spike-B gate said monolithic is NO-GO; split is the strategy). If gtfn compile is impractical locally (toolchain/time), the numpy engine (Task 3) is the sanctioned fallback — measure + document, don't block.

**Tests:** generated split operators reproduce the Task-3 numpy engine (embedded backend) to 1e-12; drift guard; per-chunk correctness. gtfn compile measured (record; may skip if toolchain absent). Commit `feat(amps): split-operator codegen for collection kernel`.

---

### Task 5: `core/breakup.py` — Low-List collisional breakup + runtime fragment consumers

**Files:** Create `$WT/.../amps/core/breakup.py`; Test `test_breakup.py`

**Interfaces:** the ibreak path in coalescence (H2): when a rain-rain collision breaks up, `cal_breakup_dis_LL` (Low-List filament/sheet/disk fragment distributions), add_fragments_col_vec, P_breakup/Q_breakup2, and the bu_fd/bu_tmass runtime indexing (the i1d_pair/kk index tying (bin_i,bin_j) to the precomputed table). Wire into coalesce_rain (Task 3) behind config.rain_collisional_breakup (micexfg 18, ON for cloudlab). Uses the fragment tables from Task 6.

**Tests:** a colliding drop pair above the breakup threshold produces fragments per Low-List distribution (vs H2 formulas); mass conserved across breakup; number increases (fragmentation); fragment-table indexing hits the right (i,j) entry. Commit `feat(amps): Low-List collisional breakup + fragment-table runtime consumers`.

---

### Task 6: port `cal_breakfragment` (real fragment tables) + validate-vs-dump instrumentation

**Files:** `$WT/.../amps/core/lookup_tables.py` (replace the BreakupFragmentTables placeholder with the real generator); Test `test_breakfragment.py`. Scale_amps: tiny dump instrumentation in `scale_atmos_phy_mp_amps.F90` (or mod_amps_lib) to dump bu_fd/bu_tmass at setup.

**Interfaces:** `make_breakup_fragment_tables(config, luts, bin_grid) -> BreakupFragmentTables` — the real cal_breakfragment (H3, Cluster-1 only): fixed air state (T=278.68, PT=850hPa, RH=100%), the bin-pair loop calling cal_Coalescence_Efficiency (token==1) + cal_breakup_dis_LL, filling bu_fd/bu_tmass; the index scalars jmin_bk/imin_bk/imax_bk/jmax_bk. CRITICAL (H3): compute liquid %len (geometric from bin mass) + %vtm DIRECTLY via `core/liquid_diag._terminal_velocity` (confirmed the right routine), NOT the monolithic diag_pq. Port cal_sig_sf/cal_Hmusig/zbrent/getznorm2/cdfnor (H3) as needed. Remove the allow_placeholder guard (real tables now).

**Instrumentation (scale_amps):** add bu_fd/bu_tmass dump at AMPS setup (they're setup-time constants) so the user gets them in existing cluster runs — the validation reference for the ported generator. Fold into the existing dump infra (Task 1/M0 style). Update the run instructions.

**Tests:** bu_tmass/bu_fd shapes + index scalars vs H3; spot fragment-distribution values recomputed from H3 formulas; the sig_sf/Hmusig Brent solves converge; a marker-gated test comparing the ported tables to dumped bu_fd/bu_tmass (skip-with-pointer if no dump). Commit (icon4py) `feat(amps): port cal_breakfragment (real Low-List fragment tables)`; commit (scale_amps) `feat(amps-dump): dump bu_fd/bu_tmass for breakfragment validation` + push.

---

### Task 7: wire collision+breakup into the warm loop

**Files:** `$WT/.../amps/implementations/warm_loop.py`; Test extend `test_warm_loop.py`

**Interfaces:** add the `_coalescence` hook to run_warm_micro_tendency's col_loop (per H1 §0 caller: `micexfg(2)==1 && flagp_r>0`, after the refresh, before repair('af_col'); it operates on per-volume state so the hook converts ×den in / ÷den out per the discipline). Wire breakup (micexfg 18) inside it. Enforce D5 (it_cl==1 diag currency — collision consumes diag). Add the n_step_cl>1 real-physics coverage (D4).

**Tests:** collision fires in the col_loop at the right point; a supersaturated multi-bin state run through the full loop (activation+vapor-dep+collision+breakup+repair) conserves total water+number appropriately (collection reduces number, conserves mass; breakup increases number); n_step_cl>1 exercised with real collision. Commit `feat(amps): wire coalescence+breakup into warm micro-tendency loop`.

---

### Task 8: M2b wrap — full suite, mapping-doc, README, push both repos

**Steps:** full amps suite green+pristine; pre-commit clean; update scale_amps `fortran-mapping.md` (M2b modules + remaining not-yet: sedimentation=M2c, all ice=M3); update package README (M2b status); commit+push both repos. Ledger + memory: M2b closed; note per-call replay + bu_fd/bu_tmass validation pending user dumps.

---

## Plan Self-Review (write-time)

- Scope: M2b = collection + breakup (the milestone def). Hardening (T1) opens per the M2a review's "before Task-3 collision" directive. Coalescence engine (T2 kernel, T3 numpy engine, T4 codegen), breakup (T5 runtime, T6 cal_breakfragment port + dump validation), wiring (T7), wrap (T8).
- Per-volume/per-mass: called out in Global Constraints + T3 + T7 with a real-magnitude conservation test — the explicit guard against repeating the M2a activation bug.
- Fixture realism: mandated in Global Constraints (the M2a lesson).
- Dependency order: hardening → kernel → engine → codegen → breakup runtime → fragment-table port → wire → wrap. T5 (runtime breakup) depends on T6 (tables) — note: T6 can produce a test-time table for T5's unit tests, or T5's tests use a small synthetic table; T6 delivers the real one. Sequence T6 before T5's full integration or let T5 use a synthetic fixture table (state in T5).
- Known risk: gtfn compile of the split collection kernel at nbins=40 (spike-B said split is the strategy; T4 measures, numpy fallback sanctioned). Per-call + bu_fd/bu_tmass dump validation pending user cluster runs (not a code blocker).
- Correction to dependency note: reorder so Task 6 (cal_breakfragment tables) comes before Task 5 (runtime breakup that consumes them), OR T5 uses a synthetic table fixture and T6 validates the real one — implementer picks the cleaner path, stated in the report.
