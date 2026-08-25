# AMPS → icon4py port — handover

Last updated: 2026-08-25. Written for a cold pickup — a new session or a
collaborator who needs to know where the work stands, what's blocked, and
where to look, without re-deriving it.

## TL;DR

Porting the AMPS spectral-bin habit-predicting microphysics from Fortran
(SCALE) to Python/GT4Py (icon4py). **M0, M1, M2a done. M2b (warm collection +
collisional breakup) is code-complete and whole-branch-reviewed — no
correctness bugs — but warm per-call validation is BLOCKED** on a Fortran
dump-side number-PPV slot bug (handed to congchia). **M2c (sedimentation) is
held** per the standing "resolve validation before M2c" rule. Ball is with
congchia on one decisive `+1`-vs-`+2` LOG check.

## Repos, branches, paths

| What | Where |
|---|---|
| Fortran source (reference impl) | `scale_amps`, branch **`cloudlab_port`**, pushed to `origin` = `github.com/jcanton/scale_amps` (jcanton's fork) |
| AMPS Fortran entry point | `contrib/AMPS/`, entry `ifc_cloud_micro`; SCALE glue in `scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90` |
| Port target | `icon4py`, worktree `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics`, branch **`amps_microphysics`**, pushed to `origin` = `C2SM/icon4py` |
| Port package | `model/atmosphere/subgrid_scale_physics/amps/` (package `icon4py-atmosphere-amps`) |
| SDD ledger (per-task progress) | `<port>/.superpowers/sdd/progress.md` |
| Port doc facts | `scale_amps/docs/superpowers/facts/` — `m1/fortran-mapping.md` (module map + divergences), `verification/number-ppv-slot-mismatch.md` (the blocker) |

**Do NOT** commit on the `cloudlab` branch — use `cloudlab_port`. `init.conf` /
`run.conf` are congchia's namelists — leave untracked. Both repo branches are
currently fully pushed.

Who: congchia (Ong Chia Rui) authored the Fortran and generates the reference
dumps. jcanton owns the port.

## Milestone status

| M | Scope | Status |
|---|---|---|
| M0 | Feasibility / codegen gate | ✅ done |
| M1 | Thermo + framework, Fortran mapping | ✅ done |
| M2a | Activation solver primitives (Kohler, zbrent, svsteady) | ✅ done — 24 tests |
| **M2b** | **Warm collection (stochastic coalescence) + Low-List collisional breakup + `cal_breakfragment`** | **✅ code-complete + reviewed; ⛔ validation BLOCKED** |
| M2c | Sedimentation (Euler + PPM) — full replication | ⏸ HELD (resolve validation first) |
| M3 | Ice (habit prediction, shapes) | future |
| M4–M6 | remaining processes | future |

## THE BLOCKER — number-PPV QTRC slot off-by-one (Fortran side)

Full detail + evidence: `docs/superpowers/facts/verification/number-ppv-slot-mismatch.md`.

**Symptom:** warm reference dumps show the droplet-number / aerosol liquid PPVs
(`rcon_q`/`rmat_q`/`rmas_q`) as uniformly trace (~1e-9 #/g, ~14 orders below
physical ~1e5) across 4.76M points, while the **mass** field `rmt_q` validates
EXACT (5 sig figs) against the dump's own `qcvm`. The run LOG shows realistic
clouds.

**Ruled out:** (1) dump instrumentation — gather copies all `ipr` incl. `rcon`,
`AMPS_DUMP_micro` writes verbatim, all after `moistthermo2_scale`; mass reads
back exact so reader/PPV-axis mapping is correct. (2) pack/unpack round-trip —
since `moist_denv = DENS*factor_mxr1`, a no-op microphysics gives `RHOQ_t = 0`
exactly; number is conserved through the round-trip.

**Leading hypothesis:** two conventions in `scale_atmos_phy_mp_amps.F90`
disagree by one slot on where the per-bin number lives (stride `numberPPVL`):

- **Microphysics pack** (L1684-1693) + **unpack** (L2723-2739): number `rcon` at
  `I_QPPVL + 2 + numberPPVL*(ibr-1)` (block order `[rmat, rmas, rcon]`, canonical
  qrpv `rmt=1,rmat=2,rmas=3,rcon=4` from `mod_amps_utility.F90:1619-1636`).
- **Init `qhyd2qtrc`** (`lcon_index`, L4507) + **diagnostic/output `qtrc2qhyd`**
  (`liqConc_index = I_QPPVL+2-1`, L3979/L4751): number at `I_QPPVL + 1 + …`, with
  `lamt_index = I_QPPVL-1` (writes before the liquid-PPV region — the red flag).

**Why it reconciles everything:** LOG realistic (diagnostic reads +1 where init
wrote number); dump trace (micro pack reads +2, which init left ~0); mass real
(`I_QL`, unaffected); `rmas>rmt` and `rmat` tracks `rcon` (scrambled aerosol/
number slots = inert-residue signature).

**Decisive check for congchia:** in one known-cloudy `(k,i,j)`, LOG both
`QTRC(I_QPPVL+1)` and `QTRC(I_QPPVL+2)` (bin 1) beside `qcvm(k)` right before the
L2097 dump. Real number in +1 / ~0 in +2 confirms pack(+2) is out of step with
init/diag(+1) → reconcile to one sub-order; also confirm `num_h_moments(1)` /
`numberPPVL` stride and check `lamt_index=I_QPPVL-1` isn't corrupting the last
ice PPV.

**Impact on port:** NONE implied. The port reads its reference from the DUMP
(pack convention). Once Fortran is reconciled and re-dumped with real number,
warm per-call validation (and the held fast-follows below) can proceed against
trustworthy data.

## Held fast-follows (all wait on clean number data)

These cluster on degenerate-`con` handling and would calibrate to garbage if
fixed against the current trace-number dumps. Tracked in the SDD ledger:

- **F2** — `mean_mass <= binb[-1]` ceiling in coalescence is NOT literal Fortran
  (`cal_meanmass_vec` is floor-only) and is likely redundant in the real
  pipeline (`liquid_diag` pre-zeroes `mean_mass` before `coalesce_rain`; unit
  tests reach it only by bypassing the pipeline). Resolve: remove or document.
- **BLOWUP** — degenerate-`con` case: port gives 2.78e5 vs Fortran phase-2 dump
  8.2e-14 (repro rank22 / t2100 / i3,j3 / pt11 / bin3).
- **F1** — CLOSED (commit `7a2ddb86f`): counter-gate isolation test reworked to
  `mean_mass=1e-3` + sole occupied bin, mutation-verified.

Minors (tracked, no action): `breakfragment.py:789` `mrat` div0 matches
(also-unguarded) Fortran, never fires on cloudlab (136/136 nonzero); breakfragment
ground-truth test skipped pending `amps_dump_setup.bin`; `vapor_deposition.py:533`
unguarded errstate (pre-existing, harmless).

## Architecture facts worth not re-deriving

- **Units:** AMPS is CGS internally (g, cm, dyn/cm²), float64. SI↔CGS boundary is
  ThermoState CGS-canonical at 2 producers (`ptotv` ×10 Pa→dyn/cm², `moist_denv`
  ×1e-3 kg/m³→g/cm³).
- **Per-mass vs per-volume:** coalescence works per-VOLUME (×den in / ÷den out);
  qrpv mixing ratios are per-mass.
- **Operator split:** `cal_micro_tendency` = `col_loop`×`n_step_cl` /
  `vap_loop`×`n_step_vp`; `af_col`/`af_vap` repair steps.
- **Codegen (M0 gate):** GT4Py split-operator, ≤8 bins/op — monolithic 40-bin
  gtfn is a NO-GO. Collection rate-matrix is generated (25 ops, ≤64 pairs/op;
  1600 pairs reproduce numpy to 1e-12). LUTs tiled (Cell,K).
- **Reference dumps:** big-endian binary streams; `MicroRecord` wraps pre/post
  around `ifc_cloud_micro` with `trpv_thil`/`trpv_qtp`. Reader is
  `<port>/amps/driver/ref_data.py`; warm replay harness in `driver/box.py`
  (uses dumped `trpv_thil`/`trpv_qtp` when present).
- **Validation status is honest:** "active-bin fidelity" for M2b is VACUOUS —
  max `rcon` 1e-9 << 1e-6 active threshold, so no bin qualifies active. This is
  NOT a passing reproduction; it is blocked on the number-PPV fix.

## Next steps (in order)

1. **congchia:** run the decisive `+1`-vs-`+2` LOG check; reconcile the number-PPV
   slot; re-dump warm reference with real number.
2. **Port:** re-run warm per-call validation (M2a + M2b) against clean dumps.
3. **Port:** resolve held fast-follows (F2 ceiling necessity, BLOWUP degenerate-con)
   against trustworthy data.
4. **User checkpoint:** M2c go/no-go, then implement sedimentation (Euler + PPM,
   full replication).

## Workflow conventions

- SDD: fresh implementer + reviewer subagents per task; ledger at
  `.superpowers/sdd/progress.md`; whole-branch review per milestone.
- Checkpoint with the user per sub-milestone; work autonomously between.
- Commit footer per `~/.claude/CLAUDE.md`: `🤖 Written by an agent on behalf of
  @jcanton` as the last line (replaces default Co-Authored-By / Generated-with
  trailers). CSCS CI trigger comments stay bare (`cscs-ci run …`).
- Push port progress to `C2SM/icon4py` origin; push scale_amps docs to jcanton's
  fork so congchia can see them.
