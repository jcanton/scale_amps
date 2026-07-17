# M0 Spike Report — GT4Py feasibility for the AMPS port

Date: 2026-07-17. Spikes executed 2026-07-16 (overnight run) in the icon4py
worktree `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics`,
branch `amps_microphysics`, commits `7ff3702e2..c6fcbe0b5`. Environment:
gt4py 1.1.11, Python 3.12, backends embedded + gtfn_cpu (`run_gtfn_cached`),
Apple Silicon macOS (compile times are indicative, not cluster-representative).

## Verdict table

Spike A measured four `as_offset` table-gather variants (not the two the
plan anticipated); Spike D measured three `esat` variants (analytic plus two
table forms). Rows below reflect what was actually run.

| Pattern | Spike | embedded | gtfn_cpu first/steady | Verdict |
|---|---|---|---|---|
| `as_offset` self-gather, per-point shift (`gather_shift`, bin remap) | A | first=0.003s, steady=2.02ms | first=4.298s, steady=0.71ms | GO |
| K-only table gather, direct return (`table_gather`) | A | FAILED — `DSLError` at `@gtx.field_operator` decoration (return-type dims deduced as K-only), backend-independent | same `DSLError`, occurs before backend dispatch | NO-GO |
| K-only table, broadcast to (Cell,K) before gather (`table_gather_broadcast`) | A | FAILED — runtime `IndexError` in the embedded gather kernel | FAILED — bare `NotImplementedError` in GTFN's `fuse_as_fieldop` ITIR transform | NO-GO (different failure per backend) |
| K-only table, pre-tiled to full (Cell,K) shape (`table_gather_tiled`) | A | first=0.002s, steady=2.27ms | first=0.181s, steady=0.16ms (warm-cache re-run; no independent cold-compile number was captured for this exact spike-A variant — see Spike D's `esat_table_tiled` below, same idiom, genuine cold first=2.36s) | GO — adopted LUT idiom |
| Generated collection kernel, nbins=8 | B | first=0.0s, steady=18.7ms | first=45.3s, steady=5.1ms | GO |
| Generated collection kernel, nbins=20 | B | first=0.1s, steady=116.4ms | first=340.6s, steady=73.3ms | GO with caveat — in-process compile-cache reuse only, no cross-process reuse (see Compile-time scaling) |
| Generated collection kernel, nbins=40, monolithic single operator | B | first=0.5s, steady=497.6ms | FAILED — `RecursionError` at Python's default recursion limit; first=2578.7s, steady=606.1ms once the limit is raised to 1e5 (muphys precedent) | NO-GO — exceeds the 900s (15 min) threshold either way |
| 80-field (2x40-bin) NamedTuple scan carry, genuine sequential dependency | C | full-scale (NCELLS=4096) impractical: original (degenerate-carry) generator did not complete in ~90-95 min (process disappeared, consistent with an unconfirmed OOM kill); fixed generator's full-scale embedded deliberately not re-attempted (opt-in only, expected to cost at least as much); small-scale (16/64/256 cells) confirms correctness, ~55-56ms/cell | first=186.9s, steady=115.9ms | GO on gtfn_cpu (just past the "clean" 120s boundary, inside go-with-caching band); embedded not viable at production scale |
| Analytic Murphy-Koop `esat_analytic` | D | first=0.02s, steady=10.22ms | first=5.28s, steady=12.76ms | GO — exact vs. the analytic formula by construction, slower than the table on gtfn_cpu |
| Table + linear interp, pre-tiled `esat_table_tiled` | D | first=0.02s, steady=15.95ms | first=2.36s, steady=2.08ms; max_rel_err_vs_analytic=3.90e-03 | GO — faster than analytic on gtfn_cpu; adopted for exact-Fortran-match validation (see Findings §3) |
| Table + linear interp, K-only (`esat_table_k_only`, non-portable) | D | FAILED — runtime `ValueError`, dims ordered wrong | unexpectedly succeeds: first=2.37s, steady=1.22ms | NO-GO — backend-inconsistent, not adopted despite gtfn success |
| Counter RNG (int64 LCG, no bitwise), 1 quadratic + 3 Lehmer rounds | E | first=0.00s, steady=2.43ms (run 1) / steady=2.95ms (run 2) | first=6.18s, steady=0.79ms (run 1) / first=6.30s, steady=2.35ms (run 2) | GO |

## Compile-time scaling (Spike B)

nbins=8: 45.3s (independent repro: 45.2s). nbins=20: 340.6s (repro: 340.9s).
nbins=40: `RecursionError` at Python's default recursion limit (crashes in
~13s, inside gt4py's own `MergeLet` ITIR transform pass — a compiler-pipeline
failure, not a slow compile) → 2578.7s (~43 minutes) once the limit is raised
to 1e5, matching the recursion-limit context manager icon4py's own `muphys`
package already uses for this identical gt4py limitation.

The measured scaling is markedly super-linear: ~7.5x compile time per step
(8→20: 7.52x for 2.5x more bins; 20→40: 7.57x for 2x more bins), consistent
with roughly `nbins^2.9` (near-cubic) growth — the `eve`-visitor-style ITIR
transform appears to cost worse than the tree size itself.

Extrapolation for the real coalescence kernel (this stand-in already
generates the full `nbins^2` pair-term structure — ~1600 gain terms plus
~1600 loss terms at nbins=40 — that a real coalescence process shares, but
the real kernel adds per-pair property-transfer terms, e.g. liquid/aerosol
mass partitioning across the two receiving bins, that this Golovin-kernel
stand-in does not model): estimating conservatively ~2-4x more expression
nodes per bin pair than the stand-in (an assessment, not an independently
measured number) and combining it with the near-cubic scaling already
observed, a monolithic 40-bin real coalescence operator would very plausibly
land well past the 3600s budget that got the stand-in's nbins=40 to complete,
or hit the recursion-limit crash at a bin count below 40. Against the plan's
thresholds (120s clean go / 120s-15min go-with-caching / over 15 min or a
crash is no-go): the monolithic 40-bin formulation is already a firm NO-GO
for the stand-in alone, and the real kernel's extra complexity only pushes
it further past the boundary — the strongest evidence for the split-operator
codegen requirement below.

## Findings & decisions

1. **Codegen design consequence** (Spike B): a single, fully-unrolled
   `field_operator` per collection/coalescence process does not scale to the
   full 40-bin count on gtfn_cpu — a hard `RecursionError` at Python's
   default recursion limit, and a 2578.7s (~43 min) compile even after
   raising the limit to 1e5, both against a 900s (15 min) budget. nbins=20
   compiles in 340.6s (inside the go-with-caching band) but only if the
   architecture ensures the compile happens once per long-lived process and
   the compiled program is reused in-process thereafter; nbins=8 is a clean
   go at 45.3s. The codegen layer must therefore emit **split operators**
   (e.g. one `field_operator` per destination bin, or banded groups of bins)
   rather than one monolithic 40-bin operator, and/or lean on the scan
   formulation (Spike C) where the physics is genuinely sequential; the
   numpy fallback stays sanctioned per kernel for any split that still fails
   to compile in budget.
2. **Scan carry convention**: gt4py's `scan_operator` emits the
   **post-update** carry — the value a scan step's body returns is
   simultaneously the carry fed into the next level's call and the value
   written to the current level's output field. Codegen for implicit
   sedimentation must generate scan bodies that compute and return this
   level's carry-dependent update in one step (matching the brief's
   already-correct `numpy_reference` convention); no pre/post flip is
   needed. Verified exactly (max abs diff 0.0) at small scale and via the
   gtfn_cpu correctness assert at full scale, and now protected by a
   committed `run_perturbation_check` regression test (a level-0
   perturbation must propagate into level 1's output) — added after review
   found the spike's original commit had measured a degenerate,
   non-sequential carry (see Compile-time scaling context: the corrected
   186.9s datum, not the original 82.7s, is this spike's real number).
3. **LUT policy**: a plain K-only table field gathered via `as_offset` and
   returned directly is a NO-GO — `DSLError` at `field_operator` decoration,
   backend-independent. Broadcasting the table to (Cell,K) before the gather
   is also a NO-GO, failing differently per backend (embedded: runtime
   `IndexError`; gtfn_cpu: compile-time `NotImplementedError`). A K-only
   table combined arithmetically with another gathered value before return
   decorates successfully but is backend-inconsistent (embedded: runtime
   `ValueError` on dims ordering; gtfn_cpu: unexpectedly succeeds) and is
   therefore not adopted despite the gtfn "success." The only idiom portable
   across both backends is the table **pre-tiled to full (Cell,K) shape**
   (memory-replicated across cells) — adopted as the port's LUT idiom, at an
   NCELLS-times memory cost judged negligible for small microphysics tables
   but worth flagging for any table that isn't tiny. Because Fortran's
   `estbar` **is** the table (not the analytic formula), exact reproduction
   of legacy Fortran runs requires the port to implement TABLE semantics
   (the tiled-gather idiom), not swap in an analytic formula: the tiled
   table differs from the true analytic Murphy-Koop value by up to
   max_rel_err_vs_analytic=3.90e-03 (~0.39%, intrinsic to the table's 1 K
   resolution and linear interpolation — not a port defect, since the DSL
   tiled implementation matches its own numpy replica of the table math at
   rtol=1e-12). **This revises the design spec's "analytic Murphy-Koop
   replaces `estbar` LUT" line**: analytic is not a drop-in replacement
   whenever a per-call validation tolerance tighter than roughly 1e-3 to
   1e-2 relative against legacy Fortran output is required — in that case
   the port must implement the table idiom. Analytic remains the right
   choice only where the validation target is physical accuracy against the
   true formula rather than bit-for-bit-ish Fortran matching — note the
   table is in fact the *faster* option on gtfn_cpu (2.08ms vs analytic's
   12.76ms steady-state), so choosing analytic there is a validation-fidelity
   tradeoff, not a performance one.
4. **RNG**: the construction that passed is one nonlinear (quadratic) mixing
   round, `x = (x*x + x + 1) mod M31`, inserted immediately after the
   initial affine combine of `(cell, k, bin, step)`, followed by the
   original 3 Park-Miller Lehmer rounds (A1=16807, A2=48271, A3=69621) — 4
   total mixing rounds, no A4 needed. Quality (numpy replica,
   `common.NCELLS x common.NLEV` sample): mean=0.5001 (bar: within 0.005 of
   0.5), var=0.0833 (bar: within 0.005 of 1/12=0.0833), lag1_cell=0.0001 and
   lag1_k=-0.0024 (bar: magnitude under 0.01), hist_dev=0.0003 (bar: under
   0.01) — all pass with 4x to ~70x margin. DSL output is bit-exact
   (`np.array_equal`, not just `np.allclose`) against the numpy replica on
   both backends. gtfn_cpu additionally requires the mixing constants to be
   passed as explicit `gtx.int64` scalar parameters, not referenced as
   module-level closure constants — the latter fails at ITIR-lowering time
   with `EveValueError` ("Symbols ... not found"), a different failure mode
   than the int32/int64 literal-promotion pitfall the plan anticipated.
5. **Failures/surprises verbatim**:
   - `table_gather` (Spike A, plain K-only, both backends, decoration-time):
     `DSLError`, "Annotated return type does not match deduced return type:
     annotation is 'Field[[Cell, K], float64]', got 'Field[[K], float64]'."
   - `table_gather_broadcast` (Spike A, embedded, runtime):
     `IndexError('arrays used as indices must be of integer (or boolean) type')`
   - `table_gather_broadcast` (Spike A, gtfn_cpu, compile-time, bare, no
     message): `NotImplementedError()` — raised in gt4py's own
     `fuse_as_fieldop.py` ITIR transform.
   - nbins=40 collection kernel (Spike B, gtfn_cpu, default recursion
     limit): `RecursionError: maximum recursion depth exceeded` — inside
     gt4py's own `MergeLet` ITIR transform pass.
   - `esat_table_k_only` (Spike D, embedded, runtime):
     `ValueError("Dimensions 'K[vertical], Cell[horizontal]' are not ordered
     correctly, expected 'Cell[horizontal], K[vertical]'.")`
   - counter RNG mixing constants as bare Python int (Spike E, decoration):
     `DSLError: Could not promote 'Field[[Cell], int64]' and 'int32' to
     common type in call to '*'.`
   - counter RNG mixing constants as `np.int64`-wrapped module constant
     (Spike E, gtfn_cpu, ITIR-lowering): `EveValueError: Symbols
     {SymbolRef('M31'), SymbolRef('C_K'), SymbolRef('C_CELL'),
     SymbolRef('ONE')} not found.`
   - counter RNG, naive 3-round-Lehmer-only construction (Spike E,
     statistical, not a DSL error): lag1_cell=-0.337, lag1_k=-0.496 against
     a 0.01 magnitude bar — fails by 33x to 50x; adding a 4th Lehmer round
     (the plan's own prescribed remedy for a failing quality assert) does
     not fix it either (lag1_cell=-0.387, lag1_k=-0.442 with 4 rounds; 5-7
     rounds also fail with no trend toward zero).
   - none reported as "FAILED" for Spike C beyond the embedded
     non-completion already covered in the verdict table (no Python
     exception was ever raised there — the process simply disappeared).

   Two further gt4py 1.1.11 API quirks (Spike A) are not outright spike
   failures but are load-bearing for every later spike's `as_offset` usage:
   a bare `Dimension` value in `offset_provider` is rejected on the compiled
   path (`offset_provider={}` is the correct form for a dynamic Cartesian
   `as_offset` shift, not `{"Koff": dims.KDim}`); and a `field_operator`
   parameter literally named `shift` collides with an internal GTFN lowering
   symbol (`EveValueError: Multiple definitions of symbol 'shift'`).

## Seeded fallback list (numpy-escape candidates going into M1+)

| Kernel | Reason | Trigger measured in |
|---|---|---|
| Monolithic 40-bin collection/coalescence operator (one `field_operator`, all bin pairs unrolled) | gtfn_cpu: `RecursionError` at default recursion limit; 2578.7s (~43 min, exceeds the 900s/15-min threshold) once the limit is raised to 1e5 | Spike B |
| Full-scale (NCELLS=4096) embedded execution of a wide (80-field) sequential scan carry | Does not complete in practical time (~90-95 min then the process is gone, consistent with an OOM kill, not independently confirmed via OS logs) | Spike C |
| Large (non-tiny) value-indexed lookup tables, where NCELLS-times tiled memory replication is prohibitive | The only portable `as_offset`-based LUT idiom is the pre-tiled (Cell,K) form; plain K-only and broadcast-before-gather forms are NO-GO on at least one backend each | Spike A, Spike D |

Escape hatch for row 1: split codegen (per-destination-bin or banded
operators) first; numpy fallback only if a split still exceeds the 900s
budget. Escape hatch for row 2: gtfn_cpu is the only backend exercised at
production scale; embedded is reserved for small-scale (16-256 cell)
correctness calibration. Escape hatch for row 3: numpy lookup for any table
whose tiled-memory cost is judged too large, on a case-by-case basis — no
such table was encountered in these spikes (all measured tables were small).

## Gate decision

M1 may start with Approach A (codegen unrolling) as planned: **yes**, with
the amendment that the collection-kernel codegen targets split operators
(per-destination-bin or banded groups of bins) rather than a single
monolithic 40-bin `field_operator`, since the monolithic form is a confirmed
NO-GO at gtfn_cpu compile time (Spike B: `RecursionError` at the default
recursion limit, 2578.7s/~43 min even with the limit raised). Deviations
from the spec required:

- Collection-kernel codegen must emit split (per-destination-bin or banded)
  operators, not one monolithic 40-bin operator; numpy fallback remains
  sanctioned per kernel for any split that still exceeds the 900s budget
  (Spike B).
- Any codegen path producing deep generated expression trees must wrap
  program compilation in a raised Python recursion limit (at least 1e5,
  matching the `muphys` driver precedent) — necessary but not sufficient on
  its own for the monolithic case (Spike B).
- LUT codegen must use the tiled (Cell,K) pre-replicated table-gather idiom;
  plain K-only and broadcast-before-gather forms are NO-GO (Spike A, D).
- The design spec's "analytic Murphy-Koop replaces `estbar` LUT" line (spec
  §4, `core/thermo.py`) is revised: the port must implement TABLE semantics
  (tiled table-gather) wherever validation requires reproduction of legacy
  Fortran output at a tolerance tighter than ~1e-3 to 1e-2 relative;
  analytic remains viable, and is faster on gtfn_cpu, only where the
  validation target is physical accuracy rather than bit-for-bit-ish Fortran
  matching (Spike D).
- Counter-based RNG construction must include at least one nonlinear
  (quadratic) mixing round in addition to the plan's Lehmer rounds, and must
  pass its mixing constants as explicit int64 scalar parameters rather than
  module-level closure constants (Spike E).
- Sedimentation scan codegen must treat the `scan_operator`'s returned
  NamedTuple as the post-update carry, matching this level's output field
  value directly — no pre/post flip needed (Spike C).
