# Adversarial Verification: "cloudlab result is not affected by ihabit_gm_random because temperature is not low enough"

## 1. Load-bearing facts, independently checked

**Fact A — RNG fires only inside the flag, and nowhere else in the codebase consumes the stream.** Verified: `rand2_ty` is called at exactly three sites, all inside `if(ihabit_gm_random.eq.1)` blocks — `contrib/AMPS/mod_amps_utility.F90:782,783` (`cal_growth_mode_inl_vec`) and `:1016` (`cal_growth_mode_hex_inl_vec`). A repo-wide grep for `rand2_ty`/`rdsd` shows **no other consumer** of the generator state (the `mod_amps_lib.F90:1861` local is a diagnostic path). This is stronger than report 1 established and it neutralizes report 1's own desync caveat (see §2).

**Fact B — every consumption of the RNG-drawn result is guarded at 253.16 K (−20.0 °C).** Verified at all six sinks:
- `contrib/AMPS/class_Group.F90:9115-9117` — `growth_mode=igm` only if size <10 µm **and** `T<253.16`; the else branch (`:9124-9139`) derives `growth_mode` from the crystal's existing `habit`, never from `igm` (checked — closes a gap none of the reports explicitly closed).
- `contrib/AMPS/mod_amps_core.F90:16640-16643` — hex `igm` used only if `T<=253.16 .and. init_growth==1`; and `init_growth=1` itself is only set under the `T<253.16` guard above.
- `mod_amps_core.F90:17188-17197` — nucleation `gamma_d` from `nuc_gmode`/`igm` only if `T<=253.16`; warm branch uses `gamma(1,n)` with the in-code comment "gamma is set up to be the same in T>-20C for planar and columnar growths."
- `mod_amps_core.F90:17224-17231` — `growth_mode(n)=nuc_gmode` only if `T<=253.16`; else deterministic axis-length ratio.
- `mod_amps_core.F90:22568-22569` (`assign_Qp_v3p`) and `:22740-22741` (vec variant) — the `nuc_gmode` habit branch is gated on `(level==3|5|7) .and. T<253.16_PS`; the warm else branch ignores `nuc_gmode`.
- Deterministic `nuc_gmode` overwrites exist (`class_Group.F90:9578`, `class_Thermo_Var.F90:303,364` via `get_growth_mode_max`) — these only *reduce* RNG influence.

**Fact C — warmest multi-habit table row is −20 °C and clamps warmer.** Verified: `pol/pla/col_frq.dat` are 51×101; data row 51 (warm edge, `itmp=nint((Tc+70))+1` clamped to 51) has pol≈0.639, pla≈0.406, col≈0.0052 at column 1 — genuinely multi-habit, so a draw *would* diverge from argmax there. But per Fact B it is discarded above −20 °C. Report 1's numbers reproduce.

**Fact D — cloudlab never gets close to −20 °C.** Verified from `scale-rm/test/case/cloudlab/env.txt`: domain (KMAX=60 × DZ=20 m = 1200 m top, `scripts/run.conf:20,34`) min T = **266.09 K (−7.06 °C) at 600 m**, max 272.16 K. Reaching 253.16 K would need ~13 K of cooling within the 1 h spin-up + 20 min restart (`run.conf:40`, `restart_run.conf:40`) — physically implausible under the configured weak forcings. *This is the one link in the chain that is a physical judgment, not a code proof.*

**Fact E — seeding injects aerosol, not ice, with zero RNG involvement.** Verified: `scale-rm/src/user/mod_user.F90:506-524` writes only `QS_MP+I_QPPVA+...` (aerosol mass / number / soluble mass), category `ica==2` only, one k-level, first J-row; no `I_QPPVI` ice PPV, no call into any habit routine. Ice enters solely via deposition nucleation later (`mod_amps_core.F90:17150-17155`: needs `T<TF`, `Si≥0`, dust category `con≥1e-22 & mass≥1e-22`), whose habit assignment is the T-guarded path in Fact B. So "seeded ice bypasses the RNG" is false in the direct sense (nothing to bypass — no ice is injected) and moot in the indirect sense (nucleated ice habit is deterministic above −20 °C).

**Fact F — the flag is actually 1 in the cloudlab config**: `scale-rm/test/case/cloudlab/AMPSTASK.F:85` (`ihabit_gm_random = 1`). No report stated this explicitly. Because of Facts A+B, with the domain warm, flag=1 vs flag=0 should be **bit-for-bit identical**: draws are made but discarded, and the perturbed `rdsd` state feeds nothing else.

## 2. Contradictions / gaps between the reports

1. **Report 1's RNG-desync caveat is wrong in its implication.** It warns that toggling the flag "shifts the entire downstream RNG stream — so results can differ even in cells that are single-habit... independent of the habit-choice divergence." Since `rand2_ty` has no consumers outside the two T-guarded habit routines (Fact A), desync can only manifest through those same guarded sinks. In cloudlab (never ≤ −20 °C) desync has exactly zero effect. This internally contradicts report 1's own −20 °C conclusion and, taken at face value, would refute the colleague's claim — it should be struck.
2. **Report 3's ICE_FF claim is numerically wrong.** It says the sigmoid is "~0 above ≈ −9 °C." The formula (`mod_amps_core.F90:17324`: `ICE_FF = 0.97 − 0.97/(1+exp(−0.88(T−263.95)))`) gives **0.485 at −9.2 °C**, ≈0.14 at cloudlab's coldest level (−7 °C), ≈0.01 at −4 °C. So seeded INP do activate a non-trivial fraction at cloudlab temperatures. Doesn't touch the RNG verdict, but the sub-claim as written would falsely suggest the seeding run makes almost no ice.
3. **Report 2's radiative-drift bound is estimate, not evidence.** "Does not drop below ~262 K" is hand-waved (flux-divergence arithmetic). Fine as a bound argument for a 13 K margin, but the honest closure is checking min(T) in the actual history output — none of the reports did.
4. Minor confirmations, no conflicts: spin-up seeding window overlap (window `[18:00:00, 19:01:30)`, `run.conf:213-218`, spin-up ends 18:00:00) — reports 2 and 3 agree and configs verify; restart window = first 60 s (`restart_run.conf:214-219`), RATE=10, z-limit 300 m all verified. Reports 1 and 3 agree that RNG draws occur (and are wasted) even in the warm/ice-free spin-up — consistent with Fact A.

## 3. Verdict

**SUPPORTED, with two conditions.** The claim's mechanism is exactly right and code-verified: every path by which `ihabit_gm_random` can influence model state is guarded by `T ≤ 253.16 K` (−20 °C), the RNG stream has no other consumers, seeding injects only INP aerosol, and the cloudlab domain starts at −7.1 °C minimum. Conditions: (i) it holds only while no grid cell ever reaches −20 °C — a ~13 K margin that is physically near-certain but should be confirmed once from history-file min(T); (ii) it holds for this code state — any future enabling of another `rdsd` consumer (or of immersion/contact/homogeneous freezing paths that might use it) reopens the desync channel.

## 4. Vetted list of other run-affecting quirks

| Quirk | Status | Confidence |
|---|---|---|
| `MP_do_precipitation=.false.` is mandatory for AMPS; AMPS-internal sedimentation still runs (`l_sediment=.true.`) — sed dumps unaffected | configs verified (`run.conf:137,151`; `restart_run.conf:138,152`); driver-abort detail not re-checked | **High** |
| Aerosols fully prognostic (`l_fix_aerosols=.false.`, `run.conf:150`) — port must not regenerate from initial profiles | config verified | **High** |
| Spin-up is warm only because IN category is empty (`N_ap_ini = 317.0, 0.0, 0.0, 0.0`, `AMPSTASK.F:198`) and drop-freezing modes off (`micexfg`, `AMPSTASK.F:286`); ice code paths (incl. RNG draws) still execute and no-op | verified, incl. nucleation gate `mod_amps_core.F90:17150-17155` | **High** — but its "ICE_FF≈0 above −9 °C" sub-claim is **wrong** (0.485 at −9 °C) |
| `RESTART_SKIP_READING_AMPS_ICE=.true.` (`restart_run.conf:62`) absorbs 20→40 ice-bin change by skipping all ice tracers; skipped QTRC slices formally uninitialized | flag + bin change verified (`run.conf:145` vs `restart_run.conf:146`); uninit caveat plausible, not re-traced | **High** (flag) / **Medium** (uninit hazard) |
| 13 spin-up + 5 seeding dump times (not 12); seeding records ~2× larger (nbi=40) | arithmetic consistent with verified durations; MP-call-count code path not re-traced | **Medium** |
| Hard-coded `/cluster/scratch/congchia/...` data paths (`AMPSTASK.F:203,209,212`) break any other user's run | verified | **High** |
| Dump files clobbered between runs (`status='replace'`, fixed names) | not re-verified | **Medium** |
| Dumps absent in `l_no_ice_heat=.true.` code path (both configs set `.false.`, `run.conf:152`) | config verified; call-site claim not re-traced | **Medium-High** |
| Spin-up seeding window opens exactly at spin-up end; tendency computed then discarded (restart written first) | window verified (`run.conf:209-218`); write-order claim not re-traced | **Medium** |

Key file:line evidence (verified this session): `mod_amps_utility.F90:778-784,1012-1017` (RNG-gated draws); `class_Group.F90:9115-9139` (guarded write + habit-based else); `mod_amps_core.F90:16640-16652, 17188-17197, 17224-17231, 22568-22569, 22740-22741` (T-guarded sinks); `mod_amps_core.F90:17324` (ICE_FF); `AMPSTASK.F:85,198,203,286`; `mod_user.F90:506-524` (aerosol-only seeding); `env.txt` min T 266.09 K @ 600 m; `scripts/run.conf` / `scripts/restart_run.conf` flags as tabulated.