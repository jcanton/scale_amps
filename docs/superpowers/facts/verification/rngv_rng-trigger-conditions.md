I have all the evidence needed. Here is the report.

---

# AMPS Stochastic Ice-Habit Selection (`ihabit_gm_random=1`): When It Can Change Results

Repo: `/Users/jcanton/projects/scale_amps`, branch `cloudlab_port` (confirmed checked out).

## (a) Every RNG-gated code path

Only two routines actually call the RNG (`rand2_ty`), and both do so **only inside an `if(ihabit_gm_random.eq.1)` block**. When the flag is `0`, no `rand2_ty` call executes at all — the `else` branches use pure max-frequency logic.

### Path 1 — `cal_growth_mode_inl_vec` (planar/columnar/polycrystal + ppo/ros sub-habit)
`contrib/AMPS/mod_amps_utility.F90:750`, RNG loop at **lines 778–784** (two draws `rnum1`, `rnum2` **per bin×level, unconditionally** — no ice-presence or ambiguity guard):
```
778   if(ihabit_gm_random.eq.1) then
781     do in=1,nbin*L
782       call rand2_ty(rnum1(in),rdsd)
783       call rand2_ty(rnum2(in),rdsd)
```
Decision (lines 810–831): `x1<=pol_frq → igm=1` (then `x2<ppo_frq & Tc>-45 → 5`; `x2<ros_frq & Tc<=-45 → 4`); `elseif x1<=pol_frq+pla_frq → igm=2`; `else igm=3`. Deterministic counterpart (lines 846–867) picks `argmax(pol,pla,col)`.

Called from:
- `class_AirGroup.F90:253` — `ini_AirGroup`, `nbin=1`; result → `a%TV(n)%nuc_gmode` (line 259). **Nucleation growth mode.**
- `class_AirGroup.F90:331` — `update_AirGroup`, `nbin=1`; result → `a%TV(n)%nuc_gmode` (line 335). **Nucleation growth mode, refreshed every step.**
- `class_Group.F90:9104` — `cal_growth_mode_vec` (full `N_bin×L`); result `igm(i,n)` → `g%IS(i,n)%growth_mode` **only** for bins with `icond1==0`, size `<10 µm`, and `ag%TV(n)%T<253.16` (lines 9113–9118). **Depositional growth mode of newly-nucleated small crystals.**

### Path 2 — `cal_growth_mode_hex_inl_vec` (planar-vs-columnar only)
`contrib/AMPS/mod_amps_utility.F90:986`, RNG loop at **lines 1012–1017** (one draw `rnum` **per bin×level, unconditionally**):
```
1012  if(ihabit_gm_random.eq.1) then
1015    do in=1,nbin*L
1016      call rand2_ty(rnum(in),rdsd)
```
Decision (lines 1032–1036): `x<=pla_frq → igm=1 else 2`. Deterministic (lines 1049–1054): `pla_frq>col_frq → 1 else 2`.

Called from:
- `mod_amps_core.F90:16622` — depositional hex growth (`N_bin×L`); `igm(i,n)` consumed at **line 16643** (`gamma_d=gamma(n,igm(i,n))`) **only when `ag%TV(n)%T<=253.16` AND `g%IS(i,n)%init_growth==1`** (lines 16640–16641).
- `mod_amps_core.F90:17162` — nucleation hex (`nbin=1`); `igm(n)` consumed at **line 17190** only when `ag%TV(n)%T<=253.16` (line 17188).

## (b) Warmest T at which stochastic vs deterministic can diverge → **−20.0 °C (253.16 K)**

### Table geometry (from reader + dims)
Reader `RDCETB`, `mod_amps_utility.F90:189–223`: each `*_frq.dat` = 1 comment line, then `nrow ncol`, then `nrow` rows each of `ncol` values, read as `frq(i,j)`, `i=1..nrow` (rows), `j=1..ncol` (cols). All five files are **`51 101`** (verified: 51 data lines × 101 fields each).

Axis mapping (identical in `cal_growth_mode_inl_vec:764–793`, `_hex:999–1026`, `get_growth_mode:680–682`):
- **Temperature axis = rows `itmp`**: `sttmp=-70`, `dtmp=1`, `itmp=nint((Tc+70)/1)+1`. So `itmp=1 ↔ Tc=-70 °C` (cold end), `itmp=51 ↔ Tc=-20 °C` (**warm end**), 1 °C resolution. `itmp` is clamped to `[1,51]`, so **any T warmer than −20 °C reuses row 51**.
- **Supersaturation axis = cols `isi`**: `stsi=0`, `dsi=0.01`, `isi=nint(|si|/0.01)+1`, range si=0.0…1.0, 101 points.

Normalization (`mod_amps_utility.F90:226–242`): pol/pla/col clamped `≥0` then rescaled so `pol+pla+col=1` per cell; `ros`/`ppo` only clamped `≥0` (used as conditional sub-probabilities).

### Table evidence
Numerically scanning the actual data files in `.../AMPS_DATA/collision_data/`:
- **Every one of the 51 rows (−70…−20 °C) has ≥2 nonzero of {pol,pla,col} in all 101 columns.** There is *no* warm sub-region where a single habit dominates deterministically.
- The **warmest row, `itmp=51` = −20.0 °C**, is ambiguous: normalized values e.g. `si=0.00 → pol 0.608 / pla 0.387 / col 0.005`; `si=1.00 → pol 0.812 / pla 0.034 / col 0.154` (three nonzero). For the hex routine (pla vs col), 69/101 columns at −20 °C have both nonzero.

So max-frequency pick (e.g. `pol` at −20 °C, si=0) versus a random draw (≈38.7 % of the time picks `pla`, 0.5 % picks `col`) genuinely diverge **at the table's warm edge −20 °C**, and would diverge for any warmer T too because of row clamping.

### The binding threshold
Divergence is capped not by the table (which stays ambiguous to its warm edge) but by the **consumer guards**, all keyed to **253.16 K = −20.0 °C**:
- `class_Group.F90:9115` — random `growth_mode` applied only if `T<253.16`.
- `mod_amps_core.F90:16640` — random hex `igm` applied only if `T<=253.16`.
- `mod_amps_core.F90:17188` / `:17224` — nucleation `nuc_gmode` applied only if `T<=253.16`; warmer, `growth_mode` is set from the deterministic axis-length ratio (lines 17227–17231).

**Conclusion: the warmest temperature at which stochastic vs deterministic habit selection can change the simulation is exactly −20.0 °C (253.16 K)** — set jointly by (i) the frq table remaining multi-habit up to its warm edge −20 °C, and (ii) the `T ≤ 253.16 K` application guards. Above −20 °C the dice is still rolled but the result is discarded.

## (c) Nucleation-time habit assignment — **YES, it uses the RNG**

The nucleation growth mode `nuc_gmode` is produced by the RNG, not by a separate deterministic path:
- Set from the `cal_growth_mode_inl_vec` dice roll in `ini_AirGroup` (`class_AirGroup.F90:253→259`) and re-set every step in `update_AirGroup` (`:331→335`).
- Consumed at nucleation: `growth_mode(n)=ag%TV(n)%nuc_gmode` (`mod_amps_core.F90:17225`) → `assign_Qpini_v3_vec(...,growth_mode)` (`:17244`); also `assign_Qp_v3p_vec` reads `nuc_gmode=ag%TV(n)%nuc_gmode` (`:22730`) and `assign_Qp_v3p(...,nuc_gmode,...)` (`:21752`, `:22522`) branches habit on `nuc_gmode` values 1/4/5 (`:22575–22644`).

So the dice roll happens **both at nucleation (via `nuc_gmode`) and during depositional growth (via `IS%growth_mode` / hex `igm`)** — but in every case the stochastic result is only *applied* when `T ≤ −20 °C`.

## (d) RNG call-site guards

The dice roll is **unconditional per bin per substep** whenever `ihabit_gm_random==1` — it is **not** gated on ice presence or on growth-mode ambiguity:
- `mod_amps_utility.F90:781` and `:1015` loop `do in=1,nbin*L` with no per-bin conditional; `rand2_ty` fires for every bin×level slot regardless of whether that bin holds ice or whether the (T,si) cell is single-habit.
- Ambiguity/ice filtering happens only *downstream*, when the returned `igm`/`nuc_gmode` is written into a crystal's state (the `icond1==0`, size `<10 µm`, `T<253.16`/`<=253.16` guards cited above).

Consequence worth flagging: because `rand2_ty` advances the shared generator state `rdsd` a fixed 2 (or 1) times per bin×level *only* in the `==1` branch and *zero* times in the `==0` branch, toggling `ihabit_gm_random` shifts the entire downstream RNG stream — so results can differ even in cells that are single-habit, via RNG-stream desynchronization, independent of the habit-choice divergence analyzed above.