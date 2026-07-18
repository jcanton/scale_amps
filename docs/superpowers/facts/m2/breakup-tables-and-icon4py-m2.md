I have everything needed. Writing up the deliverable.

---

# M2 Planning Extraction — Two Parts

---

# PART A — Breakup-fragment table generation (`cal_breakfragment`): port vs. dump

## A.1 Full verbatim of `cal_breakfragment` (`contrib/AMPS/mod_amps_lib.F90:1831–2017`)

```fortran
    subroutine cal_breakfragment(NRBIN,estbar,esitbar)
      use scale_prc, only: &
         PRC_abort
      use com_amps
      use class_AirGroup, only: &
         AirGroup, &
         make_AirGroup_2
      use class_Group, only: &
         Group, &
         make_group, &
         ini_group_MP
      use mod_amps_utility, only: &
         random_genvar
      use mod_amps_core, only: &
         cal_Coalescence_Efficiency, &
         cal_breakup_dis_LL, &
         diag_pq
      implicit none

  integer,intent(in) :: NRBIN

  ! saturation vapor lookup table
!tmp  real  :: estbar(150),esitbar(111)
  real(DS),intent(in)  :: estbar(150),esitbar(111)

  ! steady growth problem: phase over which saturation is calculated.
  ! phase 1: water, 2: ice
  integer    :: phase

  ! random generator vars (not used)
  type (random_genvar) :: rdsd

  ! steady growth problem: relative humidity
  real(PS)    :: RH

  ! define the vapor pressure
  real(PS)              :: e

  ! definition for reading command line
  character(len=20)      :: chbuf
  real(PS)               :: temp ! ambient temperature in Celcius


  real(PS),dimension(1) :: rv,den,pt,T,W,z


  type (Group)               :: liquid
  type (AirGroup)            :: steady


  ! coalescence efficiency
  real (PS)                 :: E_coal
  real (PS) :: D,D_L,D_S,S_T,S_C,DS_S,CKE
  ! low-diameter cut off related to the resolution of the experiments (cm)
  real(DS),parameter :: D_0=0.01

  integer                                :: i, j,n



  ! message from reality_check
!tmp  integer,pointer,dimension(:)   :: mes_rc
  integer,dimension(1)   :: mes_rc
  integer,dimension(1)   :: ID,JD,KD

  ! optional status check
  integer  :: var_Status,em

  ! coefficient: sqrt(2*pi), pi/6.0
  real(DS),parameter  :: coef1=2.506628275,coef2=0.523598776

  ! maximum number of elements for lookup table, bu_tmass and bu_fd
  integer :: i1d_pair_max,kk_max

  integer :: ierr

!tp  allocate(mes_rc(1))


  !if ( IsMaster ) then
  !  write(fid_alog,*) "Now calculating collisional breakup fragments."
  !end if

  bu_fd=0.0_PS
  bu_tmass=0.0_PS

  T=278.6795
  PT=850.0e+2_PS
  W=0.0
  phase=1
  RH=100.0_PS

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! construct liquid hydrometeor
!tmp  liquid = make_Group ( 1, 0.0, 1, NRBIN, &
!tmp       srat_r, sadd_r, minmass_r, 1,1.0,1.0,1.0,1.0,binbr)
  call make_Group (liquid,1, 0.0_PS, 1, NRBIN, &
       srat_r, sadd_r, minmass_r, 1,1.0_PS,den_aps(1),den_api(1),eps_ap(1),binbr)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! construt thermo_var object
  if ( IsMaster .and. debug ) then
     write(fid_alog,*) "ck estbar",estbar(1:10)
  end if

  steady=make_AirGroup_2(1,estbar,esitbar,RV,DEN,PT,T,W,phase,RH)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! initialize the mass and con of liquid
  call ini_group_MP(liquid)


  ! diagnose the property
  call diag_pq(liquid,steady,1,em, mes_rc,ID,JD,KD,rdsd,ihabit_gm_random &
              ,eps_ap(1),nu_aps,phi_aps,m_aps)


  ! find bin that has the minimum size for possible breakup.
  do i=1,NRBIN
     if(liquid%MS(i,1)%len>=D_0) then
        jmin_bk=i
        exit
     end if
  end do

  imin_bk=jmin_bk+1
  imax_bk=NRBIN
  jmax_bk=NRBIN-1

  i1d_pair_max=(imax_bk-1)-jmin_bk+1+(imax_bk-imin_bk)*(1+imax_bk-imin_bk)/2
  kk_max=i1d_pair_max*liquid%N_BIN

  if ( IsMaster .and. debug ) then
    write(fid_alog,*) "i1d_pair_max,kk_max",i1d_pair_max,kk_max

    write(fid_alog,*) "size of bu_tmass, bu_fd",size(bu_tmass),size(bu_fd)/2
  end if

  if(size(bu_tmass)<i1d_pair_max) then
    LOG_ERROR("cal_breakfragment",*) "Error.  Increase the bu_tmass array"
    call PRC_abort
  endif
  if(size(bu_fd)/2<kk_max) then
    LOG_ERROR("cal_breakfragment",*) "Error.  Increase the bu_fd array"
    call PRC_abort
  endif



  do i=imin_bk,imax_bk
     do j=jmin_bk,i-1
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! determine the coalescence efficiency
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call cal_Coalescence_Efficiency(liquid,i,liquid,j,1,steady%TV(1),E_coal,&
             D_L,D_S,S_T,S_C,DS_S,CKE)

        if(CKE<=1.0e-20) cycle
!!c             call cal_breakup_dis(liquid,i,j,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd,&
!!c                                  1.0,D_L,D_S,S_T,S_C,dS_S,CKE)
        call cal_breakup_dis_LL(liquid,i,j,imin_bk,jmin_bk,bu_tmass,bu_fd,&
                                D_L,D_S,S_T,S_C,CKE)


     end do

     if(D_L*1.0e+2>0.4.and.D_L*1.0e+2<0.5) then
        call print_bufd(i,liquid,"check_bkup0p4.dat")
     elseif(D_L*1.0e+2>0.3.and.D_L*1.0e+2<0.4) then
        call print_bufd(i,liquid,"check_bkup0p3.dat")
     end if

  end do

!!c  write(fid_alog,*) "cal_break bu_fd"
!!c  write(fid_alog,*) bu_fd
!!c  write(fid_alog,*) "cal_break bu_tmass"
!!c  write(fid_alog,*) bu_tmass


  ! +++ deallocate memory used for the Group Object +++
!tmp  call delete_Group(liquid)
!tmp  call delete_AirGroup (steady)

!tmp  deallocate(mes_rc)
end subroutine cal_breakfragment
```

## A.2 Transitive call chain — line-count table

Two distinct dependency clusters emerge. **Cluster 1** is the actual breakup math (small, pure, deterministic). **Cluster 2** is the `Group`/`AirGroup` construction + `diag_pq` property diagnosis (large, drags in the full microphysics diagnostic stack).

| Routine | File | Lines | Size | Role in table generation | Cluster |
|---|---|---|---|---|---|
| `cal_breakfragment` | mod_amps_lib.F90 | 1831–2017 | **187** | driver: sets fixed state, loops bin-pairs `(i,j)`, fills `bu_fd`/`bu_tmass` | — |
| `cal_breakup_dis_LL` | mod_amps_core.F90 | 12019–12447 | **429** | THE core: Low-List (1982) filament/sheet/disk fragment distributions | 1 |
| `cal_Coalescence_Efficiency` | mod_amps_core.F90 | 11796–12017 | **222** | per-pair `E_coal`, `CKE`, `S_T`, `S_C`, `D_L`, `D_S` (only token==1 branch used) | 1 |
| `cal_sig_sf` | mod_amps_core.F90 | 27254–27302 | **49** | normal-dist σ via Brent root-find | 1 |
| `cal_Hmusig` | mod_amps_core.F90 | 27304–27354 | **51** | lognormal H/μ/σ via Brent root-find | 1 |
| `zbrent` | mod_amps_utility.F90 | 12778–12907 | **130** | Brent root-finder with 3 internal objective functions | 1 |
| `getznorm2` | mod_amps_utility.F90 | 369–375 | **7** | standard-normal CDF wrapper | 1 |
| `cdfnor` (+ CDFLIB helpers `cumnor`, `dinvnr`, `stvaln`, `devlpl`…) | mod_amps_utility.F90 | 10033–10230 | **~198 (+~200 helpers)** | special-function impl behind `getznorm2` | 1 |
| **Cluster-1 subtotal** | | | **~1470** | pure breakup math + normal CDF | |
| `make_Group` | class_Group.F90 | 201–577 | **377** | builds bin grid (`binb`, `MS%len`) — only grid part needed | 2 |
| `ini_group_MP` | class_Group.F90 | 10991–11063 | **73** | zero-init mass/concentration | 2 |
| `make_AirGroup_2` | class_AirGroup.F90 | 98–128 | **31** | wraps thermo var + est/esit LUTs | 2 |
| `make_thermo_var3_2` | class_Thermo_Var.F90 | 309–369 | **60** | one thermo state → `sig_wa` etc. | 2 |
| `diag_pq` | mod_amps_core.F90 | 12552–13199 | **648** | diagnoses `%len`, `%vtm`, densities; dispatches ↓ | 2 |
| `cal_meanmass_vec` | class_Mass_Bin.F90 | 1840–1849 | 9 | mean mass per bin | 2 |
| `cal_den_aclen_vec` | class_Group.F90 | 11222–11300 | 78 | density / a-c lengths | 2 |
| `cal_terminal_vel_vec` | class_Group.F90 | 7988–8289 | **301** | **terminal velocity `%vtm`** (needed by `E_coal`) | 2 |
| `cal_ventilation_coef_vec` | class_Group.F90 | 7575–7986 | **411** | ventilation (not used by breakup, but called unconditionally) | 2 |
| `cal_capacitance_vec` | class_Group.F90 | 7288–7407 | 119 | capacitance (unused by breakup) | 2 |
| `cal_coef_vapdep2_vec` | class_Group.F90 | 9348–9506 | 158 | vapor-dep coeffs (unused by breakup) | 2 |
| `cal_surface_temp2_vec` | class_Group.F90 | 4370–4571 | 201 | surface temp (unused by breakup) | 2 |
| `cal_Ice_Shape_v4`, `diag_habit_v4`, `ini_Ice_Shape_v4` | class_Ice_Shape.F90 | — | 111+133+… | ice-only; entered under `icond1(i,n)==0` mask | 2 |
| `cal_cs_spheroid3_vec`, `cal_bulk_density3_vec`, `diag_sh_type_v6` | mod_amps_core.F90 / class_Mass_Bin.F90 | — | 320+459+74 | ice geometry/density (unused for liquid) | 2 |
| **Cluster-2 subtotal (liquid-relevant path)** | | | **~2100+** | Group build + `diag_pq` diagnostic stack | |
| `print_bufd` | mod_amps_lib.F90 | 2019+ | ~40 | debug-file writer only; **not needed** | — |

## A.3 What inputs does the table actually need?

The critical finding: **`cal_breakfragment` runs on a single, hardcoded, fixed air state** — it is not a runtime routine. Inside the body:

```
T = 278.6795 K ;  PT = 850 hPa ;  W = 0 ;  phase = 1 (water) ;  RH = 100%
```

These are literals, not advected fields. `rv/den` are the length-1 dummies at that fixed state. The only "external" inputs are:

1. **Bin grid parameters** (`NRBIN`, `srat_r`, `sadd_r`, `minmass_r`, `binbr` from `com_amps`) — deterministic; define `liquid%binb` and `liquid%MS%len`.
2. **Saturation-vapor LUTs** `estbar(150)`, `esitbar(111)` — themselves setup-time LUTs, consumed only to build the single `steady` thermo state (which yields one scalar `sig_wa`, the surface tension of water at 278.68 K).
3. **Low-List fit coefficients** — all hardcoded constants inside `cal_breakup_dis_LL` (`app`, `bpp`, `CKE0`, `W0`, and the dozens of Low-List polynomial coefficients).
4. **The normal-CDF LUT/function** `getznorm2`→`cdfnor` — self-contained special-function code, no state.

What `cal_breakup_dis_LL` genuinely reads from the `Group`: `g_1%binb(:)` (bin mass boundaries) and `g_1%N_BIN`. What `cal_Coalescence_Efficiency`'s **token==1** (liquid-liquid) branch reads: `%len` (diameter, geometric from bin mass), `%vtm` (drop terminal velocity), and `th_var%sig_wa` (a function of T only). Nothing depends on runtime `Q`, transport, or timestep.

**Conclusion on infrastructure need:** The *breakup physics itself* (Cluster 1, ~1470 lines incl. CDFLIB) needs only bin boundaries + Low-List constants + the normal CDF — genuinely standalone. BUT the current code obtains `%len` and `%vtm` by calling the full `diag_pq` (Cluster 2, ~2100+ lines), which unconditionally invokes terminal-velocity, ventilation, capacitance, vapor-deposition, surface-temperature and (conditionally) ice-shape/bulk-density vec-routines. So *as written*, generating the table pulls in the full property-diagnosis infrastructure — even though the breakup math only needs `%len` (trivially geometric) and liquid-drop `%vtm` (one fall-speed formula). The ventilation/capacitance/vapor-dep/surface-temp/ice calls are dead weight for the table but are called anyway because `diag_pq` is monolithic.

## A.4 Decision: DUMP `bu_fd` / `bu_tmass`, do not port the generator in M2

`bu_fd` and `bu_tmass` are **pure setup-time constants** for a fixed bin configuration: the air state is hardcoded, the bins are deterministic, the coefficients are literals. They do not change at runtime and do not depend on the simulation.

**Recommended for M2: dump `bu_fd`, `bu_tmass`, and the four index scalars (`jmin_bk`, `imin_bk`, `imax_bk`, `jmax_bk`) plus the bin boundaries `binb` from one reference AMPS run, ship them as a static LUT (the same pattern the other LUT generators already use — cf. `amps/core/lookup_tables.py` and `codegen/convert_luts.py`), and load at setup.**

Evidence it is cheaper and safer:

- **Cheaper.** Porting the generator faithfully means reproducing Cluster 1 (~1470 lines, incl. the CDFLIB `cdfnor`/`cumnor`/`dinvnr` special-function chain) **and** the liquid path of `diag_pq` (mean-mass + density + `cal_terminal_vel_vec` at ~300 lines, since `%vtm` feeds `CKE`), realistically **>1800 lines** of new numeric code. Dumping is a handful of Fortran `write` statements against arrays that already exist in `com_amps`, plus a NumPy loader.
- **Safer.** The generator is FP-fragile: `cal_sig_sf`/`cal_Hmusig` iterate a **Brent root-finder** (`zbrent`) whose objective mixes `getznorm2` CDF evaluations with `min(0.99999…, …)` clamps and `-999.9` non-convergence fallbacks; `cal_breakup_dis_LL` sums lognormal/normal moment integrals with `dexp(-x²/2)` cancellation. Reproducing this bit-comparably in Python/GT4Py is exactly the catastrophic-cancellation / summation-order class of divergence that `spike_b_collection_codegen.py` already documented (its nbins=40 case needed rtol loosened to 1e-6 purely from reordered summation). A dumped table is bit-reproducible by construction and sidesteps all of it.
- **Setup-time constant.** There is zero runtime benefit to on-the-fly generation — the values are identical every run for a given `NRBIN`/bin grid.

**The only condition that flips the decision to "port":** if M2 must let the rain bin grid (`NRBIN`, `srat_r`, `minmass_r`, `binbr`) be reconfigurable at setup *and* the table must track arbitrary new grids without a reference-run regeneration step. If the target bin configuration is fixed (as in the cloudlab reference config all M0 dumps come from), dumping is strictly better. Recommend: **dump for M2; revisit porting only in a later milestone if grid-reconfigurability becomes a requirement**, and if so, port only Cluster 1 + a direct liquid-drop `%len`/`%vtm` computation — never the whole `diag_pq`.

---

# PART B — icon4py physics-scheme composition template for M2 `warm_phase.py`

Source: muphys `graupel.py` (the runnable-granule analog) + `run_graupel_only.py` (host loop) + the proven amps `spikes/` idioms.

## B.1 muphys `graupel` program-composition pattern

**Layered field_operator composition, one `@gtx.program` at the top.** The granule is built bottom-up (`graupel.py`):

1. **Leaf field_operators** (pure, per-cell transition math): `symmetric`, `cond_symmetric`, `sink_saturation`, and the big `_q_t_update` (the process "matrix": all vapor/cloud/rain/snow/ice/graupel transitions, sink-saturation limiting, then water-content + temperature update). These call the imported `core.transitions`/`core.properties`/`core.thermo` scalar helpers.
2. **A scan_operator + its wrappers** for sedimentation (see B.2).
3. **A mid-level field_operator** `_precipitation_effects` that computes the level-shifted `t_kp1 = concat_where(KDim < last_lev, t(KDim+1), t)` and calls the scan.
4. **The composing field_operator** `graupel(...)`: computes the qmin masks, runs `q, t = where(mask, _q_t_update(...), (q, te))`, then `_precipitation_effects(...)`, and returns the 8-tuple `(t, Q(...), pflx, pr, ps, pi, pg, pre)`.
5. **The `@gtx.program`** `graupel_run(...)`: the single entry point. It calls `graupel(...)` **once** with `out=(...)` and an explicit per-output `domain=(...)` tuple. Note the domain trick: prognostic fields use the full `(vertical_start, vertical_end)` K-range, whereas surface-flux diagnostics (`pr/ps/pi/pg/pre`) are written only on the bottom level `(vertical_end-1, vertical_end)`.

Chaining is **static composition of field_operators inside one program**, not a Python-level pipeline. State between stages is threaded as `NamedTuple`s (`Q`, `PrecipStateQx`, `TempState`, `IntegrationState`).

## B.2 The sedimentation scan_operator (`_precip_and_t`)

The vertical-implicit sedimentation is a `@gtx.scan_operator` — the pattern M2's bin-sedimentation must follow (proven wide by `spike_c_wide_scan.py`):

- **Decorator / init state:** `@gtx.scan_operator(axis=dims.KDim, forward=True, init=IntegrationState(...))`. The `init` is a fully-populated `NamedTuple` carry: four `PrecipStateQx` (r/s/i/g, each `x,p,vc,activated`), a `TempState`, plus `rho` and `pflx_tot`. `activated=False` flags gate the first-touch level.
- **Per-level body:** `def _precip_and_t(previous_level: IntegrationState, t, t_kp1, rho, q, mask_r,…, dt, dz)`. It computes `zeta = dt/(2*dz)`, per-class fall-speed scale factors, an `any_mask | previous_level_activated` gate, then calls the per-class helper `precip_qx_level_update(previous_level.r, previous_level.rho, prefactor, exponent, offset, zeta, vc, q, rho, mask)` for each hydrometeor. Each update reads `previous_level_q.p`/`.x`/`.vc` (the carry) → genuine downward recurrence. Returns a new `IntegrationState`.
- **Scan-output convention** (from `spike_c_wide_scan.py`, measured not assumed): gt4py emits the *post-update* carry at each level — the returned `NamedTuple` is simultaneously the next level's input and this level's output field value.
- **Wrapping:** the scan is called inside `_precipitation_effects` (a field_operator), whose result tuple is unpacked into the exposed CellK fields. This "scan inside a field_operator that unpacks the carry to output fields" is exactly `spike_c`'s `_sed_{n}` wrapper shape.

## B.3 `setup_program` constant-arg binding (`run_graupel_only.py::setup_graupel`)

```python
graupel_run_program = model_options.setup_program(
    backend=backend,
    program=graupel.graupel_run,
    constant_args={"dt": ta.wpfloat(dt), "qnc": ta.wpfloat(qnc),
                   "enable_masking": enable_masking},
    horizontal_sizes={"horizontal_start": 0, "horizontal_end": ncells},
    vertical_sizes={"vertical_start": 0, "vertical_end": nlev},
    offset_provider={},
)
gtx.wait_for_compilation()
```

Key points M2 must copy:
- **`constant_args`** bakes scalars (`dt`, `qnc`, config flags) into the compiled program so the per-step call site passes only fields.
- The whole setup is wrapped in **`utils.recursion_limit(10**4)`** — mandatory: gt4py's ITIR transform pipeline recurses per-node and overflows Python's default limit on large unrolled trees. This is the same `muphys_driver_utils.recursion_limit` the amps spikes (`spike_b`, `spike_c`) already reuse at `10**5`.
- `setup_program` returns a **bound callable**; horizontal/vertical extents are fixed at setup, not per call.

## B.4 The host time-loop driver (`run_graupel_only.py::main`)

Box/column driver shape:
1. Resolve backend + allocator (`model_backends.BACKENDS[...]`, `get_allocator`), pick dtype from `ta.precision`.
2. Load input into an in/out field bundle (`common.GraupelInput.load(...)`, `GraupelOutput.allocate(...)`). Note the `use_inout_buffers` idiom: same buffers passed as input and output (safe because no offset reads).
3. Call `setup_graupel(...)` once (compile).
4. **Time loop:** `for _x in range(itime+1): graupel_run_program(dz=…, te=…, p=…, rho=…, q_in=…, t_out=…, q_out=…, pflx=…, pr=…, …)` — only fields at the call site (scalars already bound). Timing starts on iteration 1 (skip cold-compile); `device_utils.sync(allocator)` brackets the timed region.
5. Write output.

The amps `driver/box.py` is the single-column (`ncells=1`) analog and already stubs this: `run_box(case)` is the declared **M2 wiring point** — its docstring states it must assemble "the ported DSL/numpy microphysics process kernels … into a per-timestep driver loop over `case.n_steps`". `case_from_micro_record` already builds a runnable `BoxCase` from a dumped `MicroRecord` (phase=1 "pre") for replay/diff against the phase=2 "post" record.

## B.5 Proven amps codegen + scan idioms M2 must reuse

From `spikes/` docstrings and `codegen/`:

- **SPLIT operators are BINDING** (`codegen/generate.py` + `templates.py` docstrings, M0-gate amendment): every codegen builder MUST emit **one `@gtx.field_operator` per bin-group of ≤8 consecutive bins**, never one monolithic 40-bin operator. Evidence: `spike_b` measured a single unrolled 40-bin collection kernel at **~2579 s (~43 min) gtfn_cpu compile**. Use the shared `templates.chunk_bins(nbins, chunk_size=8)` helper; the first real (non-demo) consumer named is "M2's collection kernel," chunked **by destination bin**.
- **Wide scan is validated** (`spike_c`): an 80-field NamedTuple carry (2×40 bins) compiles gtfn_cpu in **186.9 s** and runs 115.9 ms/call at NCELLS=4096/NLEV=61 — a GO. This de-risks a 40-bin sedimentation scan. Includes an executable **perturbation check** (perturb level-0 q, assert it propagates to level-1 output) proving the carry is a genuine sequential recurrence, not degenerate level-parallel — M2's sed scan should carry an equivalent test.
- **Codegen source hygiene:** every generated float literal must be wrapped `float(...)` (numpy≥2.0 reprs `np.float64` as `"np.float64(0.1)"`, which emits an undefined `np` symbol into generated source). Load generated operators via `common.load_generated_operator` (`importlib.util.spec_from_file_location`, fresh module each time; needs a real on-disk file because gt4py's decorator uses `inspect.getsource`). Generated files are committed under `core/generated/` with a drift-guard (`check_regenerated`).
- **Cache validity for any timing:** clear `<worktree>/.gt4py_cache` before cold-compile measurements (`GT4PY_BUILD_CACHE_LIFETIME=PERSISTENT` in this env).

## B.6 Concrete composition template for `warm_phase.py`

```
warm_phase.py
├── NamedTuple carries
│     Qwarm(...)                     # warm-phase specific masses / bin moments
│     WarmScanState(...)             # sedimentation carry (per-bin, ≤80 fields; spike_c-proven)
│
├── leaf @gtx.field_operator's       # pure per-cell warm-rain process math
│     _collision_coalescence(...)    # CODEGEN'd, SPLIT ≤8 bins/operator, chunked by dest bin
│                                     #   (uses dumped bu_fd/bu_tmass LUT from Part A)
│     _condensation_evaporation(...)
│     _q_update(...)                 # analog of muphys _q_t_update: sink-saturation + q/t update
│
├── @gtx.scan_operator(axis=KDim, forward=True, init=WarmScanState(...))
│     _warm_sedimentation(previous_level, <per-bin q args>, rho, zeta) -> WarmScanState
│         # per-level implicit fall, muphys precip_qx_level_update shape, spike_c wide carry
│
├── mid @gtx.field_operator
│     _sedimentation_effects(last_level, kmin_*, q, t, rho, dz, dt)
│         # computes t_kp1 via concat_where(KDim < last_lev, t(KDim+1), t); calls the scan; unpacks
│
├── composing @gtx.field_operator
│     warm_phase(last_level, dz, te, p, rho, q, dt, <config scalars>, enable_masking) -> tuple
│         # qmin masks; q,t = where(mask, _q_update(...), (q, te)); then _sedimentation_effects
│
└── @gtx.program(grid_type=UNSTRUCTURED)
      warm_phase_run(<in fields>, dt, <config scalars>, <out fields>,
                     horizontal_start/end, vertical_start/end, enable_masking) -> None
          # single call to warm_phase(..., out=(...), domain=(... per-output ranges ...))
```

Host side (mirror `setup_graupel` + `run_box`):

```python
def setup_warm_phase(*, dt, backend, horizontal_start, horizontal_end,
                     vertical_start, vertical_end, **config):
    with utils.recursion_limit(10**5):                     # MANDATORY (ITIR recursion)
        prog = model_options.setup_program(
            backend=backend, program=warm_phase.warm_phase_run,
            constant_args={"dt": ta.wpfloat(dt), **config},  # bake scalars
            horizontal_sizes={...}, vertical_sizes={...}, offset_provider={})
    gtx.wait_for_compilation()
    return prog

# time loop (fill in driver/box.py::run_box):
prog = setup_warm_phase(dt=case.dt, backend=backend, ...)
for _ in range(case.n_steps):
    prog(dz=..., te=..., p=..., rho=..., q_in=..., q_out=..., t_out=..., <flux outs>)
```

Non-negotiables carried from the spikes: (1) collision-coalescence kernel emitted as **SPLIT ≤8-bin field_operators** via `codegen`, never monolithic; (2) sedimentation as a **single wide `scan_operator`** with a NamedTuple carry + a perturbation test; (3) `float(...)`-wrap every codegen literal; (4) `recursion_limit` around setup; (5) `warm_phase_run` is the sole `@gtx.program` entry point, with per-output `domain` ranges (full K for prognostics, bottom level for surface fluxes).