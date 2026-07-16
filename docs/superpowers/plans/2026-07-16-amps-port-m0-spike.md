# AMPS→icon4py Port — M0 (Spike + Scaffolding + Instrumentation) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver M0 of the AMPS port (spec: `docs/superpowers/specs/2026-07-16-amps-icon4py-port-design.md`): SCALE instrumentation patch + run instructions for reference data, the new `amps` icon4py package skeleton registered in the workspace, and five GT4Py feasibility spikes with a written go/no-go report.

**Architecture:** Two repos. `scale_amps` (branch `cloudlab_port`, already checked out at `/Users/jcanton/projects/scale_amps`) gets a binary-stream dump instrumentation inside `scale_atmos_phy_mp_amps.F90` plus a Python reader. `icon4py` gets a new worktree `.worktrees/amps_microphysics` (branch `amps_microphysics` off `main`) with the scaffolded `amps` package and spike scripts under `spikes/`.

**Tech Stack:** Fortran (SCALE), Python 3.12, uv workspace, gt4py==1.1.11, numpy, pytest.

## Global Constraints

- icon4py: run `uv sync --extra all` and all `uv run` commands **from the worktree root** (`/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics`), never from a subpackage.
- icon4py: gt4py pinned `==1.1.11`; Python 3.12; ruff line length 100; every new file starts with the 7-line license header (content of `HEADER.txt` as `#` comments — see Task 6 for the exact block).
- icon4py: package name `icon4py-atmosphere-amps`, version `0.2.0`, import path `icon4py.model.atmosphere.subgrid_scale_physics.amps`, depends only on `icon4py.model.common` (tach rule).
- GT4Py DSL limits (from recon, drive the spikes): no `for`/`while`; no bitwise ops or shifts on integers (only `+ - * / // %` — `%` requires integral right operand); `as_offset` is gather-only; `if` requires scalar condition.
- SCALE: dumps are float64 regardless of build (`real(a,8)` on write); all new namelist vars go into `PARAM_ATMOS_PHY_MP_AMPS_bin`; nothing inside the OpenMP parallel region may lexically reference a new module variable unless it is added to the `!$omp shared(...)` list (`default(none)` is in force) — the patch avoids this by referencing module state only inside called module procedures.
- SCALE cannot be built locally: Fortran verification = careful diff self-review; compile happens on the user's cluster.
- Commits: conventional style. scale_amps commits go to `cloudlab_port` and are pushed to `origin`. icon4py commits go to `amps_microphysics`, local only (do NOT push; remote choice is the user's).
- Reference for muphys conventions being mirrored: `model/atmosphere/subgrid_scale_physics/muphys/` on icon4py `main`.

---

### Task 1: SCALE dump infrastructure (namelist, files, writer helpers)

**Files:**
- Modify: `scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90` (module vars ~line 146, namelist ~line 363, open-hook ~line 1426, flush ~line 3719, new private procedures before `end module`)

**Interfaces:**
- Produces (used by Tasks 2–3): module vars `l_amps_dump`, `amps_dump_step_stride`, `amps_dump_is/ie/js/je`, `amps_dump_dir`; procedures `AMPS_DUMP_open(nsect)`, `AMPS_DUMP_flush()`, `logical function AMPS_DUMP_active(i,j)`, writers `AMPS_DUMP_w_i0/w_i1/w_r0/w_r1/w_r2/w_r3/w_r4(fid, a)`; per-thread unit array `amps_dump_fid(:)`.
- Record byte conventions (used by Task 4 reader): all ints written as int32 (`int(x,4)`), all reals as float64 (`real(x,8)`), stream access, native (little) endian; every array is prefixed by its int32 dimension sizes (1 per rank); scalars via `w_i0`/`w_r0` have no prefix.

- [ ] **Step 1: Add module variables**

In `scale_atmos_phy_mp_amps.F90`, directly after line 146 (`logical  :: fix_aerosol_type(4) = ...`), insert:

```fortran
  !--- reference-data dump instrumentation (AMPS->icon4py port, M0)
  logical            :: l_amps_dump           = .false. ! master switch for binary dumps
  character(len=256) :: amps_dump_dir         = "."     ! output directory (must exist)
  integer            :: amps_dump_step_stride = 300     ! dump every Nth MP step (TIME_AMPS)
  integer            :: amps_dump_is          = 0       ! local i range of dumped columns
  integer            :: amps_dump_ie          = -1      ! (ie<is disables; halo-inclusive local indices)
  integer            :: amps_dump_js          = 0
  integer            :: amps_dump_je          = -1
  integer, allocatable :: amps_dump_fid(:)              ! one unit per OpenMP thread (isect)
  logical            :: amps_dump_opened      = .false.
```

- [ ] **Step 2: Extend the namelist**

At the namelist (line 346–363), change the last two entries and append (keep trailing comment style):

```fortran
       amps_debug,         & ! debugging on or off
       amps_ignore,        & ! ignore amps microphysics or not
       l_amps_dump,        & ! write binary reference dumps around ifc_cloud_micro/sedimentation
       amps_dump_dir,      & ! directory for dump files
       amps_dump_step_stride, & ! dump every Nth MP step
       amps_dump_is,       & ! local i-range start of dumped columns
       amps_dump_ie,       & ! local i-range end   (ie<is disables)
       amps_dump_js,       & ! local j-range start
       amps_dump_je          ! local j-range end
```

(Note: `amps_ignore` was previously the last item without a comma — add the comma.)

- [ ] **Step 3: Add the dump procedures**

Immediately before `end module scale_atmos_phy_mp_amps` (end of file), insert:

```fortran
  !-----------------------------------------------------------------------------
  ! Reference-data dump instrumentation (AMPS->icon4py port).
  ! Format: per-rank-per-thread unformatted stream files. Ints are int32,
  ! reals are float64. Arrays are prefixed by one int32 size per rank.
  !-----------------------------------------------------------------------------
  subroutine AMPS_DUMP_open(nsect)
    use scale_prc, only: PRC_myrank, PRC_abort
    use scale_io, only: IO_get_available_fid
    integer, intent(in) :: nsect
    character(len=512) :: fname
    integer :: is_, ierr
    allocate(amps_dump_fid(nsect))
    do is_ = 1, nsect
       write(fname,'(A,A,I6.6,A,I3.3,A)') trim(amps_dump_dir), '/amps_dump_r', PRC_myrank, '_t', is_, '.bin'
       amps_dump_fid(is_) = IO_get_available_fid()
       open( unit   = amps_dump_fid(is_), &
             file   = trim(fname),        &
             form   = 'unformatted',      &
             access = 'stream',           &
             status = 'replace',          &
             iostat = ierr )
       if ( ierr /= 0 ) then
          LOG_ERROR("AMPS_DUMP_open",*) "cannot open dump file: ", trim(fname)
          call PRC_abort
       end if
    end do
    amps_dump_opened = .true.
  end subroutine AMPS_DUMP_open

  subroutine AMPS_DUMP_flush()
    integer :: is_
    if ( .not. amps_dump_opened ) return
    do is_ = 1, size(amps_dump_fid)
       flush(amps_dump_fid(is_))
    end do
  end subroutine AMPS_DUMP_flush

  logical function AMPS_DUMP_active(i, j)
    integer, intent(in) :: i, j
    AMPS_DUMP_active = l_amps_dump                                   &
                 .and. mod(TIME_AMPS, amps_dump_step_stride) == 0    &
                 .and. i >= amps_dump_is .and. i <= amps_dump_ie     &
                 .and. j >= amps_dump_js .and. j <= amps_dump_je
  end function AMPS_DUMP_active

  subroutine AMPS_DUMP_w_i0(fid, a)
    integer, intent(in) :: fid, a
    write(fid) int(a,4)
  end subroutine AMPS_DUMP_w_i0

  subroutine AMPS_DUMP_w_i1(fid, a)
    integer, intent(in) :: fid
    integer, intent(in) :: a(:)
    write(fid) int(size(a,1),4)
    write(fid) int(a,4)
  end subroutine AMPS_DUMP_w_i1

  subroutine AMPS_DUMP_w_r0(fid, a)
    integer, intent(in) :: fid
    real(RP), intent(in) :: a
    write(fid) real(a,8)
  end subroutine AMPS_DUMP_w_r0

  subroutine AMPS_DUMP_w_r1(fid, a)
    integer, intent(in) :: fid
    real(RP), intent(in) :: a(:)
    write(fid) int(size(a,1),4)
    write(fid) real(a,8)
  end subroutine AMPS_DUMP_w_r1

  subroutine AMPS_DUMP_w_r2(fid, a)
    integer, intent(in) :: fid
    real(RP), intent(in) :: a(:,:)
    write(fid) int(size(a,1),4), int(size(a,2),4)
    write(fid) real(a,8)
  end subroutine AMPS_DUMP_w_r2

  subroutine AMPS_DUMP_w_r3(fid, a)
    integer, intent(in) :: fid
    real(RP), intent(in) :: a(:,:,:)
    write(fid) int(size(a,1),4), int(size(a,2),4), int(size(a,3),4)
    write(fid) real(a,8)
  end subroutine AMPS_DUMP_w_r3

  subroutine AMPS_DUMP_w_r4(fid, a)
    integer, intent(in) :: fid
    real(RP), intent(in) :: a(:,:,:,:)
    write(fid) int(size(a,1),4), int(size(a,2),4), int(size(a,3),4), int(size(a,4),4)
    write(fid) real(a,8)
  end subroutine AMPS_DUMP_w_r4
```

Note: `LOG_ERROR` requires the module's existing `use scale_io` macros — this file already uses them (`LOG_INFO` at line 370). `TIME_AMPS` is the module step counter (line 130). `RP` comes from `use scale_precision` (line 19).

- [ ] **Step 4: Hook open + flush into the tendency routine**

After line 1424 (`!$ nsect = omp_get_max_threads()`), insert:

```fortran
    if ( l_amps_dump .and. .not. amps_dump_opened ) call AMPS_DUMP_open(nsect)
```

Immediately before line 3720 (`TIME_AMPS = TIME_AMPS + 1`), insert:

```fortran
    if ( l_amps_dump ) call AMPS_DUMP_flush()
```

(Line numbers shift as you edit — locate by content, not number.)

- [ ] **Step 5: Self-review the diff**

Run: `git -C /Users/jcanton/projects/scale_amps diff` and check: (a) no new variable is lexically referenced inside the `!$omp parallel do` region (the hooks in Step 4 are outside it; `AMPS_DUMP_active`/writer calls come in Tasks 2–3 and only reference module state inside procedures); (b) namelist syntax valid (commas, continuation `&`); (c) writer helpers compile-plausible (intent, kind conversions).

- [ ] **Step 6: Commit**

```bash
cd /Users/jcanton/projects/scale_amps
git add scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90
git commit -m "feat(amps-dump): dump infrastructure for port reference data"
```

---

### Task 2: SCALE per-call micro dumps around `ifc_cloud_micro`

**Files:**
- Modify: `scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90` (new procedure before `end module`; call-site edits around line 2071)

**Interfaces:**
- Consumes: Task 1 writers, `AMPS_DUMP_active`, `amps_dump_fid`.
- Produces: **micro record layout v1** (Task 4 reader must match exactly):
  header int32s in order: `magic=1095586131, version=1, phase (1=pre|2=post), TIME_AMPS, i, j, isect, nmic, npr, nbr, ncr, npi, nbi, nci, npa, nba, nca, mxnbin, istrt, jseed, ifrst, isect_seed, nextn` (23 × int32), then float64 `dt`;
  then arrays in order: `kmicvm` (i1), `qcvm, v3v, qvvm, moist_denvm, ptotvm, tvm, wbvm, trpv_thil, trpv_qtp` (each r1, length nmic), `qrpvm` (r4: npr,nbr,ncr,nmic), `qipvm` (r4), `qapvm` (r4);
  phase 2 only, appended: `dmtendlm` (r3: 10,2,nmic), `dcontendlm` (r3), `dbintendlm` (r4: 7,2,mxnbin,nmic).

- [ ] **Step 1: Add the micro-dump procedure**

Before `end module` (next to the Task 1 procedures), insert:

```fortran
  subroutine AMPS_DUMP_micro(phase, isect, i, j, nmic, istrt, dt,           &
                             jseed_t, ifrst_t, isect_seed_t, nextn_t,       &
                             kmicvm, qcvm, v3v, qvvm, moist_denvm, ptotvm,  &
                             tvm, wbvm, trpvm, qrpvm, qipvm, qapvm,         &
                             dmtendlm, dcontendlm, dbintendlm)
    integer,  intent(in) :: phase, isect, i, j, nmic, istrt
    real(RP), intent(in) :: dt
    integer,  intent(in) :: jseed_t, ifrst_t, isect_seed_t, nextn_t
    integer,  intent(in) :: kmicvm(:)
    real(RP), intent(in) :: qcvm(:), v3v(:), qvvm(:), moist_denvm(:), ptotvm(:), tvm(:), wbvm(:)
    real(RP), intent(in) :: trpvm(:,:)
    real(RP), intent(in) :: qrpvm(:,:,:,:), qipvm(:,:,:,:), qapvm(:,:,:,:)
    real(RP), intent(in) :: dmtendlm(:,:,:), dcontendlm(:,:,:), dbintendlm(:,:,:,:)
    integer :: fid
    fid = amps_dump_fid(isect)
    write(fid) int(1095586131,4), int(1,4)
    write(fid) int(phase,4), int(TIME_AMPS,4), int(i,4), int(j,4), int(isect,4), int(nmic,4)
    write(fid) int(npr,4), int(nbr,4), int(ncr,4), int(npi,4), int(nbi,4), int(nci,4), &
               int(npa,4), int(nba,4), int(nca,4), int(mxnbin,4)
    write(fid) int(istrt,4), int(jseed_t,4), int(ifrst_t,4), int(isect_seed_t,4), int(nextn_t,4)
    call AMPS_DUMP_w_r0(fid, dt)
    call AMPS_DUMP_w_i1(fid, kmicvm(1:nmic))
    call AMPS_DUMP_w_r1(fid, qcvm(1:nmic))
    call AMPS_DUMP_w_r1(fid, v3v(1:nmic))
    call AMPS_DUMP_w_r1(fid, qvvm(1:nmic))
    call AMPS_DUMP_w_r1(fid, moist_denvm(1:nmic))
    call AMPS_DUMP_w_r1(fid, ptotvm(1:nmic))
    call AMPS_DUMP_w_r1(fid, tvm(1:nmic))
    call AMPS_DUMP_w_r1(fid, wbvm(1:nmic))
    call AMPS_DUMP_w_r1(fid, trpvm(1:nmic,1))   ! thil
    call AMPS_DUMP_w_r1(fid, trpvm(1:nmic,2))   ! qtp
    call AMPS_DUMP_w_r4(fid, qrpvm(:,:,:,1:nmic))
    call AMPS_DUMP_w_r4(fid, qipvm(:,:,:,1:nmic))
    call AMPS_DUMP_w_r4(fid, qapvm(:,:,:,1:nmic))
    if ( phase == 2 ) then
       call AMPS_DUMP_w_r3(fid, dmtendlm(:,:,1:nmic))
       call AMPS_DUMP_w_r3(fid, dcontendlm(:,:,1:nmic))
       call AMPS_DUMP_w_r4(fid, dbintendlm(:,:,:,1:nmic))
    end if
  end subroutine AMPS_DUMP_micro
```

`npr..nca` are module publics (lines 291–299); `mxnbin` comes from `use maxdims` (line 23) — all visible to module procedures without OMP clauses.

- [ ] **Step 2: Wrap call site #1**

At the call site (line 2071 context: inside `if(level.ge.5) then` ... `call PROF_rapstart("amps_micro",3)`), insert immediately BEFORE `call ifc_cloud_micro(`:

```fortran
             if ( AMPS_DUMP_active(i,j) ) then
                call AMPS_DUMP_micro(1, isect, i, j, nmic, istrt, dt,             &
                                     jseed(isect), ifrst(isect), isect_seed(isect), nextn(isect), &
                                     kmicvm, qcvm, v3v, qvvm, moist_denvm, ptotvm, &
                                     tvm, wbvm, trpvm, qrpvm, qipvm, qapvm,        &
                                     dmtendlm, dcontendlm, dbintendlm)
             end if
```

and immediately AFTER `call PROF_rapend("amps_micro",3)`:

```fortran
             if ( AMPS_DUMP_active(i,j) ) then
                call AMPS_DUMP_micro(2, isect, i, j, nmic, istrt, dt,             &
                                     jseed(isect), ifrst(isect), isect_seed(isect), nextn(isect), &
                                     kmicvm, qcvm, v3v, qvvm, moist_denvm, ptotvm, &
                                     tvm, wbvm, trpvm, qrpvm, qipvm, qapvm,        &
                                     dmtendlm, dcontendlm, dbintendlm)
             end if
```

All actual arguments (`qcvm`, `qrpvm`, `jseed`, `ifrst`, `isect_seed`, `nextn`, `dmtendlm`, ...) are already in the OMP private/shared lists — no `!$omp` edits needed (module state is referenced only inside `AMPS_DUMP_*` procedures). Note the pre-call dump captures `istrt` before `micro_io_strt(isect)=.false.` resets it (that happens after the call, line 2089). The dead second call site (line 3524, inside `if (.false.)`) is NOT instrumented.

- [ ] **Step 3: Self-review the diff**

`git diff` — check argument order matches the procedure signature exactly, and both wrappers sit inside the `if(level.ge.5)` block.

- [ ] **Step 4: Commit**

```bash
cd /Users/jcanton/projects/scale_amps
git add scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90
git commit -m "feat(amps-dump): per-call in/out dumps around ifc_cloud_micro"
```

---

### Task 3: SCALE sedimentation dumps

**Files:**
- Modify: `scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90` (new procedure; call-site edits around lines 2364–2427)

**Interfaces:**
- Consumes: Task 1 writers.
- Produces: **sed record layout v1** (Task 4 reader must match):
  header int32s: `magic=1095586132, version=1, phase (3=pre|4=post), TIME_AMPS, i, j, isect, isn (0=liquid|1=ice), iadvv, np, nb, nc, k1, k2, k1m, k2m` (16 × int32), float64 `dt`;
  then arrays: `k1b` (i1: nb*nc flattened via reshape), `k2b` (i1), `qpv` (r4: np,nb,nc,nzh-length K), `q_this` (r1), `q_other` (r1), `qcv` (r1), `qtp` (r1), `moist_denv` (r1), `thetav` (r1), `qvv` (r1), `tv` (r1), `dens_col` (r1), `momz_col` (r1), `u_col` (r1), `v_col` (r1), `cz_col` (r1), `fz_col` (r1), `dzzmv` (r1), `dzvmv` (r1), `mmass` (r3: nb,nc,K), `den_t` (r1), `momz_t` (r1), `rhou_t` (r1), `rhov_t` (r1), `rhoe_t` (r1), `sflx` (r0).
  Same layout for phase 4 (post), arrays reflecting post-call state.

- [ ] **Step 1: Add the sed-dump procedure**

Before `end module`, insert:

```fortran
  subroutine AMPS_DUMP_sed(phase, isect, i, j, isn, iadvv_l, dt,             &
                           np, nb, nc, k1, k2, k1m, k2m, k1b, k2b,           &
                           qpv, q_this, q_other, qcv, qtp, moist_denv,       &
                           thetav, qvv, tv, dens_col, momz_col, u_col, v_col,&
                           cz_col, fz_col, dzzmv, dzvmv, mmass,              &
                           den_t, momz_t, rhou_t, rhov_t, rhoe_t, sflx)
    integer,  intent(in) :: phase, isect, i, j, isn, iadvv_l, np, nb, nc, k1, k2, k1m, k2m
    real(RP), intent(in) :: dt
    integer,  intent(in) :: k1b(:,:), k2b(:,:)
    real(RP), intent(in) :: qpv(:,:,:,:), mmass(:,:,:)
    real(RP), intent(in) :: q_this(:), q_other(:), qcv(:), qtp(:), moist_denv(:), thetav(:), qvv(:), tv(:)
    real(RP), intent(in) :: dens_col(:), momz_col(:), u_col(:), v_col(:), cz_col(:), fz_col(:)
    real(RP), intent(in) :: dzzmv(:), dzvmv(:)
    real(RP), intent(in) :: den_t(:), momz_t(:), rhou_t(:), rhov_t(:), rhoe_t(:)
    real(RP), intent(in) :: sflx
    integer :: fid
    fid = amps_dump_fid(isect)
    write(fid) int(1095586132,4), int(1,4)
    write(fid) int(phase,4), int(TIME_AMPS,4), int(i,4), int(j,4), int(isect,4), int(isn,4)
    write(fid) int(iadvv_l,4), int(np,4), int(nb,4), int(nc,4)
    write(fid) int(k1,4), int(k2,4), int(k1m,4), int(k2m,4)
    call AMPS_DUMP_w_r0(fid, dt)
    call AMPS_DUMP_w_i1(fid, reshape(k1b, (/size(k1b)/)))
    call AMPS_DUMP_w_i1(fid, reshape(k2b, (/size(k2b)/)))
    call AMPS_DUMP_w_r4(fid, qpv)
    call AMPS_DUMP_w_r1(fid, q_this)
    call AMPS_DUMP_w_r1(fid, q_other)
    call AMPS_DUMP_w_r1(fid, qcv)
    call AMPS_DUMP_w_r1(fid, qtp)
    call AMPS_DUMP_w_r1(fid, moist_denv)
    call AMPS_DUMP_w_r1(fid, thetav)
    call AMPS_DUMP_w_r1(fid, qvv)
    call AMPS_DUMP_w_r1(fid, tv)
    call AMPS_DUMP_w_r1(fid, dens_col)
    call AMPS_DUMP_w_r1(fid, momz_col)
    call AMPS_DUMP_w_r1(fid, u_col)
    call AMPS_DUMP_w_r1(fid, v_col)
    call AMPS_DUMP_w_r1(fid, cz_col)
    call AMPS_DUMP_w_r1(fid, fz_col)
    call AMPS_DUMP_w_r1(fid, dzzmv)
    call AMPS_DUMP_w_r1(fid, dzvmv)
    call AMPS_DUMP_w_r3(fid, mmass)
    call AMPS_DUMP_w_r1(fid, den_t)
    call AMPS_DUMP_w_r1(fid, momz_t)
    call AMPS_DUMP_w_r1(fid, rhou_t)
    call AMPS_DUMP_w_r1(fid, rhov_t)
    call AMPS_DUMP_w_r1(fid, rhoe_t)
    call AMPS_DUMP_w_r0(fid, sflx)
  end subroutine AMPS_DUMP_sed
```

- [ ] **Step 2: Wrap the two `sclsedprz_original` calls**

Liquid (context lines 2364–2391): immediately before `call sclsedprz_original(qrv,qiv,qcv,trpv(:,2) ...` insert:

```fortran
                if ( AMPS_DUMP_active(i,j) ) then
                   call AMPS_DUMP_sed(3, isect, i, j, 0, iadvv, dt,                     &
                                      npr, nbr, ncr, k1r, k2r, k1m, k2m, k1br, k2br,    &
                                      qrpv, qrv, qiv, qcv, trpv(:,2), moist_denv,       &
                                      thetav, qvv, tv, DENS(KS-1:KE,i,j), MOMZ(KS-1:KE,i,j), &
                                      U(KS-1:KE,i,j), V(KS-1:KE,i,j),                   &
                                      CZ(KS-1:KE,i,j), FZ(KS-1:KE,i,j), dzzmv, dzvmv, mmassrv, &
                                      den_t, MOMZ_t(KS-1:KE,i,j), RHOU_t(KS-1:KE,i,j),  &
                                      RHOV_t(KS-1:KE,i,j), RHOE_t(KS-1:KE,i,j), SFLX_rain(i,j))
                end if
```

and immediately after `SFLX_rain(i,j) = abs(SFLX_rain(i,j))` the same call with `phase=4` (first argument `4` instead of `3`, everything else identical).

Ice (context lines 2399–2426): same pattern around the second `sclsedprz_original` call, with arguments `4→isn=1`: `call AMPS_DUMP_sed(3, isect, i, j, 1, iadvv, dt, npi, nbi, nci, k1i, k2i, k1m, k2m, k1bi, k2bi, qipv, qiv, qrv, qcv, trpv(:,2), moist_denv, thetav, qvv, tv, DENS(...), MOMZ(...), U(...), V(...), CZ(...), FZ(...), dzzmv, dzvmv, mmassiv, den_t, MOMZ_t(...), RHOU_t(...), RHOV_t(...), RHOE_t(...), SFLX_snow(i,j))` before, and `phase=4` after `SFLX_snow(i,j) = abs(...)`.

**Rank check (required):** before finalizing, grep the declarations of `dzzmv`, `dzvmv`, `k1br/k2br/k1bi/k2bi`, `k1m/k2m` in this file (declared in the locals block ~lines 1272–1330 and nearby). Expected: `dzzmv`, `dzvmv` rank-1 columns; `k1br(nbr,ncr)` etc. rank-2 integer; `k1m,k2m` scalars. If any differs, adjust the corresponding writer call (`w_r1`↔`w_r2`, scalar↔array) AND record the change in the reader field table (Task 4) — layouts must stay in sync.

- [ ] **Step 3: Self-review diff + commit**

```bash
cd /Users/jcanton/projects/scale_amps
git diff   # verify argument orders, phase codes, isn codes
git add scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90
git commit -m "feat(amps-dump): sedimentation in/out dumps"
```

---

### Task 4: Python dump reader + self-test + run instructions

**Files:**
- Create: `scripts/amps_dump_reader.py`
- Create: `scripts/test_amps_dump_reader.py`
- Create: `docs/superpowers/specs/2026-07-16-ref-data-run-instructions.md`

**Interfaces:**
- Consumes: record layouts v1 from Tasks 2–3 (field order fixed there).
- Produces: `read_dump_file(path) -> list[dict]` (each dict: `kind` ("micro"|"sed"), header fields by name, arrays by name as numpy arrays with Fortran dim order restored); CLI `python scripts/amps_dump_reader.py DIR -o OUT.npz` writing keys `micro_t{T}_i{I}_j{J}_{pre|post}_{field}` and `sed_t{T}_i{I}_j{J}_s{ISN}_{pre|post}_{field}`.

- [ ] **Step 1: Write the failing self-test**

`scripts/test_amps_dump_reader.py`:

```python
#!/usr/bin/env python3
"""Round-trip test for amps_dump_reader: build synthetic records byte-identically
to the Fortran writers, then parse and compare. Run: python scripts/test_amps_dump_reader.py"""
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from amps_dump_reader import MAGIC_MICRO, MAGIC_SED, read_dump_file


def w_i0(b, v):
    b += struct.pack("<i", v)
    return b


def w_i1(b, a):
    a = np.asarray(a, dtype="<i4")
    return b + struct.pack("<i", a.size) + a.tobytes()


def w_r0(b, v):
    return b + struct.pack("<d", v)


def w_rn(b, a):
    """Array with per-rank int32 size prefixes, Fortran order."""
    a = np.asarray(a, dtype="<f8")
    for s in a.shape:
        b += struct.pack("<i", s)
    return b + a.tobytes(order="F")


def build_micro(phase, nmic=3, npr=6, nbr=4, ncr=1, npi=18, nbi=2, nci=1,
                npa=5, nba=1, nca=4, mxnbin=4):
    rng = np.random.default_rng(phase)
    b = b""
    for v in (MAGIC_MICRO, 1, phase, 42, 3, 4, 1, nmic,
              npr, nbr, ncr, npi, nbi, nci, npa, nba, nca, mxnbin,
              1, 12345, 0, 7, 99):
        b = w_i0(b, v)
    b = w_r0(b, 1.0)
    b = w_i1(b, np.arange(2, 2 + nmic))
    fields = {}
    for name in ("qcvm", "v3v", "qvvm", "moist_denvm", "ptotvm", "tvm", "wbvm",
                 "trpv_thil", "trpv_qtp"):
        fields[name] = rng.uniform(size=nmic)
        b = w_rn(b, fields[name])
    fields["qrpvm"] = rng.uniform(size=(npr, nbr, ncr, nmic))
    fields["qipvm"] = rng.uniform(size=(npi, nbi, nci, nmic))
    fields["qapvm"] = rng.uniform(size=(npa, nba, nca, nmic))
    for name in ("qrpvm", "qipvm", "qapvm"):
        b = w_rn(b, fields[name])
    if phase == 2:
        fields["dmtendlm"] = rng.uniform(size=(10, 2, nmic))
        fields["dcontendlm"] = rng.uniform(size=(10, 2, nmic))
        fields["dbintendlm"] = rng.uniform(size=(7, 2, mxnbin, nmic))
        for name in ("dmtendlm", "dcontendlm", "dbintendlm"):
            b = w_rn(b, fields[name])
    return b, fields


def build_sed(phase, isn=0, np_=6, nb=4, nc=1, nzh=8):
    rng = np.random.default_rng(100 + phase)
    b = b""
    for v in (MAGIC_SED, 1, phase, 42, 3, 4, 1, isn, 1, np_, nb, nc, 2, 5, 2, 6):
        b = w_i0(b, v)
    b = w_r0(b, 1.0)
    b = w_i1(b, np.full(nb * nc, 2))
    b = w_i1(b, np.full(nb * nc, 5))
    fields = {"qpv": rng.uniform(size=(np_, nb, nc, nzh))}
    b = w_rn(b, fields["qpv"])
    r1names = ("q_this", "q_other", "qcv", "qtp", "moist_denv", "thetav", "qvv",
               "tv", "dens_col", "momz_col", "u_col", "v_col", "cz_col", "fz_col",
               "dzzmv", "dzvmv")
    for name in r1names:
        fields[name] = rng.uniform(size=nzh)
        b = w_rn(b, fields[name])
    fields["mmass"] = rng.uniform(size=(nb, nc, nzh))
    b = w_rn(b, fields["mmass"])
    for name in ("den_t", "momz_t", "rhou_t", "rhov_t", "rhoe_t"):
        fields[name] = rng.uniform(size=nzh)
        b = w_rn(b, fields[name])
    b = w_r0(b, 3.5)
    fields["sflx"] = 3.5
    return b, fields


def main():
    blob1, f1 = build_micro(1)
    blob2, f2 = build_micro(2)
    blob3, f3 = build_sed(3)
    blob4, f4 = build_sed(4, isn=1)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "amps_dump_r000000_t001.bin"
        p.write_bytes(blob1 + blob3 + blob4 + blob2)
        recs = read_dump_file(p)
    assert len(recs) == 4, f"expected 4 records, got {len(recs)}"
    assert [r["kind"] for r in recs] == ["micro", "sed", "sed", "micro"]
    assert recs[0]["phase"] == 1 and recs[3]["phase"] == 2
    assert recs[0]["TIME_AMPS"] == 42 and recs[0]["i"] == 3 and recs[0]["j"] == 4
    for expected, rec in ((f1, recs[0]), (f2, recs[3]), (f3, recs[1]), (f4, recs[2])):
        for name, val in expected.items():
            got = rec[name]
            assert np.allclose(got, val), f"{rec['kind']} field {name} mismatch"
    assert recs[1]["isn"] == 0 and recs[2]["isn"] == 1
    print("OK: 4 records round-tripped")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 scripts/test_amps_dump_reader.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'amps_dump_reader'`

- [ ] **Step 3: Write the reader**

`scripts/amps_dump_reader.py`:

```python
#!/usr/bin/env python3
"""Parse AMPS reference-data dump files written by the instrumented
scale_atmos_phy_mp_amps.F90 (branch cloudlab_port).

Binary layout v1: sequence of records; ints int32 LE, reals float64 LE,
arrays prefixed by one int32 size per rank, Fortran (column-major) order.
Record types: micro (magic 1095586131) and sed (magic 1095586132) — field
order defined in the M0 plan Tasks 2-3 and mirrored here.

Usage: python amps_dump_reader.py DUMP_DIR -o out.npz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

MAGIC_MICRO = 1095586131
MAGIC_SED = 1095586132

MICRO_HEADER = ("phase", "TIME_AMPS", "i", "j", "isect", "nmic",
                "npr", "nbr", "ncr", "npi", "nbi", "nci", "npa", "nba", "nca",
                "mxnbin", "istrt", "jseed", "ifrst", "isect_seed", "nextn")
MICRO_R1 = ("qcvm", "v3v", "qvvm", "moist_denvm", "ptotvm", "tvm", "wbvm",
            "trpv_thil", "trpv_qtp")
MICRO_R4 = ("qrpvm", "qipvm", "qapvm")
MICRO_POST_EXTRA = (("dmtendlm", 3), ("dcontendlm", 3), ("dbintendlm", 4))

SED_HEADER = ("phase", "TIME_AMPS", "i", "j", "isect", "isn",
              "iadvv", "np", "nb", "nc", "k1", "k2", "k1m", "k2m")
SED_R1 = ("q_this", "q_other", "qcv", "qtp", "moist_denv", "thetav", "qvv",
          "tv", "dens_col", "momz_col", "u_col", "v_col", "cz_col", "fz_col",
          "dzzmv", "dzvmv")
SED_R1_TAIL = ("den_t", "momz_t", "rhou_t", "rhov_t", "rhoe_t")


class _Cursor:
    def __init__(self, buf: bytes):
        self.buf = buf
        self.pos = 0

    def eof(self) -> bool:
        return self.pos >= len(self.buf)

    def i4(self, n: int = 1):
        out = np.frombuffer(self.buf, dtype="<i4", count=n, offset=self.pos)
        self.pos += 4 * n
        return int(out[0]) if n == 1 else out.copy()

    def f8(self, n: int = 1):
        out = np.frombuffer(self.buf, dtype="<f8", count=n, offset=self.pos)
        self.pos += 8 * n
        return float(out[0]) if n == 1 else out.copy()

    def i1_arr(self) -> np.ndarray:
        n = self.i4()
        return self.i4(n) if n > 0 else np.empty(0, dtype="<i4")

    def rn_arr(self, ndims: int) -> np.ndarray:
        shape = tuple(self.i4() for _ in range(ndims))
        n = int(np.prod(shape))
        flat = self.f8(n) if n > 0 else np.empty(0)
        return np.asarray(flat).reshape(shape, order="F")


def _read_micro(c: _Cursor) -> dict:
    rec: dict = {"kind": "micro"}
    version = c.i4()
    assert version == 1, f"unsupported micro record version {version}"
    for name in MICRO_HEADER:
        rec[name] = c.i4()
    rec["dt"] = c.f8()
    rec["kmicvm"] = c.i1_arr()
    for name in MICRO_R1:
        rec[name] = c.rn_arr(1)
    for name in MICRO_R4:
        rec[name] = c.rn_arr(4)
    if rec["phase"] == 2:
        for name, ndims in MICRO_POST_EXTRA:
            rec[name] = c.rn_arr(ndims)
    return rec


def _read_sed(c: _Cursor) -> dict:
    rec: dict = {"kind": "sed"}
    version = c.i4()
    assert version == 1, f"unsupported sed record version {version}"
    for name in SED_HEADER:
        rec[name] = c.i4()
    rec["dt"] = c.f8()
    rec["k1b"] = c.i1_arr()
    rec["k2b"] = c.i1_arr()
    rec["qpv"] = c.rn_arr(4)
    for name in SED_R1:
        rec[name] = c.rn_arr(1)
    rec["mmass"] = c.rn_arr(3)
    for name in SED_R1_TAIL:
        rec[name] = c.rn_arr(1)
    rec["sflx"] = c.f8()
    return rec


def read_dump_file(path: str | Path) -> list[dict]:
    c = _Cursor(Path(path).read_bytes())
    records = []
    while not c.eof():
        magic = c.i4()
        if magic == MAGIC_MICRO:
            records.append(_read_micro(c))
        elif magic == MAGIC_SED:
            records.append(_read_sed(c))
        else:
            raise ValueError(f"{path}: bad magic {magic} at byte {c.pos - 4}")
    return records


def _key(rec: dict) -> str:
    if rec["kind"] == "micro":
        ph = "pre" if rec["phase"] == 1 else "post"
        return f"micro_t{rec['TIME_AMPS']}_i{rec['i']}_j{rec['j']}_{ph}"
    ph = "pre" if rec["phase"] == 3 else "post"
    return f"sed_t{rec['TIME_AMPS']}_i{rec['i']}_j{rec['j']}_s{rec['isn']}_{ph}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump_dir", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=Path("amps_ref.npz"))
    args = ap.parse_args()
    out: dict = {}
    files = sorted(args.dump_dir.glob("amps_dump_r*_t*.bin"))
    if not files:
        raise SystemExit(f"no amps_dump_r*_t*.bin files in {args.dump_dir}")
    for f in files:
        for rec in read_dump_file(f):
            base = _key(rec)
            for name, val in rec.items():
                if name == "kind":
                    continue
                out[f"{base}_{name}"] = val
    np.savez_compressed(args.output, **out)
    print(f"wrote {args.output} ({len(out)} arrays from {len(files)} files)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 scripts/test_amps_dump_reader.py`
Expected: `OK: 4 records round-tripped`

- [ ] **Step 5: Write the run instructions**

`docs/superpowers/specs/2026-07-16-ref-data-run-instructions.md`:

```markdown
# AMPS Reference-Data Runs — Cluster Instructions (M0)

Branch: `cloudlab_port` of jcanton/scale_amps (instrumented
`scale_atmos_phy_mp_amps.F90`). Build exactly as usual (same SCALE_SYS,
same Makefile flags); the instrumentation is inert unless enabled by
namelist.

## 1. Enable dumps

Add to the `&PARAM_ATMOS_PHY_MP_AMPS_bin` group of the run config
(`scale-rm/test/case/cloudlab/scripts/run.conf` for the warm spin-up;
`restart_run.conf` for the ice seeding run):

    l_amps_dump          = .true.,
    amps_dump_dir        = "./amps_dump",
    amps_dump_step_stride = 300,
    amps_dump_is         = 3,
    amps_dump_ie         = 4,
    amps_dump_js         = 3,
    amps_dump_je         = 4,

Create the directory next to the run: `mkdir -p amps_dump`.

Notes:
- Indices are LOCAL (per-rank, halo-inclusive; first interior point is 3
  with IHALO=JHALO=2). Every rank dumps its own 2x2 box -> 32 ranks x 4
  columns.
- `amps_dump_step_stride=300` with TIME_DT_ATMOS_PHY_MP=1.0 s dumps every
  300 s: 12 dump times over the 1 h spin-up.
- Volume estimate, warm run (npr=6,nbr=40, npi=18,nbi=20, npa=5,nba=1,nca=4):
  <= ~0.35 MB per record pair per cloudy column-call; worst case
  12 x 4 x 2 x 0.35 MB ~ 34 MB per rank, ~1 GB total. Usually much less
  (only cloudy levels are in the packed vectors; micro dumps only fire
  when the column has cloudy points).

## 2. Runs wanted (in order of usefulness)

1. **Warm spin-up** (`run.conf`, ice off in AMPSTASK sense — spin-up config
   as-is): primary M2 validation target.
2. **Seeding run** (`restart_run.conf`, ice on, 40 ice bins): M3/M4 target.
3. **Seeding run with deterministic habit**: copy `AMPSTASK.F`, set
   `ihabit_gm_random=0`, rerun 2. Gives the tight ice comparison target
   (counter-based RNG in the port cannot match the Fortran LCG stream).

## 3. Collect results

    tar czf amps_dump_<runname>.tar.gz amps_dump/

Copy the tarballs back. Locally they are parsed with:

    python3 scripts/amps_dump_reader.py amps_dump/ -o amps_ref_<runname>.npz

## 4. Sanity checks after the run

- `amps_dump/` contains `amps_dump_r<rank>_t<thread>.bin` files with
  nonzero size for ranks whose local 2x2 box had cloud.
- The reader completes without "bad magic" errors — if it fails, the
  record layout in the Fortran patch and `scripts/amps_dump_reader.py`
  have diverged; do not proceed, report back.
```

- [ ] **Step 6: Commit and push**

```bash
cd /Users/jcanton/projects/scale_amps
git add scripts/amps_dump_reader.py scripts/test_amps_dump_reader.py \
        docs/superpowers/specs/2026-07-16-ref-data-run-instructions.md
git commit -m "feat(amps-dump): python dump reader, round-trip test, run instructions"
git push origin cloudlab_port
```

---

### Task 5: icon4py worktree + environment baseline

**Files:**
- Create: worktree `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics` (branch `amps_microphysics` off `main`)

**Interfaces:**
- Produces: working uv environment at the worktree root; all later tasks run there. `WT=/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics` used below.

- [ ] **Step 1: Create worktree and branch**

```bash
git -C /Users/jcanton/projects/icon4py worktree add .worktrees/amps_microphysics -b amps_microphysics main
```

Expected: `Preparing worktree (new branch 'amps_microphysics')` ... `HEAD is now at 91604b053 ...`

- [ ] **Step 2: Install environment**

```bash
cd /Users/jcanton/projects/icon4py/.worktrees/amps_microphysics
uv sync --extra all
```

Expected: resolves and installs (several minutes first time); `gt4py==1.1.11` in the output.

- [ ] **Step 3: Baseline test — muphys stencil test passes**

```bash
cd /Users/jcanton/projects/icon4py/.worktrees/amps_microphysics
uv run --group test --frozen pytest -n0 --benchmark-disable \
  model/atmosphere/subgrid_scale_physics/muphys/tests/muphys/stencil_tests/test_cloud_to_graupel.py -v
```

Expected: PASS (embedded backend default). If this fails, the environment is broken — fix before proceeding (no commit for this task; worktree creation is not a repo change).

---

### Task 6: Scaffold + register the `amps` package

**Files (all under `$WT` = the worktree):**
- Create: `model/atmosphere/subgrid_scale_physics/amps/pyproject.toml`
- Create: `model/atmosphere/subgrid_scale_physics/amps/README.md`
- Create: `model/atmosphere/subgrid_scale_physics/amps/src/icon4py/model/atmosphere/subgrid_scale_physics/amps/__init__.py`
- Create: `model/atmosphere/subgrid_scale_physics/amps/tests/__init__.py`
- Create: `model/atmosphere/subgrid_scale_physics/amps/tests/amps/__init__.py`
- Create: `model/atmosphere/subgrid_scale_physics/amps/tests/amps/unit_tests/__init__.py`
- Create: `model/atmosphere/subgrid_scale_physics/amps/tests/amps/unit_tests/test_package.py`
- Modify: `pyproject.toml` (root: dependencies, tool.uv.sources, workspace members, pytest pythonpath, mypy_path)
- Modify: `tach.toml` (source_roots + module entry)
- Modify: `noxfile.py` (ModelSubpackagePath Literal)

**Interfaces:**
- Produces: importable `icon4py.model.atmosphere.subgrid_scale_physics.amps` with `__version__`; spikes (Tasks 7–11) live in `model/atmosphere/subgrid_scale_physics/amps/spikes/`.

License header block used in EVERY new Python file (exact 7 lines):

```python
# ICON4Py - ICON inspired code in Python and GT4Py
#
# Copyright (c) 2022-2024, ETH Zurich and MeteoSwiss
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
```

- [ ] **Step 1: Write the failing smoke test**

`model/atmosphere/subgrid_scale_physics/amps/tests/amps/unit_tests/test_package.py` (with license header, then):

```python
def test_package_imports_and_has_version():
    import icon4py.model.atmosphere.subgrid_scale_physics.amps as amps

    assert amps.__version__ == "0.2.0"
```

Also create `tests/__init__.py` (license header + pkgutil line, copy muphys):

```python
# Build on-the-fly a (legacy) namespace package for 'tests' using pkgutil
__path__ = __import__("pkgutil").extend_path(__path__, __name__)
```

and `tests/amps/__init__.py`, `tests/amps/unit_tests/__init__.py` (license header only).

- [ ] **Step 2: Run test to verify it fails**

```bash
cd $WT && uv run --group test --frozen pytest -n0 model/atmosphere/subgrid_scale_physics/amps/tests/ -v
```

Expected: FAIL/ERROR — `ModuleNotFoundError` (package not installed) or "file or directory not found".

- [ ] **Step 3: Create the package**

`model/atmosphere/subgrid_scale_physics/amps/pyproject.toml` — copy the muphys one verbatim (`model/atmosphere/subgrid_scale_physics/muphys/pyproject.toml`) with exactly these substitutions:
- `name = "icon4py-atmosphere-amps"`
- `description = "AMPS spectral bin (habit-predicting) microphysics parameterization."`
- bumpversion message: `'Bump icon4py-atmosphere-amps version: {current_version} → {new_version}'`
- bumpversion files filename: `"src/icon4py/model/atmosphere/subgrid_scale_physics/amps/__init__.py"`
- keep `version = "0.2.0"`, `dependencies` (`icon4py-common[io]~=0.2.0`, `numpy>=1.23.3`, `gt4py==1.1.11`, `packaging>=20.0`), ruff `extend = "../../../../pyproject.toml"`, setuptools namespaces find.

`src/icon4py/model/atmosphere/subgrid_scale_physics/amps/__init__.py` — copy the muphys `__init__.py` verbatim (license header, `__version__: Final = "0.2.0"` etc.). No other `__init__.py` under `src/` (PEP 420 namespace packages).

`README.md`:

```markdown
# icon4py-atmosphere-amps

GT4Py port of the AMPS spectral bin (habit-predicting) microphysics scheme
(Hashino & Tripoli), ported from the Fortran implementation embedded in
SCALE (contrib/AMPS). Work in progress; see the design spec in the
scale_amps repository: docs/superpowers/specs/2026-07-16-amps-icon4py-port-design.md.
```

- [ ] **Step 4: Register in root config files**

Root `pyproject.toml`:
- `dependencies` list (line ~90): add `"icon4py-atmosphere-amps~=0.2.0",` after the advection entry (keep alphabetical).
- `[tool.uv.sources]` (line ~419): add `icon4py-atmosphere-amps = {workspace = true}` (alphabetical).
- `[tool.uv.workspace] members` (line ~432): add `"model/atmosphere/subgrid_scale_physics/amps",` before the microphysics entry (alphabetical).
- `[tool.pytest.ini_options] pythonpath` (line ~267): add `"model/atmosphere/subgrid_scale_physics/amps",` in the same position.
- mypy `mypy_path` (line ~188): add `$MYPY_CONFIG_FILE_DIR/model/atmosphere/subgrid_scale_physics/amps/src:` in the same position.

`tach.toml`:
- `source_roots`: add `"model/atmosphere/subgrid_scale_physics/amps/src",` before the microphysics line.
- Add module entry after the muphys one:

```toml
[[modules]]
path = "icon4py.model.atmosphere.subgrid_scale_physics.amps"
depends_on = [{ path = "icon4py.model.common" }]
```

`noxfile.py` `ModelSubpackagePath` Literal (line ~34): add `"atmosphere/subgrid_scale_physics/amps",` before the microphysics entry.

- [ ] **Step 5: Sync and run test to verify it passes**

```bash
cd $WT && uv sync --extra all
uv run --group test --frozen pytest -n0 model/atmosphere/subgrid_scale_physics/amps/tests/ -v
```

Expected: `1 passed`.

- [ ] **Step 6: Run pre-commit**

```bash
cd $WT && uv run --group dev --frozen --isolated pre-commit run --all-files
```

Expected: all hooks pass (license headers auto-inserted if any were missed; rerun until clean; `git add` any hook-modified files).

- [ ] **Step 7: Commit**

```bash
cd $WT
git add model/atmosphere/subgrid_scale_physics/amps pyproject.toml tach.toml noxfile.py uv.lock
git commit -m "feat(amps): scaffold icon4py-atmosphere-amps package"
```

---

### Task 7: Spike harness + Spike A (`as_offset` gather remap + table-gather idiom)

**Files (under `$WT/model/atmosphere/subgrid_scale_physics/amps/spikes/`):**
- Create: `common.py`
- Create: `spike_a_remap_gather.py`
- Create: `.gitignore` (content: `_generated/`)

**Interfaces:**
- Consumes: Task 6 package registration (spikes import `icon4py.model.common`).
- Produces (used by Tasks 8–11): `common.py` API — `NCELLS: int = 4096`, `NLEV: int = 61`, `make_field(rng_or_array) -> gtx.Field` on (CellDim, KDim), `backends() -> dict[str, Backend|None]` (`{"embedded": None, "gtfn_cpu": run_gtfn_cached}`), `time_first_and_steady(fn, n_steady=10) -> tuple[float, float]`, `load_generated_operator(source: str, module_name: str, attr: str)`, `GENDIR` path. Every spike script prints a `RESULT ...` table line per measurement and hard-asserts numerical correctness vs numpy.

All spike files start with the license header (Task 6 block).

- [ ] **Step 1: Write `common.py`**

```python
"""Shared helpers for M0 feasibility spikes. Not part of the amps package API."""

from __future__ import annotations

import importlib.util
import pathlib
import sys
import time

import gt4py.next as gtx
import numpy as np
from gt4py.next.program_processors.runners.gtfn import run_gtfn_cached

from icon4py.model.common import dimension as dims


NCELLS = 4096
NLEV = 61
GENDIR = pathlib.Path(__file__).parent / "_generated"


def backends() -> dict:
    return {"embedded": None, "gtfn_cpu": run_gtfn_cached}


def make_field(array: np.ndarray) -> gtx.Field:
    dim_map = {2: (dims.CellDim, dims.KDim), 1: (dims.KDim,)}
    return gtx.as_field(dim_map[array.ndim], array)


def zeros_field(shape: tuple[int, ...] = (NCELLS, NLEV)) -> gtx.Field:
    return make_field(np.zeros(shape))


def time_first_and_steady(fn, n_steady: int = 10) -> tuple[float, float]:
    """Return (first_call_seconds, steady_state_seconds). First call includes
    toolchain compilation for compiled backends."""
    t0 = time.perf_counter()
    fn()
    first = time.perf_counter() - t0
    t0 = time.perf_counter()
    for _ in range(n_steady):
        fn()
    steady = (time.perf_counter() - t0) / n_steady
    return first, steady


def load_generated_operator(source: str, module_name: str, attr: str):
    """Write generated DSL source to a real file (gt4py parses inspect.getsource)
    and import the named attribute from it."""
    GENDIR.mkdir(exist_ok=True)
    path = GENDIR / f"{module_name}.py"
    path.write_text(source)
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return getattr(mod, attr)
```

- [ ] **Step 2: Write `spike_a_remap_gather.py`**

```python
"""Spike A: gather-formulated bin remap and value-indexed table lookup via
as_offset. Answers: does as_offset(Koff, int_field) work on embedded and
gtfn_cpu, for (a) field self-gather with a computed shift, (b) a K-stored
table gathered at an index computed from field values?

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_a_remap_gather.py
"""

from __future__ import annotations

import gt4py.next as gtx
import numpy as np
from gt4py.next import astype, maximum, minimum
from gt4py.next.experimental import as_offset

import common
from icon4py.model.common import dimension as dims, field_type_aliases as fa, type_alias as ta
from icon4py.model.common.dimension import Koff


@gtx.field_operator
def _gather_shift(
    q: fa.CellKField[ta.wpfloat], shift: fa.CellKField[gtx.int32]
) -> fa.CellKField[ta.wpfloat]:
    return q(as_offset(Koff, shift))


@gtx.field_operator
def _table_gather(
    t: fa.CellKField[ta.wpfloat],
    table: gtx.Field[gtx.Dims[dims.KDim], ta.wpfloat],
    k_index: gtx.Field[gtx.Dims[dims.KDim], gtx.int32],
) -> fa.CellKField[ta.wpfloat]:
    idx = astype(maximum(0.0, minimum(t, 59.0)), gtx.int32)
    return table(as_offset(Koff, idx - k_index))


def run_gather_shift(backend) -> None:
    rng = np.random.default_rng(42)
    q_np = rng.uniform(size=(common.NCELLS, common.NLEV))
    k = np.arange(common.NLEV)
    shift_np = rng.integers(-k, common.NLEV - 1 - k, size=(common.NCELLS, common.NLEV))
    q = common.make_field(q_np)
    shift = gtx.as_field((dims.CellDim, dims.KDim), shift_np.astype(np.int32))
    out = common.zeros_field()
    op = _gather_shift.with_backend(backend) if backend is not None else _gather_shift

    def call():
        op(q, shift, out=out, offset_provider={"Koff": dims.KDim})

    first, steady = common.time_first_and_steady(call)
    expected = np.take_along_axis(q_np, k[None, :] + shift_np, axis=1)
    assert np.allclose(out.asnumpy(), expected), "gather_shift mismatch vs numpy"
    print(f"RESULT gather_shift backend={'embedded' if backend is None else 'gtfn_cpu'} "
          f"first={first:.3f}s steady={steady * 1e3:.2f}ms")


def run_table_gather(backend) -> None:
    rng = np.random.default_rng(43)
    t_np = rng.uniform(0.0, 59.0, size=(common.NCELLS, common.NLEV))
    table_np = np.linspace(100.0, 200.0, common.NLEV)
    t = common.make_field(t_np)
    table = common.make_field(table_np)
    k_index = gtx.as_field((dims.KDim,), np.arange(common.NLEV, dtype=np.int32))
    out = common.zeros_field()
    op = _table_gather.with_backend(backend) if backend is not None else _table_gather

    def call():
        op(t, table, k_index, out=out, offset_provider={"Koff": dims.KDim})

    first, steady = common.time_first_and_steady(call)
    idx = np.clip(t_np, 0.0, 59.0).astype(np.int32)
    assert np.allclose(out.asnumpy(), table_np[idx]), "table_gather mismatch vs numpy"
    print(f"RESULT table_gather backend={'embedded' if backend is None else 'gtfn_cpu'} "
          f"first={first:.3f}s steady={steady * 1e3:.2f}ms")


if __name__ == "__main__":
    for name, backend in common.backends().items():
        run_gather_shift(backend)
        run_table_gather(backend)
    print("SPIKE A: PASS")
```

- [ ] **Step 3: Run it**

```bash
cd $WT && uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_a_remap_gather.py
```

Expected: four `RESULT` lines + `SPIKE A: PASS`. If gtfn_cpu fails (compile error on `as_offset` with a K-only table): capture the exact error into the report (Task 12) and mark the table-gather idiom no-go (fallback: analytic fits / numpy escape); the field self-gather and analytic paths are independent.

- [ ] **Step 4: Commit**

```bash
cd $WT && git add model/atmosphere/subgrid_scale_physics/amps/spikes
git commit -m "feat(amps): spike A - as_offset gather remap + table lookup idiom"
```

---

### Task 8: Spike B — generated 40-bin collection kernel (THE gate)

**Files:**
- Create: `$WT/model/atmosphere/subgrid_scale_physics/amps/spikes/spike_b_collection_codegen.py`

**Interfaces:**
- Consumes: `common.py` (Task 7).
- Produces: measured gtfn_cpu compile times for nbins ∈ {8, 20, 40}; the go/no-go datum for Approach A codegen (spec §9 risk 1).

- [ ] **Step 1: Write the spike**

```python
"""Spike B: code-generate an unrolled bin-pair collection (coalescence) kernel
and measure embedded/gtfn_cpu compile + run cost at nbins = 8, 20, 40.

Physics stand-in: Golovin kernel K_ij = k0*(m_i+m_j), mass-doubling bins
m_b = 2^b, Kovetz-Olund style two-bin deposit of coalesced mass, plus loss
terms. Structure (nbins^2 pair terms scattered to <=2 destination bins each,
plus nbins x nbins loss sums) matches the real AMPS coalescence shape.

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_b_collection_codegen.py
"""

from __future__ import annotations

import time

import numpy as np

import common


HEADER = '''\
import gt4py.next as gtx

from icon4py.model.common import field_type_aliases as fa, type_alias as ta

'''


def pair_weights(nbins: int) -> dict:
    """(i, j) -> list of (dest_bin, number_fraction). Mass-doubling grid."""
    mass = 2.0 ** np.arange(nbins)
    weights = {}
    for i in range(nbins):
        for j in range(i, nbins):
            m_new = mass[i] + mass[j]
            d = min(int(np.floor(np.log2(m_new))), nbins - 1)
            if d >= nbins - 1:
                weights[(i, j)] = [(nbins - 1, 1.0)]
            else:
                frac = (mass[d + 1] - m_new) / (mass[d + 1] - mass[d])
                weights[(i, j)] = [(d, frac), (d + 1, 1.0 - frac)]
    return weights


def kernel_table(nbins: int, k0: float = 1.5e-3) -> np.ndarray:
    mass = 2.0 ** np.arange(nbins)
    return k0 * (mass[:, None] + mass[None, :])


def gen_source(nbins: int) -> str:
    kern = kernel_table(nbins)
    weights = pair_weights(nbins)
    args = ",\n    ".join(f"n_{b:02d}: fa.CellKField[ta.wpfloat]" for b in range(nbins))
    rets = ", ".join("fa.CellKField[ta.wpfloat]" for _ in range(nbins))
    gains: dict[int, list[str]] = {b: [] for b in range(nbins)}
    for (i, j), deps in weights.items():
        sym = 1.0 if i == j else 1.0  # symmetric factor folded into K for the spike
        for d, frac in deps:
            gains[d].append(f"{frac * sym * kern[i, j]!r} * n_{i:02d} * n_{j:02d}")
    lines = []
    for b in range(nbins):
        loss = " + ".join(f"{kern[b, j]!r} * n_{j:02d}" for j in range(nbins))
        gain = " + ".join(gains[b]) if gains[b] else "0.0"
        lines.append(f"    dn_{b:02d} = ({gain}) - n_{b:02d} * ({loss})")
    outs = ", ".join(f"dn_{b:02d}" for b in range(nbins))
    return (
        HEADER
        + "@gtx.field_operator\n"
        + f"def _collection_{nbins}(\n    {args},\n) -> tuple[{rets}]:\n"
        + "\n".join(lines)
        + f"\n    return {outs}\n"
    )


def numpy_reference(n: np.ndarray, nbins: int) -> np.ndarray:
    kern = kernel_table(nbins)
    weights = pair_weights(nbins)
    dn = np.zeros_like(n)
    for (i, j), deps in weights.items():
        prod = kern[i, j] * n[i] * n[j]
        for d, frac in deps:
            dn[d] += frac * prod
    dn -= n * np.einsum("bj,j...->b...", kern, n)
    return dn


def run(nbins: int) -> None:
    t0 = time.perf_counter()
    src = gen_source(nbins)
    op = common.load_generated_operator(src, f"gen_collection_{nbins}", f"_collection_{nbins}")
    gen_s = time.perf_counter() - t0
    print(f"RESULT collection nbins={nbins} gen+parse={gen_s:.2f}s "
          f"source_lines={len(src.splitlines())}")

    rng = np.random.default_rng(7)
    n_np = rng.uniform(0.0, 1.0, size=(nbins, common.NCELLS, common.NLEV))
    inputs = [common.make_field(n_np[b]) for b in range(nbins)]
    expected = numpy_reference(n_np, nbins)

    for name, backend in common.backends().items():
        outs = tuple(common.zeros_field() for _ in range(nbins))
        bound = op.with_backend(backend) if backend is not None else op

        def call():
            bound(*inputs, out=outs, offset_provider={})

        first, steady = common.time_first_and_steady(call, n_steady=5)
        got = np.stack([o.asnumpy() for o in outs])
        assert np.allclose(got, expected, rtol=1e-12), f"collection nbins={nbins} {name} wrong"
        print(f"RESULT collection nbins={nbins} backend={name} "
              f"first={first:.1f}s steady={steady * 1e3:.1f}ms")


if __name__ == "__main__":
    for nb in (8, 20, 40):
        run(nb)
    print("SPIKE B: PASS")
```

- [ ] **Step 2: Run it (generous timeout — gtfn compile of the 40-bin kernel is the datum)**

```bash
cd $WT && timeout 3600 uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_b_collection_codegen.py
```

Expected: RESULT lines for nbins=8, 20, 40 on both backends + `SPIKE B: PASS`. Record every number. Interpretation thresholds for the report: gtfn first-call < 120 s for nbins=40 → clean go; 120 s–15 min → go with per-kernel compile caching (run_gtfn_cached persists across runs); > 15 min or compiler failure → no-go for single-operator formulation, note split strategies (per-destination-bin operators) and numpy fallback. If the 40-bin case dies (OOM/compiler crash), keep the 8/20 numbers, record the failure mode verbatim, and try a split variant ONLY if time permits — otherwise report as-is.

- [ ] **Step 3: Commit**

```bash
cd $WT && git add model/atmosphere/subgrid_scale_physics/amps/spikes/spike_b_collection_codegen.py
git commit -m "feat(amps): spike B - generated 40-bin collection kernel measurements"
```

---

### Task 9: Spike C — wide scan carry (sedimentation shape)

**Files:**
- Create: `$WT/model/atmosphere/subgrid_scale_physics/amps/spikes/spike_c_wide_scan.py`

**Interfaces:**
- Consumes: `common.py`.
- Produces: compile/run cost of a scan_operator with a 40-field NamedTuple carry; go/no-go for single-scan sedimentation.

- [ ] **Step 1: Write the spike**

```python
"""Spike C: scan_operator with a 40-field NamedTuple carry — the shape of
40-bin implicit sedimentation (muphys does 4 hydrometeor classes; AMPS needs
40 bins). Generated source (explicit 40-field NamedTuple + 40 scalar args).

Numerics stand-in (per level, forward/downward scan, muphys-style):
  flux_b_new = v_b * q_b * rho ; q_out_b = q_b + zeta*(carry.f_b - flux_b_new)
carry = fluxes. Correctness checked vs numpy recurrence.

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_c_wide_scan.py
"""

from __future__ import annotations

import numpy as np

import common


NBINS = 40
VTS = np.linspace(0.1, 8.0, NBINS)  # per-bin fall speeds, folded as constants


def gen_source(nbins: int) -> str:
    header = (
        "from typing import NamedTuple\n\n"
        "import gt4py.next as gtx\n\n"
        "from icon4py.model.common import dimension as dims, type_alias as ta\n\n\n"
    )
    carry_fields = "\n".join(f"    f_{b:02d}: ta.wpfloat" for b in range(nbins))
    carry = f"class Carry{nbins}(NamedTuple):\n{carry_fields}\n\n\n"
    init = ", ".join(f"f_{b:02d}=0.0" for b in range(nbins))
    qargs = ",\n    ".join(f"q_{b:02d}: ta.wpfloat" for b in range(nbins))
    body_flux = "\n".join(
        f"    fl_{b:02d} = {VTS[b]!r} * q_{b:02d} * rho" for b in range(nbins)
    )
    body_q = "\n".join(
        f"    qo_{b:02d} = q_{b:02d} + zeta * (carry.f_{b:02d} - fl_{b:02d})"
        for b in range(nbins)
    )
    ret_carry = ", ".join(f"f_{b:02d}=fl_{b:02d}" for b in range(nbins))
    outs_tuple = ", ".join(f"s.f_{b:02d}" for b in range(nbins))
    rets = ", ".join("gtx.Field[gtx.Dims[dims.CellDim, dims.KDim], ta.wpfloat]" for _ in range(nbins))
    fargs = ",\n    ".join(
        f"q_{b:02d}: gtx.Field[gtx.Dims[dims.CellDim, dims.KDim], ta.wpfloat]"
        for b in range(nbins)
    )
    fcall = ", ".join(f"q_{b:02d}" for b in range(nbins))
    return (
        header
        + carry
        + f"@gtx.scan_operator(axis=dims.KDim, forward=True, init=Carry{nbins}({init}))\n"
        + f"def _sed_scan_{nbins}(\n    carry: Carry{nbins},\n    {qargs},\n"
        + "    rho: ta.wpfloat,\n    zeta: ta.wpfloat,\n"
        + f") -> Carry{nbins}:\n"
        + body_flux + "\n" + body_q + "\n"
        + f"    return Carry{nbins}({ret_carry})\n\n\n"
        + "@gtx.field_operator\n"
        + f"def _sed_{nbins}(\n    {fargs},\n"
        + "    rho: gtx.Field[gtx.Dims[dims.CellDim, dims.KDim], ta.wpfloat],\n"
        + "    zeta: ta.wpfloat,\n"
        + f") -> tuple[{rets}]:\n"
        + f"    s = _sed_scan_{nbins}({fcall}, rho, zeta)\n"
        + f"    return {outs_tuple}\n"
    )


def numpy_reference(q: np.ndarray, rho: np.ndarray, zeta: float) -> np.ndarray:
    out = np.empty_like(q)
    carry = np.zeros((q.shape[0], q.shape[1]))
    for k in range(q.shape[2]):
        fl = VTS[:, None] * q[:, :, k] * rho[None, :, k]
        out[:, :, k] = fl  # scan output = carry AFTER update = fluxes
        carry = fl
    return out


def main() -> None:
    src = gen_source(NBINS)
    op = common.load_generated_operator(src, f"gen_sed_{NBINS}", f"_sed_{NBINS}")
    rng = np.random.default_rng(11)
    q_np = rng.uniform(0.0, 1e-3, size=(NBINS, common.NCELLS, common.NLEV))
    rho_np = rng.uniform(0.8, 1.2, size=(common.NCELLS, common.NLEV))
    inputs = [common.make_field(q_np[b]) for b in range(NBINS)]
    rho = common.make_field(rho_np)
    expected = numpy_reference(q_np, rho_np, 0.5)

    for name, backend in common.backends().items():
        outs = tuple(common.zeros_field() for _ in range(NBINS))
        bound = op.with_backend(backend) if backend is not None else op

        def call():
            bound(*inputs, rho, 0.5, out=outs, offset_provider={})

        first, steady = common.time_first_and_steady(call, n_steady=5)
        got = np.stack([o.asnumpy() for o in outs])
        assert np.allclose(got, expected, rtol=1e-12), f"wide scan {name} wrong"
        print(f"RESULT wide_scan nbins={NBINS} backend={name} "
              f"first={first:.1f}s steady={steady * 1e3:.1f}ms")
    print("SPIKE C: PASS")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

```bash
cd $WT && timeout 3600 uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_c_wide_scan.py
```

Expected: two RESULT lines + `SPIKE C: PASS`. Note: the scan output-vs-carry convention (does the scan emit the pre- or post-update carry?) is itself a spike finding — if the numpy reference mismatches on embedded, flip the reference recurrence (emit carry BEFORE update), record which convention gt4py uses, and keep the corrected assert.

- [ ] **Step 3: Commit**

```bash
cd $WT && git add model/atmosphere/subgrid_scale_physics/amps/spikes/spike_c_wide_scan.py
git commit -m "feat(amps): spike C - 40-field scan carry measurements"
```

---

### Task 10: Spike D — saturation table vs analytic Murphy–Koop

**Files:**
- Create: `$WT/model/atmosphere/subgrid_scale_physics/amps/spikes/spike_d_esat.py`

**Interfaces:**
- Consumes: `common.py`, table-gather idiom from Spike A.
- Produces: accuracy + speed comparison driving the spec decision "analytic Murphy–Koop replaces estbar LUT".

- [ ] **Step 1: Write the spike**

```python
"""Spike D: saturation vapor pressure over liquid — analytic Murphy-Koop in
DSL vs Fortran-style table + linear interpolation (estbar: 150 entries,
T index = int(T)-163, clamped). Decides whether the port replaces the
Fortran LUTs with analytic formulas (spec 4, core/thermo.py).

Fortran reference: estbar(i) tabulates Murphy-Koop at T = 163+i K
(i = 1..150), linear interp between entries; see QSPARM2 in
mod_amps_utility.F90 (scale_amps repo).

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_d_esat.py
"""

from __future__ import annotations

import gt4py.next as gtx
import numpy as np
from gt4py.next import astype, exp, log, maximum, minimum, tanh
from gt4py.next.experimental import as_offset

import common
from icon4py.model.common import dimension as dims, field_type_aliases as fa, type_alias as ta
from icon4py.model.common.dimension import Koff


def murphy_koop_np(t: np.ndarray) -> np.ndarray:
    return np.exp(
        54.842763 - 6763.22 / t - 4.210 * np.log(t) + 0.000367 * t
        + np.tanh(0.0415 * (t - 218.8))
        * (53.878 - 1331.22 / t - 9.44523 * np.log(t) + 0.014025 * t)
    )


@gtx.field_operator
def _esat_analytic(t: fa.CellKField[ta.wpfloat]) -> fa.CellKField[ta.wpfloat]:
    return exp(
        54.842763 - 6763.22 / t - 4.210 * log(t) + 0.000367 * t
        + tanh(0.0415 * (t - 218.8))
        * (53.878 - 1331.22 / t - 9.44523 * log(t) + 0.014025 * t)
    )


@gtx.field_operator
def _esat_table(
    t: fa.CellKField[ta.wpfloat],
    table: gtx.Field[gtx.Dims[dims.KDim], ta.wpfloat],
    k_index: gtx.Field[gtx.Dims[dims.KDim], gtx.int32],
) -> fa.CellKField[ta.wpfloat]:
    ti = maximum(1.0, minimum(t - 163.0, 149.0))
    i0 = astype(ti, gtx.int32)
    frac = ti - astype(i0, ta.wpfloat)
    e0 = table(as_offset(Koff, i0 - 1 - k_index))
    e1 = table(as_offset(Koff, i0 - k_index))
    return e0 + frac * (e1 - e0)


def main() -> None:
    nlev = 200  # table needs >= 150 K-levels for the gather idiom
    rng = np.random.default_rng(5)
    t_np = rng.uniform(180.0, 310.0, size=(common.NCELLS, nlev))
    t = gtx.as_field((dims.CellDim, dims.KDim), t_np)
    table_np = np.zeros(nlev)
    table_np[:150] = murphy_koop_np(np.arange(1, 151) + 163.0)
    table = gtx.as_field((dims.KDim,), table_np)
    k_index = gtx.as_field((dims.KDim,), np.arange(nlev, dtype=np.int32))
    exact = murphy_koop_np(t_np)

    for name, backend in common.backends().items():
        out_a = gtx.as_field((dims.CellDim, dims.KDim), np.zeros_like(t_np))
        op_a = _esat_analytic.with_backend(backend) if backend is not None else _esat_analytic
        first_a, steady_a = common.time_first_and_steady(
            lambda: op_a(t, out=out_a, offset_provider={})
        )
        assert np.allclose(out_a.asnumpy(), exact, rtol=1e-12)
        print(f"RESULT esat_analytic backend={name} first={first_a:.2f}s "
              f"steady={steady_a * 1e3:.2f}ms")

        out_t = gtx.as_field((dims.CellDim, dims.KDim), np.zeros_like(t_np))
        op_t = _esat_table.with_backend(backend) if backend is not None else _esat_table
        first_t, steady_t = common.time_first_and_steady(
            lambda: op_t(t, table, k_index, out=out_t, offset_provider={"Koff": dims.KDim})
        )
        rel = np.abs(out_t.asnumpy() - exact) / exact
        print(f"RESULT esat_table backend={name} first={first_t:.2f}s "
              f"steady={steady_t * 1e3:.2f}ms max_rel_err_vs_analytic={rel.max():.2e}")
    print("SPIKE D: PASS")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

```bash
cd $WT && uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_d_esat.py
```

Expected: four RESULT lines + `SPIKE D: PASS`. The `max_rel_err_vs_analytic` number quantifies what the Fortran table interpolation itself loses — important for setting per-call validation tolerances (the Fortran uses the TABLE, so exact matching of Fortran requires the table path; the report must state this trade-off explicitly).

- [ ] **Step 3: Commit**

```bash
cd $WT && git add model/atmosphere/subgrid_scale_physics/amps/spikes/spike_d_esat.py
git commit -m "feat(amps): spike D - esat table vs analytic accuracy/speed"
```

---

### Task 11: Spike E — counter-based RNG in DSL integer ops

**Files:**
- Create: `$WT/model/atmosphere/subgrid_scale_physics/amps/spikes/spike_e_counter_rng.py`

**Interfaces:**
- Consumes: `common.py`. Constraint (Global Constraints): only `+ - * // %` on integers — no bitwise, no shifts.
- Produces: a working `hash01(cell, k, bin_id, step) -> uniform [0,1)` DSL field_operator + statistical quality numbers; go/no-go for counter-based habit RNG (spec decision "Counter-based RNG").

- [ ] **Step 1: Write the spike**

```python
"""Spike E: counter-based RNG expressible in the GT4Py DSL, which allows only
+ - * // % on integers (no bitwise ops, no shifts — verified against gt4py
1.1.11 func_to_foast). Construction: Lehmer/LCG rounds modulo the Mersenne
prime 2^31-1 over an int64 counter mixed from (cell, k, bin, step).

Quality bar (for habit selection, NOT crypto): mean ~ 0.5, var ~ 1/12,
lag-1 correlation across each axis < 0.01, coarse 16-bin histogram flat
within 1%. All checked vs numpy replica (bit-identical integers).

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_e_counter_rng.py
"""

from __future__ import annotations

import gt4py.next as gtx
import numpy as np
from gt4py.next import astype

import common
from icon4py.model.common import dimension as dims, field_type_aliases as fa, type_alias as ta


M31 = 2147483647  # 2^31 - 1 (Mersenne prime, Park-Miller modulus)
A1 = 16807        # Park-Miller multipliers for successive rounds
A2 = 48271
A3 = 69621
# mixing multipliers for the counter (distinct odd primes, < 2^20 to keep
# products of int32-range counters safely inside int64)
C_CELL = 999983
C_K = 424243
C_BIN = 786433
C_STEP = 611953


@gtx.field_operator
def _hash01(
    cell_id: gtx.Field[gtx.Dims[dims.CellDim], gtx.int64],
    k_id: gtx.Field[gtx.Dims[dims.KDim], gtx.int64],
    bin_id: gtx.int64,
    step: gtx.int64,
) -> fa.CellKField[ta.wpfloat]:
    x = (cell_id * C_CELL + k_id * C_K + bin_id * C_BIN + step * C_STEP + 1) % M31
    x = (x * A1) % M31
    x = (x * A2) % M31
    x = (x * A3) % M31
    return astype(x, ta.wpfloat) / 2147483647.0


def numpy_replica(ncells: int, nlev: int, bin_id: int, step: int) -> np.ndarray:
    cell = np.arange(ncells, dtype=np.int64)[:, None]
    k = np.arange(nlev, dtype=np.int64)[None, :]
    x = (cell * C_CELL + k * C_K + bin_id * C_BIN + step * C_STEP + 1) % M31
    for a in (A1, A2, A3):
        x = (x * a) % M31
    return x.astype(np.float64) / float(M31)


def main() -> None:
    cell_id = gtx.as_field((dims.CellDim,), np.arange(common.NCELLS, dtype=np.int64))
    k_id = gtx.as_field((dims.KDim,), np.arange(common.NLEV, dtype=np.int64))

    for name, backend in common.backends().items():
        op = _hash01.with_backend(backend) if backend is not None else _hash01
        out = common.zeros_field()

        def call():
            op(cell_id, k_id, gtx.int64(7), gtx.int64(1234), out=out, offset_provider={})

        first, steady = common.time_first_and_steady(call)
        got = out.asnumpy()
        assert np.allclose(got, numpy_replica(common.NCELLS, common.NLEV, 7, 1234)), \
            f"hash01 {name} != numpy replica"
        print(f"RESULT counter_rng backend={name} first={first:.2f}s "
              f"steady={steady * 1e3:.2f}ms")

    # statistical quality on a larger sample across steps (numpy replica,
    # bit-identical to the DSL version as just asserted)
    sample = np.concatenate(
        [numpy_replica(common.NCELLS, common.NLEV, b, s).ravel()
         for b in range(4) for s in range(4)]
    )
    mean, var = sample.mean(), sample.var()
    grid = numpy_replica(common.NCELLS, common.NLEV, 7, 1234)
    lag1_cell = np.corrcoef(grid[:-1].ravel(), grid[1:].ravel())[0, 1]
    lag1_k = np.corrcoef(grid[:, :-1].ravel(), grid[:, 1:].ravel())[0, 1]
    hist = np.histogram(sample, bins=16, range=(0.0, 1.0))[0] / sample.size
    hist_dev = np.abs(hist - 1.0 / 16).max()
    print(f"RESULT rng_quality mean={mean:.4f} var={var:.4f} "
          f"lag1_cell={lag1_cell:.4f} lag1_k={lag1_k:.4f} hist_dev={hist_dev:.4f}")
    assert abs(mean - 0.5) < 0.005, "mean off"
    assert abs(var - 1.0 / 12.0) < 0.005, "variance off"
    assert abs(lag1_cell) < 0.01 and abs(lag1_k) < 0.01, "lag-1 correlation too high"
    assert hist_dev < 0.01, "histogram not flat"
    print("SPIKE E: PASS")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

```bash
cd $WT && uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_e_counter_rng.py
```

Expected: RESULT lines + `SPIKE E: PASS`. If a quality assert fails, add one more LCG round (`A4 = 40692`) and rerun; record final round count. If int64 scalars/fields hit a DSL typing error, record the exact error and test the `int32` variant (counters < 2^15 per axis keep products in range) — the finding matters more than the pass.

- [ ] **Step 3: Commit**

```bash
cd $WT && git add model/atmosphere/subgrid_scale_physics/amps/spikes/spike_e_counter_rng.py
git commit -m "feat(amps): spike E - counter-based RNG in DSL integer ops"
```

---

### Task 12: Spike report + go/no-go + fallback list

**Files:**
- Create: `/Users/jcanton/projects/scale_amps/docs/superpowers/specs/2026-07-16-m0-spike-report.md`
- Modify (if needed): nothing in icon4py — report lives with the spec in scale_amps.

**Interfaces:**
- Consumes: every `RESULT`/failure line from Tasks 7–11 (they are in the terminal transcripts; re-run any script whose numbers were lost).
- Produces: the M0 gate document that M1's plan will cite.

- [ ] **Step 1: Write the report**

Structure (fill every `<...>` with measured values — placeholders may not survive into the commit):

```markdown
# M0 Spike Report — GT4Py feasibility for the AMPS port

Date: <run date>. Environment: gt4py 1.1.11, Python 3.12, backends
embedded + gtfn_cpu (run_gtfn_cached), Apple Silicon macOS (compile times
are indicative, not cluster-representative).

## Verdict table

| Pattern | Spike | embedded | gtfn_cpu first/steady | Verdict |
|---|---|---|---|---|
| as_offset self-gather (bin remap) | A | <ms> | <s>/<ms> | <go/no-go> |
| K-table value-indexed gather (LUTs) | A | <ms> | <s>/<ms> | <go/no-go> |
| Generated collection kernel, nbins=8/20/40 | B | <ms> | <s>/<ms> per size | <go/conditional/no-go> |
| 40-field NamedTuple scan carry | C | <ms> | <s>/<ms> | <go/no-go> |
| Analytic Murphy-Koop vs table | D | <ms both> | <s>/<ms both> | <which one + why> |
| Counter RNG (int64 LCG, no bitwise) | E | <ms> | <s>/<ms> | <go/no-go> |

## Compile-time scaling (Spike B)

nbins=8: <s>; nbins=20: <s>; nbins=40: <s>. Extrapolation comment for the
real coalescence kernel (which adds per-pair property transfer terms —
expect ~<factor>x the spike's expression count): <assessment against the
120 s / 15 min thresholds from the plan>.

## Findings & decisions

1. Codegen design consequence: <one-operator-per-process vs
   per-destination-bin split, based on B>.
2. Scan carry convention: gt4py scan emits <pre|post>-update carry
   (Spike C finding) — codegen for sedimentation must <consequence>.
3. LUT policy: <analytic vs table-gather per table class; note that exact
   Fortran matching requires table semantics — tolerance implication for
   per-call validation>.
4. RNG: <construction that passed, number of LCG rounds, quality numbers>.
5. Failures/surprises verbatim: <every error message hit, or "none">.

## Seeded fallback list (numpy-escape candidates going into M1+)

| Kernel | Reason | Trigger measured in |
|---|---|---|
| <e.g. coalescence 40-bin single operator> | <compile time> | Spike B |
| <none, if all clean> | | |

## Gate decision

M1 may start with Approach A (codegen unrolling) as planned: <yes/no + one
sentence>. Deviations from the spec required: <list or "none">.
```

- [ ] **Step 2: Verify no placeholders remain**

Run: `grep -n '<' /Users/jcanton/projects/scale_amps/docs/superpowers/specs/2026-07-16-m0-spike-report.md | grep -v '^.*http'`
Expected: no output (all angle-bracket placeholders replaced).

- [ ] **Step 3: Commit and push (scale_amps), final icon4py commit**

```bash
cd /Users/jcanton/projects/scale_amps
git add docs/superpowers/specs/2026-07-16-m0-spike-report.md
git commit -m "docs: M0 spike report - GT4Py feasibility results and gate decision"
git push origin cloudlab_port
```

If any spike scripts were modified during report writing (corrected asserts, extra rounds):

```bash
cd /Users/jcanton/projects/icon4py/.worktrees/amps_microphysics
git add model/atmosphere/subgrid_scale_physics/amps/spikes
git commit -m "fix(amps): spike adjustments found during report writing"
```

---

## Plan Self-Review (done at write time)

- **Spec coverage:** M0 items from spec §8: worktree+branch (T5), package skeleton registered (T6), spike a-collection (T8), b-wide-scan (T9→T:C), c-as_offset remap (T7), d-LUT vs analytic (T10), e-counter RNG (T11), spike report + fallback list (T12), instrumentation patch + run instructions (T1–T4). All covered.
- **Layout sync risk:** record layouts are defined once in Task 2/3 Interfaces and mirrored in Task 4 code; Task 3 Step 2 carries an explicit rank-verification instruction with the sync rule.
- **Type consistency:** writer names `AMPS_DUMP_w_*` consistent across T1–T3; `common.py` API names consistent across T7–T11 (`make_field`, `zeros_field`, `backends`, `time_first_and_steady`, `load_generated_operator`); reader function `read_dump_file` matches test import.
- **Known uncertainty, by design:** spikes are measurements — several steps define what to record on failure instead of pretending success is guaranteed (B step 2, C step 2, E step 2, A step 3). That is their purpose; executors must record, not fix silently.
```
