# AMPS Warm-Phase Orchestration — Verbatim Extraction

Repo: `/Users/jcanton/projects/scale_amps`. All paths absolute. Line numbers from `cat -n`.

Warm-phase relevance is annotated inline with `[WARM]`, `[ICE-ONLY skip]`, or `[MIXED]` markers next to each process call.

---

## 1. `cal_micro_tendency` — `class_Cloud_Micro.F90:788-1423` (FULL)

This is the operator-split substep driver. Structure: outer `col_loop1` over `n_step_cl` collision substeps (dt = `dt_cl`), inner `vap_loop` over `n_step_vp` vapor substeps (dt = `dt_vp`). Each substep opens with the **refresh preamble** (`update_mesrc` → `diag_t` → `update_airgroup` → `diag_pq` for ice/rain/aerosol) which is repeated 4× (col-loop it_cl>1 block at 906-916, the `mv_ice2liq` block at 932-953, the vap-loop head at 1206-1214, and the final post-loop refresh at 1374-1381).

```fortran
  subroutine cal_micro_tendency(CM,ID,JD,KD,qtp,thil,iproc,istrt &
                               ,LL,dmtendl,dcontendl,dbintendl &
! <<< 2014/10 T. Hashino added for KiD
!                 ,dM_auto_liq,dM_accr_liq &
!                 ,dM_auto_ice,dM_accr_ice &
!                 ,dM_auto_rim,dM_accr_rim)
!                 ,cptime_mtd)
                               )

    use scale_prc, only: &
       PRC_abort
    use class_AirGroup, only: &
       update_AirGroup
    use class_Group, only: &
       update_group_all, &
       ini_tendency, &
       ini_tendency_ag, &
       update_mesrc, &
       update_modelvars_all
    use mod_amps_core, only: &
       coalescence, &
       melting_shedding, &
       hydrodyn_breakup, &
       ice_nucleation1, &
       ice_nucleation2, &
       cal_aptact_var8_kc04dep, &
       cal_aptact_var8_vec, &
       vapor_deposition, &
       mv_ice2liq, &
       diag_t, &
       diag_pq
    use mod_amps_check, only: &
       repair
    implicit none
! >>> 2014/10 T. Hashino added for KiD
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate tendencies by cloud microphysical processes
    !
    ! This version marches the collision processes and vapor-related processes
    ! with seperate time step, and update the particle properties as it proceeds.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Cloud_Micro), intent(inout)    :: CM
    integer :: ID(*),JD(*),KD(*)
    real(MP_KIND)    :: qtp(*),thil(*)
    integer,intent(in) :: iproc,istrt,LL
    real(MP_KIND), intent(inout) :: dmtendl(10,2,LL), dcontendl(10,2,LL), dbintendl(7,2,mxnbin,LL)
! <<< 2014/10 T. Hashino added for KiD
    real(PS), dimension(LMAX) :: &
                  dM_auto_liq,dM_accr_liq &
                 ,dM_auto_ice,dM_accr_ice &
                 ,dM_auto_rim,dM_accr_rim
! >>> 2014/10 T. Hashino added for KiD
    !real,dimension(20) :: cptime_mtd
    !
    !real,dimension(20) :: cptime
    !real :: cpsum
    !character(len=20),dimension(20)  :: proname
    integer,dimension(mxntend) :: iupdate_gr,iupdate_gs
    integer,dimension(mxntend,CM%ncat_a) :: iupdate_ga
    integer :: i,it_cl,it_vp!,n

    !
    real(PS),dimension(LMAX) :: T_a_r, RV_r
    ! error message
    integer   :: em
    !real(PS) :: s1,r1

    !integer :: i_gg, j_gg ! CHIARUI
    !type (Group), save :: ggg1, ggg2 ! CHIARUI

    integer :: ierr

    ! initialize
    em=0
    !cptime=0.0


!dbg    write(fid_alog,'("current time",F12.2)') CM%cur_time
    if(debug .and. CM%cur_time.eq.0.0) then
       write(fid_alog,'("micexfg",20I2)') CM%micexfg(1:20)
    endif

    RV_r(1:CM%L)=CM%air%TV(1:CM%L)%rv

!!    n=178
!!    i=7
!!    call denchk1(CM%solid_hydro,i,n,id,jd,kd,'in_caltend')
!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfcol')

    col_loop1: do it_cl=1,CM%n_step_cl
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 1. Collection processes
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! set time step size
      CM%rain%dt=CM%dt_cl
      CM%solid_hydro%dt=CM%dt_cl
      do i=1,CM%ncat_a
        CM%aerosol(i)%dt=CM%dt_cl
      enddo
      ! set the tendency mask for repair and update
      !                           1 2 3 4 5 6 7 8 9 10 11 12
          iupdate_gr(1:mxntend)=(/0,1,0,1,0,1,0,1,1, 1, 1, 1/)
          iupdate_gs(1:mxntend)=(/0,1,0,1,0,1,0,1,1, 1, 1, 1/)
      do i=1,CM%ncat_a
        iupdate_ga(1:mxntend,i)=(/0,1,0,1,0,0,0,1,1, 1, 1, 1/)
      enddo
      ! initialize all tendencies that are being updated.
      call ini_tendency(CM%rain,iupdate_gr,1)
      call ini_tendency(CM%solid_hydro,iupdate_gs,1)
      do i=1,CM%ncat_a
        call ini_tendency(CM%aerosol(i),iupdate_ga(1,i),1)
      enddo
      call ini_tendency_ag(CM%air,3)

      !call tic(s1,r1)
!      call PROF_rapstart('amps_update_diagnose',3)
      if(it_cl>1) then
        ! update hydrometeor flags
        call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                          CM%mes_rc,CM%ncat_a)
        ! re-analyze temperature field because T is changed by qr and qi.
        call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
            ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
            )

        ! update air group variables
        call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

      endif
      !cptime(1)=cptime(1)+toc(s1,r1)
    ! diagnose physical quantities of each group object
!!c    write(fid_alog,*) "diag"
!!c      write(fid_alog,*) "ck rdsd01",CM%rdsd
      !call tic(s1,r1)
      if( CM%flagp_s > 0 ) &
        call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                      CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                      CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps, &
                      CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
      !cptime(2)=cptime(2)+toc(s1,r1)


!!c      write(fid_alog,*) "ck rdsd02",CM%rdsd
      !call tic(s1,r1)
      if(CM%level_comp>=6.and.it_cl==1) then
!         ! move 75% melt ice particles
        call mv_ice2liq(CM%solid_hydro,CM%rain,CM%aerosol,CM%air,CM%mes_rc)

        ! update hydrometeor flags
        call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                          CM%mes_rc,CM%ncat_a)

        ! re-analyze temperature field because T is changed by qr and qi.
        call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
            ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
            )

        ! update air group variables
        call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

        if( CM%flagp_s > 0 ) &
          call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps, &
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
      endif
      !cptime(3)=cptime(3)+toc(s1,r1)

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfdiagpq2')

      !call tic(s1,r1)
      if( CM%flagp_r > 0 ) &
        call diag_pq( CM%rain, CM%air, CM%level_comp, em, &
                      CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                      CM%eps_ap0(1),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                      CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
      !cptime(4)=cptime(4)+toc(s1,r1)

      !call tic(s1,r1)
      if( CM%flagp_a /= 0 ) then
        do i=1,CM%ncat_a
          call diag_pq( CM%aerosol(i), CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(i),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml,&
                        i,CM%flagp_a,CM%ap_lnsig,CM%ap_mean)
        end do
      end if
      !cptime(5)=cptime(5)+toc(s1,r1)
      if( em /= 0 ) then
        LOG_ERROR("cal_micro_tendency",*) "diag_pg > error!"
        call PRC_abort
      endif
!      call PROF_rapend('amps_update_diagnose',3)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !    PRINT OUT
      !call tic(s1,r1)
      !if(it_cl==1) then
      !  if(CM%micexfg(1)==1.and.mod(int(CM%cur_time*100), CM%T_print_period*100)==0) then

!!c        write(fid_alog,*) "print rain"
      !    if( CM%flagp_r > 0 ) &
      !      call print_out( CM%rain, CM%air, CM%cur_time, CM%output_format,&
      !                      ID,JD,KD,iproc,istrt)
!!c        write(fid_alog,*) "print aerosols"
      !    if( CM%flagp_a /= 0 ) &
      !      call print_out_ap( CM%aerosol, CM%air, CM%cur_time, CM%output_format,&
      !                         ID,JD,KD,iproc,istrt,CM%ncat_a )
!!c        write(fid_alog,*) "print ices"
      !    if( CM%flagp_s > 0 ) &
      !      call print_out( CM%solid_hydro, CM%air, CM%cur_time, CM%output_format,&
      !                      ID,JD,KD,iproc,istrt)
      !  end if
      !endif
      !cptime(20)=cptime(20)+toc(s1,r1)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfcoale1')
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! collision-coalescence process
      ! NOTE:
!!c    write(fid_alog,*) "coalescence"
      ! +++ rain and rain process +++
      !call tic(s1,r1)
!      call PROF_rapstart('amps_collisioncoal',3)
      if(CM%micexfg(2)==1.and.CM%flagp_r > 0 ) &
         call coalescence( CM%rain,CM%rain,CM%air,CM%level_comp,CM%coll_level,CM%mes_rc &
                ,CM%micexfg(18) &  ! 2014/10 T. Hashino changed for KID
                ,CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd &
                ,ID,JD,KD &
! <<< 2014/10 T. Hashino added for KiD
                ,CM%adrpdrp,CM%drpdrp &
                ,CM%ahexdrp,CM%hexdrp &
                ,CM%abbcdrp,CM%bbcdrp &
                ,CM%acoldrp,CM%coldrp &
                ,CM%agp1drp,CM%gp1drp &
                ,CM%agp4drp,CM%gp4drp &
                ,CM%agp8drp,CM%gp8drp &
                ,dM_auto_liq,dM_accr_liq)
! >>> 2014/10 T. Hashino added for KiD

!!c         call coalescence( CM%rain,CM%rain,CM%air,CM%level_comp,CM%mes_rc,CM%cur_time,0,&
!!c                 CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd)
      !cptime(6)=cptime(6)+toc(s1,r1)
!      call PROF_rapend('amps_collisioncoal',3)
```
**`coalescence(rain,rain,...)` at 1015-1028 — `[WARM]` rain–rain collision-coalescence, guard `micexfg(2)==1 .and. flagp_r>0`.**

```fortran
      ! +++ auto conversion of cloud droplet bin +++
      !call tic(s1,r1)
!tmp      if(CM%micexfg(12)==1.and.CM%flagp_r > 0 ) &
!tmp           call auto_conversion( CM%rain,CM%air,CM%level_comp,CM%mes_rc)
      !cptime(7)=cptime(7)+toc(s1,r1)


      ! +++ aggregation process +++
!!c    write(fid_alog,*) "in of agg"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_aggregation',3)
      if(CM%micexfg(3)==1.and.CM%flagp_s > 0 ) &
         call coalescence(CM%solid_hydro,CM%solid_hydro,CM%air,CM%level_comp,CM%coll_level,CM%mes_rc &
                ,0,CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd &
                ,ID,JD,KD &
! <<< 2014/10 T. Hashino added for KiD
                ,CM%adrpdrp,CM%drpdrp &
                ,CM%ahexdrp,CM%hexdrp &
                ,CM%abbcdrp,CM%bbcdrp &
                ,CM%acoldrp,CM%coldrp &
                ,CM%agp1drp,CM%gp1drp &
                ,CM%agp4drp,CM%gp4drp &
                ,CM%agp8drp,CM%gp8drp &
                ,dM_auto_ice,dM_accr_ice)
! >>> 2014/10 T. Hashino added for KiD
      !cptime(8)=cptime(8)+toc(s1,r1)
!      call PROF_rapend('amps_aggregation',3)
!!c    write(fid_alog,*) "out of agg"
!!c    call check_tendency_ap(CM%solid_hydro,CM%aerosol,CM%air,-1)
    ! +++ riming process +++
!!c    write(fid_alog,*) "in of rim rain"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_riming',3)
      if(CM%micexfg(4)==1.and.CM%flagp_s > 0 .and. CM%flagp_r > 0 ) then
         call coalescence( CM%solid_hydro, CM%rain, CM%air, CM%level_comp,CM%coll_level,CM%mes_rc &
                ,0,CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd &
                ,ID,JD,KD &
! <<< 2014/10 T. Hashino added for KiD
                ,CM%adrpdrp,CM%drpdrp &
                ,CM%ahexdrp,CM%hexdrp &
                ,CM%abbcdrp,CM%bbcdrp &
                ,CM%acoldrp,CM%coldrp &
                ,CM%agp1drp,CM%gp1drp &
                ,CM%agp4drp,CM%gp4drp &
                ,CM%agp8drp,CM%gp8drp &
                ,dM_auto_rim,dM_accr_rim)
! >>> 2014/10 T. Hashino added for KiD
      !cptime(9)=cptime(9)+toc(s1,r1)
      endif
!      call PROF_rapend('amps_riming',3)
!!c    write(fid_alog,*) "out of rim rain"

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! melting and shedding process
!!c    write(fid_alog,*) "ms"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_meltshed',3)
      if(CM%micexfg(8)==1.and. CM%flagp_s > 0 .and. CM%flagp_r > 0) &
         call melting_shedding( CM%solid_hydro, CM%air, CM%level_comp, &
         CM%rain, CM%mes_rc)
      !cptime(18)=cptime(18)+toc(s1,r1)
!      call PROF_rapend('amps_meltshed',3)
!!c    call check_tendency_ap(CM%rain,CM%aerosol,CM%air,1,KD,ID,JD)
!!c    call check_tendency_ap(CM%solid_hydro,CM%aerosol,CM%air,1)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c    call ck_mlttend(CM%rain,CM%solid_hydro,ID,JD,KD)


      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! hydrodynamic breakup process
      ! NOTE:
      ! upgrade the density for hydrodynamic breakup process
      ! calculate the volume of a spheroid circumscribing the ice crystals, and
      ! identify aggregates and grauples.
!!c    write(fid_alog,*) "hb"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_hydrobreakup',3)
      if(CM%micexfg(11)==1.and. CM%flagp_r > 0 .and. CM%hbreak_r > 0 ) &
         call hydrodyn_breakup( CM%rain, CM%hbreak_r,CM%mes_rc)
      !cptime(10)=cptime(10)+toc(s1,r1)
      !call tic(s1,r1)
!!c    write(fid_alog,*) "hb ice"
      if(CM%micexfg(9)==1.and. CM%flagp_s > 0 .and. CM%hbreak_s > 0 ) &
           call hydrodyn_breakup( CM%solid_hydro, CM%hbreak_s,CM%mes_rc)
      !cptime(11)=cptime(11)+toc(s1,r1)
!      call PROF_rapend('amps_hydrobreakup',3)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfinuc1')

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Ice nucleation 1
!!c    write(fid_alog,*) "in1"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_icenucleation',3)
      if(CM%micexfg(10)==1.and. CM%flagp_s > 0 ) &
        call Ice_Nucleation1( CM%solid_hydro, CM%rain, CM%aerosol,CM%air, CM%level_comp,CM%mes_rc,&
          CM%APSNAME,CM%nu_aps,CM%M_aps,&!,CM%phi_aps
          CM%ap_sig_cp,CM%ap_mean_cp,CM%CRIC_RN_IMM,CM%frac_dust, &!,CM%cdf_cp_180m0
          CM%micexfg(14),CM%micexfg(15),CM%micexfg(16),CM%micexfg(17),&
          ID,JD,KD, &
          CM%osm_nhs4,CM%osm_sdch)
      !cptime(13)=cptime(13)+toc(s1,r1)
!      call PROF_rapend('amps_icenucleation',3)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfrepair1')

      !call tic(s1,r1)
!      call PROF_rapstart("amps_repair",3)
      call repair(CM%level_comp, CM%air, CM%rain, CM%solid_hydro, CM%ncat_a, CM%aerosol,&
         CM%flagp_r, CM%flagp_s,CM%flagp_a,CM%mes_rc,& !,CM%act_type
         .false.,.true.,&
         iupdate_gr,iupdate_gs,iupdate_ga,&
         qtp,ID,JD,KD,'af_col',it_cl)
      !cptime(12)=cptime(12)+toc(s1,r1)

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfupdategroup1')

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! update the group vars
      !call tic(s1,r1)
      call update_group_all(CM%rain,CM%solid_hydro,CM%aerosol,CM%air, &
                 CM%level_comp,CM%ncat_a,0,CM%mes_rc, &
                 qtp,RV_r,&
                 iupdate_gr,iupdate_gs,iupdate_ga,&
                 CM%flagp_r,CM%flagp_s,CM%flagp_a,CM%eps_ap0,ID,JD,KD)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !cptime(19)=cptime(19)+toc(s1,r1)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call cal_dmtend_scale(dmtendl,dcontendl,dbintendl,CM%L,CM)

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'afupdategroup1')
!      call PROF_rapend("amps_repair",3)
```

**Collision-substep process calls, in order (all inside `col_loop1`):**
| Line | Call | Guard | Warm-phase |
|---|---|---|---|
| 1015-1028 | `coalescence(rain,rain,…)` | `micexfg(2)==1 .and. flagp_r>0` | **`[WARM]`** rain–rain coalescence |
| 1038-1039 | `auto_conversion` | `micexfg(12)==1` | commented-out (`!tmp`) |
| 1047-1059 | `coalescence(solid,solid,…)` | `micexfg(3)==1 .and. flagp_s>0` | `[ICE-ONLY skip]` aggregation |
| 1069-1082 | `coalescence(solid,rain,…)` | `micexfg(4)==1 .and. flagp_s>0 .and. flagp_r>0` | `[ICE-ONLY skip]` riming |
| 1093-1095 | `melting_shedding` | `micexfg(8)==1 .and. flagp_s>0 .and. flagp_r>0` | `[MIXED]` ice→rain melting/shedding (needs ice present) |
| 1113-1114 | `hydrodyn_breakup(rain,…)` | `micexfg(11)==1 .and. flagp_r>0 .and. hbreak_r>0` | **`[WARM]`** rain hydrodynamic breakup |
| 1118-1119 | `hydrodyn_breakup(solid,…)` | `micexfg(9)==1 .and. flagp_s>0 .and. hbreak_s>0` | `[ICE-ONLY skip]` |
| 1131-1137 | `Ice_Nucleation1` | `micexfg(10)==1 .and. flagp_s>0` | `[ICE-ONLY skip]` |
| 1146-1150 | `repair(…, .false.,.true., 'af_col')` | always | repair after collision |
| 1158-1162 | `update_group_all(…,update_vapor=0,…)` | always | group state advance (collision) |
| 1167 | `cal_dmtend_scale` | always | tendency scaling |

Then the vapor sub-loop:

```fortran
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 2. vapor growth related processes
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! set time step size
      CM%rain%dt=CM%dt_vp
      CM%solid_hydro%dt=CM%dt_vp
      do i=1,CM%ncat_a
        CM%aerosol(i)%dt=CM%dt_vp
      enddo
      ! set the tendency mask for repair and update
      !                             1 2 3 4 5 6 7 8 9 10 11 12
            iupdate_gr(1:mxntend)=(/1,0,1,0,0,0,0,0,0, 0, 0, 0/)
            iupdate_gs(1:mxntend)=(/1,0,1,0,0,0,1,0,0, 0, 0, 0/)
      do i=1,CM%ncat_a
        if(i/=2) then
          iupdate_ga(1:mxntend,i)=(/1,0,1,0,0,1,0,0,0, 0, 0, 0/)
        elseif(i==2) then
          iupdate_ga(1:mxntend,i)=(/1,0,0,0,0,0,1,0,0, 0, 0, 0/)
        endif
      enddo

      vap_loop: do it_vp=1,CM%n_step_vp
        ! initialize all tendencies that are being updated.
        call ini_tendency(CM%rain,iupdate_gr,0)
        call ini_tendency(CM%solid_hydro,iupdate_gs,0)
        do i=1,CM%ncat_a
          call ini_tendency(CM%aerosol(i),iupdate_ga(1,i),0)
        enddo
        call ini_tendency_ag(CM%air,1)
        call ini_tendency_ag(CM%air,2)

        !call tic(s1,r1)
!        call PROF_rapstart('amps_update_diagnose',3)
        ! update hydrometeor flags
        call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                          CM%mes_rc,CM%ncat_a)
        ! re-analyze temperature field because T is changed by qr and qi.
        call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
            ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
            )

        ! update air group variables
        call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

        !cptime(1)=cptime(1)+toc(s1,r1)

        !call tic(s1,r1)
        if( CM%flagp_s > 0 ) &
          call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
        !cptime(2)=cptime(2)+toc(s1,r1)

        !call tic(s1,r1)
        if( CM%flagp_r > 0 ) &
          call diag_pq( CM%rain, CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(1),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
        !cptime(4)=cptime(4)+toc(s1,r1)

        !call tic(s1,r1)
        if( CM%flagp_a /= 0 ) then
          do i=1,CM%ncat_a
            call diag_pq( CM%aerosol(i), CM%air, CM%level_comp, em, &
                          CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                          CM%eps_ap0(i),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                          CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml,&
                          i, CM%flagp_a,CM%ap_lnsig,CM%ap_mean)
          end do
        end if
        !cptime(5)=cptime(5)+toc(s1,r1)
        if( em /= 0 ) then
          LOG_ERROR("cal_micro_tendency",*) "diag_pq_vap > error!"
          call PRC_abort
        endif
!        call PROF_rapend('amps_update_diagnose',3)

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Vapor variable advancement and Activation of CCN
        !call tic(s1,r1)
!        call PROF_rapstart('amps_ccn',3)
        if(CM%act_type==2) then
          !
          ! cloud droplets activation from dry aerosols
          !   stable version
          call cal_aptact_var8_vec(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
              ,CM%flagp_a,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%mes_rc&
              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps,CM%CCNMAX&
              ,CM%snrml,ID,JD,KD)

        elseif(CM%act_type==3) then
          ! cloud droplets activation from dry aerosols with saturation adjustment
!org          call cal_aptact_var9(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
!org              ,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%mes_rc&
!org              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps&
!org              ,ID,JD,KD)

        elseif(CM%act_type==4) then
          ! haze activation from dry aerosols with saturation adjustment
!          call cal_aptact_var11(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
!              ,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%mes_rc&
!              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps&
!              ,ID,JD,KD)

        else
          !
          ! haze activation from dry aerosols and
          ! ice nucleation process based on classical nucleation process approach.
          !
          ! requires use of small time steps (less than 1 sec).
          call cal_aptact_var8_kc04dep(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
              ,CM%flagp_a,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%micexfg(19)&
              ,CM%mes_rc&
              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps&
              ,CM%ap_sig_cp,CM%ap_mean_cp,CM%CCNMAX,CM%CRIC_RN_IMM,CM%frac_dust& !,CM%cdf_cp_180m0
              ,CM%snrml,ID,JD,KD)
        endif
        !cptime(14)=cptime(14)+toc(s1,r1)
!        call PROF_rapend('amps_ccn',3)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! vapor deposition process
        ! 1. vapor deposition process on drops first
!!c    write(fid_alog,*) "vr"
        !call tic(s1,r1)
!        call PROF_rapstart('amps_vapordep',3)
        if(CM%micexfg(6)==1.and. CM%flagp_r > 0 ) &
          call vapor_deposition( CM%rain,CM%aerosol, CM%air, CM%level_comp,CM%mes_rc &
                                ,ID,JD,KD &
                                ,CM%vigp,CM%rdsd,CM%ihabit_gm_random)
!!c    call check_tendency_ap(CM%rain,CM%aerosol,CM%air,-1,KD,ID,JD)
        !cptime(15)=cptime(15)+toc(s1,r1)
!!c    write(fid_alog,*) "vs"

        !call tic(s1,r1)
        if(CM%micexfg(7)==1.and. CM%flagp_s > 0 ) &
          call vapor_deposition(CM%solid_hydro,CM%aerosol,CM%air,CM%level_comp,CM%mes_rc &
                               ,ID,JD,KD &
                               ,CM%vigp,CM%rdsd,CM%ihabit_gm_random)
        !cptime(16)=cptime(16)+toc(s1,r1)
!!c    call check_tendency_ap(CM%rain,CM%aerosol,CM%air,0,KD,ID,JD)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Ice nucleation 2
        !   This one include the mode of depositional growth
!!c    write(fid_alog,*) "in2"
        !call tic(s1,r1)
        if(CM%micexfg(10)==1.and. CM%flagp_s > 0 ) &
          call Ice_Nucleation2( CM%solid_hydro, CM%rain, CM%aerosol,CM%air, CM%level_comp,CM%mes_rc &
           ,CM%APSNAME,CM%nu_aps,CM%M_aps &
           ,CM%micexfg(13) &
           ,CM%flagp_a  &
           ,ID,JD,KD &
           ,CM%vigp,CM%rdsd,CM%ihabit_gm_random,CM%frac_dust,CM%nucleation_halflife)
        !write(fid_alog,*) CM%cur_time, "ice2"
        !cptime(17)=cptime(17)+toc(s1,r1)
!        call PROF_rapend('amps_vapordep',3)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! repair routine
!!c    write(fid_alog,*) " front repair"
        !call tic(s1,r1)
!        call PROF_rapstart("amps_repair",3)
        call repair(CM%level_comp, CM%air, CM%rain, CM%solid_hydro, CM%ncat_a, CM%aerosol,&
             CM%flagp_r, CM%flagp_s,CM%flagp_a,CM%mes_rc,& !,CM%act_type
             .true.,.false.,&
             iupdate_gr,iupdate_gs,iupdate_ga,&
             qtp,ID,JD,KD,'af_vap',it_vp)
        !cptime(12)=cptime(12)+toc(s1,r1)
!!c    write(fid_alog,*) "out of repair"
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! update the group vars
        !call tic(s1,r1)
        call update_group_all(CM%rain,CM%solid_hydro,CM%aerosol,CM%air, &
                 CM%level_comp,CM%ncat_a,1,CM%mes_rc, &
                 qtp,RV_r,&
                 iupdate_gr,iupdate_gs,iupdate_ga,&
                 CM%flagp_r,CM%flagp_s,CM%flagp_a,CM%eps_ap0,ID,JD,KD)

        !cptime(19)=cptime(19)+toc(s1,r1)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        call cal_dmtend_scale(dmtendl,dcontendl,dbintendl,CM%L,CM)
!        call PROF_rapend("amps_repair",3)

!!c        stop

      enddo vap_loop
    enddo col_loop1
```

**Vapor-substep process calls, in order (all inside `vap_loop`):**
| Line | Call | Guard | Warm-phase |
|---|---|---|---|
| 1206-1207 | `update_mesrc` | always | refresh preamble |
| 1209-1211 | `diag_t` | always | refresh T from theta_il |
| 1214 | `update_airgroup` | always | refresh air group |
| 1219-1241 | `diag_pq` ice/rain/aerosol | per flagp | refresh diagnostics |
| 1259-1262 | `cal_aptact_var8_vec` | `act_type==2` | **`[WARM]`** CCN→droplet activation (stable) |
| 1284-1289 | `cal_aptact_var8_kc04dep` | `else` branch (default) | **`[WARM]`** haze activation (+ classical ice nucleation); the requested activation entry point |
| 1301-1304 | `vapor_deposition(rain,…)` | `micexfg(6)==1 .and. flagp_r>0` | **`[WARM]`** vapor deposition on drops |
| 1310-1313 | `vapor_deposition(solid,…)` | `micexfg(7)==1 .and. flagp_s>0` | `[ICE-ONLY skip]` ice vapor dep |
| 1323-1329 | `Ice_Nucleation2` | `micexfg(10)==1 .and. flagp_s>0` | `[ICE-ONLY skip]` |
| 1340-1344 | `repair(…, .true.,.false., 'af_vap')` | always | repair after vapor |
| 1352-1356 | `update_group_all(…,update_vapor=1,…)` | always | group state advance (vapor) |
| 1361 | `cal_dmtend_scale` | always | tendency scaling |

Post-loop final diagnostic refresh (drives terminal velocity):

```fortran
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Finally, diagnose the new variables
    ! for calculation of weighted terminal velocity.
    !call tic(s1,r1)
    ! update hydrometeor flags
    call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                      CM%mes_rc,CM%ncat_a)
    ! re-analyze temperature field because T is changed by qr and qi.
    call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
         ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
         )
    ! update air group variables
    call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

    !cptime(1)=cptime(1)+toc(s1,r1)

    !call tic(s1,r1)
    if( CM%flagp_s > 0 ) &
      call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                    CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                    CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                    CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
    !cptime(2)=cptime(2)+toc(s1,r1)

    !call tic(s1,r1)
    if( CM%flagp_r > 0 ) &
      call diag_pq( CM%rain, CM%air, CM%level_comp, em, &
                    CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                    CM%eps_ap0(1),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                    CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
    !cptime(4)=cptime(4)+toc(s1,r1)

    !call tic(s1,r1)
    if( CM%flagp_a /= 0 ) then
      do i=1,CM%ncat_a
        call diag_pq( CM%aerosol(i), CM%air, CM%level_comp, em, &
                      CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                      CM%eps_ap0(i),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                      CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml,&
                      i, CM%flagp_a,CM%ap_lnsig,CM%ap_mean)
      end do
    end if
    !cptime(5)=cptime(5)+toc(s1,r1)
    if( em /= 0 ) then
      LOG_ERROR("cal_micro_tendency",*) "diag_pq_vap > error!"
      call PRC_abort
    endif
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




    !cptime_mtd(1:20)=cptime_mtd(1:20)+cptime(1:20)

  end subroutine cal_micro_tendency
```

### Tendency-mask notes (drive `ini_tendency` / `repair` / `update_group_all`)
- **Collision substep** (889-893): `iupdate_gr = iupdate_gs = (/0,1,0,1,0,1,0,1,1,1,1,1/)`, `iupdate_ga = (/0,1,0,1,0,0,0,1,1,1,1,1/)`. `repair` called with `(.false.,.true.)`; `update_group_all` with `update_vapor=0`.
- **Vapor substep** (1183-1191): `iupdate_gr = (/1,0,1,0,0,0,0,0,0,0,0,0/)`, `iupdate_gs = (/1,0,1,0,0,0,1,0,0,0,0,0/)`, `iupdate_ga(i≠2) = (/1,0,1,0,0,1,0,0,0,0,0,0/)`, `iupdate_ga(i==2) = (/1,0,0,0,0,0,1,0,0,0,0,0/)`. `repair` called with `(.true.,.false.)`; `update_group_all` with `update_vapor=1`.
- The tendency-process index j: 1=vapor, 2=collision/coal, 3=activation, 4=riming, 6=(aerosol vapor), 7=(ice vapor/ice-nuc dep), 8=melting-shedding-ish, 9=Mossop-Hallett, 10=melting_shedding (per debug comments at 2103-2111).

---

## 2. `ifc_cloud_micro` — `mod_amps_lib.F90:39-326` (top-level driver, FULL)

Top-level sequence: `ini_cloud_micro` (180-187) → `reality_check` (202) → `check_water_apmass` (211-212) → `cal_micro_tendency` (226-233) → `return_output` (252-255) → `check_water_apmass` again (261-262) → `update_terminal_vel` (271-272).

```fortran
SUBROUTINE ifc_cloud_micro( &
     CM &
    ,XC, vtcloud  &
    ,XR, NRTYPE,NRBIN, NRCAT &
    ,XS, NSTYPE,NSBIN, NSCAT &
    ,XA, NATYPE,NABIN, NACAT &
    ,RV, DEN, PT, T, W &
    ,L,ID,JD,KD,iproc_t,iproc,istrt &
    ,qtp,thil &
    ,jseed,ifrst,isect_seed &
    ,nextn &
    ,dmtendl,dcontendl,dbintendl &
!    ,mpirank_w,mpirank_e,mpirank_s,mpirank_n &
! added for sheba
!    ,dmtend,tcon,tmass,tvar,dctend_in &
!    ,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth &
!    ,dmtend_z,tcon_z,tmass_z,n1mx &
! end added for sheba
! <<< 2014/10 T. Hashino added for KiD
!    ,dM_auto_liq,dM_accr_liq &
!    ,dM_auto_ice,dM_accr_ice &
!    ,dM_auto_rim,dM_accr_rim &
!    ,cptacc_ifc,cptacc_mtd &
! >>> 2014/10 T. Hashino added for KiD
     )
  !use tic_toc
  use maxdims
  use class_Cloud_Micro, only: &
     Cloud_Micro, &
     ini_cloud_micro, &
     reality_check, &
     check_water_apmass, &
     cal_micro_tendency, &
     return_output, &
     update_terminal_vel
!  use class_group, only: denchk1
!!!  use mod_amps_utility, only: random_genvar
  implicit none

  type (Cloud_Micro),intent(inout)       :: CM
!tmp  type (Cloud_Micro)      :: CM
!tmp  COMMON /ICMPRIV/CM
!!!!$omp threadprivate(/ICMPRIV/)

  ! C: cloud_drop
  ! R: rain
  ! S: solid hydrometeor
  ! NXBIN       : number of bins
  ! NXTYPE      : number of variables that each bin contain
  ! X(IBIN,ITYPE)    : quantity for ITYPE of IBIN-th bin
  ! DXDT(ITYPE) : tendency of the quantity
  ! ITYPE : 1. mixing ratio (g/g)
  !         2. concentration (#/g)
  !         3. volume of circumscribing sphere (cm^3/g)
  !         4. con. weighted a-axis length^3 (cm^3/g)
  !         5. con. weighted c-axis length^3 (cm^3/g)
  !         6. con. weighted d-axis length^3 (cm^3/g)
  !         7. con. weighted r-axis length^3 (cm^3/g)
  !         8. mixing ratio of riming production (g/g)
  !         9. mixing ratio of ice crystals (g/g)
  !        10. mixing ratio of melt water (g/g)
  !        11. heat (erg/g)
  !
!tmp  integer, intent(in)    :: NRBIN, NSBIN,NABIN
!tmp  integer, intent(in)    :: NRTYPE, NSTYPE,NATYPE
!tmp  integer, intent(in)    :: NRCAT, NSCAT,NACAT
  integer    :: NRBIN, NSBIN,NABIN
  integer    :: NRTYPE, NSTYPE,NATYPE
  integer    :: NRCAT, NSCAT,NACAT
  ! length of microphysics grid vector, thread number, flag for initial do loop
  !  iproc_t includes mpi process and openMP threads
  !  iproc   includes only openMP threads
  integer    :: L,iproc_t,iproc,istrt
  real(MP_KIND)  :: XC(L),vtcloud(L)
  real(MP_KIND)  :: XR(NRTYPE,NRBIN,NRCAT,L), XS(NSTYPE,NSBIN,NSCAT,L),XA(NATYPE,NABIN,NACAT,L)
  real(MP_KIND)  :: RV(L), DEN(L), PT(L), T(L), W(L),qtp(L),thil(L)
  integer   :: ID(L),JD(L),KD(L)

  ! random generator vars
!!!  type (random_genvar),intent(in) :: rdsd0
  integer,intent(in) ::  jseed,ifrst,isect_seed,nextn

  ! tendencies output
  real(MP_KIND), intent(inout) :: dmtendl(10,2,L),dcontendl(10,2,L),dbintendl(3,2,mxnbin,L)

  ! process ranks for MPI parallerization
  !integer,intent(in) :: mpirank_w,mpirank_e,mpirank_s,mpirank_n

  ! current model time in second
!tmp  real,intent(in) :: c_time
! added for sheba
!  real :: dmtend(11,9,*),tcon(9,*),tmass(9,*),zstv(*),zz(*),xx(*),yy(*)
!  real :: tvar(4,*),dctend_in(7,*)
!  integer :: nfpt,nfpth,nzp,nxp,nyp
!  integer :: n1mx
!  real :: dmtend_z(15,n1mx,*),tcon_z(2,n1mx,*),tmass_z(2,n1mx,*)
! end added for sheba
! <<< 2014/10 T. Hashino added for KiD
!  real, dimension(*) :: &
!                  dM_auto_liq,dM_accr_liq &
!                 ,dM_auto_ice,dM_accr_ice &
!                 ,dM_auto_rim,dM_accr_rim
! >>> 2014/10 T. Hashino added for KiD
  !real,dimension(20,*) :: cptacc_ifc,cptacc_mtd


  !real,dimension(20) :: cptime
  !real :: cpsum
!tmp  character(len=32),dimension(20)  :: proname
  integer :: i,n

  real(PS) :: s1,r1

!tmp  proname(1)="ini_cloud_micro"
!tmp  proname(2)="reality_check"
!tmp  proname(3)="check_water_apmass 1"
!tmp  proname(4)="cal_micro_tendency"
!tmp  proname(5)="print_cloud_tendency"
!tmp  proname(6)="return_output"
!tmp  proname(7)="check_water_apmass 2"
!tmp  proname(8)="update_terminal_vel"
!tmp  proname(9)="cal_dmtend"


!!c  write(fid_alog,*) "I am inside! nbin",L,KD(1),ID(1),JD(1)
!!c  write(fid_alog,'(4ES15.8)') RV(1),RV(20),T(1),T(20)

  !cptime=0.0

  ! set random number vars
  CM%rdsd%jseed=jseed
  CM%rdsd%ifrst=ifrst
  CM%rdsd%isect_seed=isect_seed
  CM%rdsd%nextn=nextn

  ! 1. Initialize the Cloud_Micro object


  !call tic(s1,r1)
!tmp  call ini_cloud_micro( &
!tmp        XC &
  call ini_cloud_micro( &
        CM &
       ,XC &
       ,XR, NRTYPE,NRBIN, NRCAT &
       ,XS, NSTYPE,NSBIN, NSCAT &
       ,XA, NATYPE,NABIN, NACAT &
       ,RV, DEN, PT, T, W &
       ,L,qtp,ID,JD,KD)
  !cptime(1)=toc(s1,r1)

!!c       n=16
!!c       write(fid_alog,*) "ck af inicloud:",kd(n),id(n),jd(n),n,CM%air%TV(n)%rv,CM%air%TV(n)%t

!!c  write(fid_alog,*) "af make cloud"
!!  n=178
!!  i=7
!!  call denchk1(CM%solid_hydro,i,n,id,jd,kd,'af_inicloud')

  ! 2. check whether the advected prognostic variables are
  !    in realistic, physical range.
  !call tic(s1,r1)

  call reality_check( CM,qtp,ID,JD,KD)
  !cptime(2)=toc(s1,r1)

!!  n=178
!!  i=7
!!  call denchk1(CM%solid_hydro,i,n,id,jd,kd,'af_real')

  !call tic(s1,r1)

  call check_water_apmass(CM,NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
                          XA,XR,XS,qtp,RV,"make_cloud",0,ID,JD,KD)
  !cptime(3)=toc(s1,r1)

!!c  write(fid_alog,*) "af water apmass"
!!c       n=16
!!c       write(fid_alog,*) "ck af waer apmass:",kd(n),id(n),jd(n),n,CM%air%TV(n)%rv,CM%air%TV(n)%t

  ! 3. calculate the tendency of cloud microphysical variables
  !call tic(s1,r1)
!  call PROF_rapstart("amps_tendency",3)
!!  n=178
!!  i=7
!!  call denchk1(CM%solid_hydro,i,n,id,jd,kd,'bf_microtend')

  call cal_micro_tendency(CM,ID,JD,KD,qtp,thil,iproc_t,istrt &
                         ,L,dmtendl,dcontendl,dbintendl &
! <<< 2014/10 T. Hashino added for KiD
!                 ,dM_auto_liq,dM_accr_liq &
!                 ,dM_auto_ice,dM_accr_ice &
!                 ,dM_auto_rim,dM_accr_rim)
!                 ,cptacc_mtd(1,iproc+1))
                         )
! >>> 2014/10 T. Hashino added for KiD

  !cptime(4)=toc(s1,r1)

!!c  write(fid_alog,*) "af cal_mic_tend"

  ! 4. print tendencies
  !call tic(s1,r1)

  !call print_cloud_tendency(CM,ID,JD,KD,iproc_t,istrt)
  !cptime(5)=toc(s1,r1)
!  call PROF_rapend("amps_tendency",3)

!!c  write(fid_alog,*) "bf ro"

  ! 5. return the calculated tendencies to model prognostic variables.
  !call tic(s1,r1)

  call return_output(CM,L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       NATYPE,NABIN,NACAT,&
       RV,XC,XR,XS,XA,&
       qtp)
  !cptime(6)=toc(s1,r1)

!!c  write(fid_alog,*) "bf cktp1"
  !call tic(s1,r1)

  call check_water_apmass(CM,NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       XA,XR,XS,qtp,RV,"return_model",1,ID,JD,KD)
  !cptime(7)=toc(s1,r1)


!!c  write(fid_alog,*) "bf update"

  ! 6. calculate terminal velocity
  !call tic(s1,r1)

  call update_terminal_vel(CM, L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       vtcloud, XR, XS)
  !cptime(8)=toc(s1,r1)

  ! added for sheba
  ! 7. calculate domain integrated variables for analysis
  !call tic(s1,r1)

! the following is for 3D dynamical simulation
!!c  call cal_dmtend(dmtend,tcon,tmass,tvar,dctend_in,iproc, &
!!c         CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT, &
!!c         XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth, &
!!c         mpirank_w,mpirank_e,mpirank_s,mpirank_n)

! the following is for KiD
  !call cal_dmtend_kid(dmtend,tcon,tmass,iproc,&
  !     CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
  !     XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth)

  !if( mod(int(CM%cur_time), 60) == 0 ) then
  !  call cal_dmtend_z(dmtend_z,tcon_z,tmass_z,iproc,&
  !       CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
  !       dM_auto_liq,dM_accr_liq, &
  !       XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth,&
  !       n1mx,&
  !       mpirank_w,mpirank_e,mpirank_s,mpirank_n)
  !endif
  !cptime(9)=toc(s1,r1)

  !cptacc_ifc(1:20,iproc+1)=cptacc_ifc(1:20,iproc+1)+cptime(1:20)

!tmp  write(fid_alog,*) "iproc+1:",iproc+1
!tmp  write(fid_alog,*) "ifc, cptime",cptime(1:20)
!tmp  write(fid_alog,*) "ifc, cptacc",cptacc_ifc(1:20,iproc+1)

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!tmp  if(mod(CM%cur_time,600.0)==0) then
!!!  if(mod(CM%cur_time,1800.0)==0) then
!tmp    write(fid_alog,*) "CPU time (s)"
!tmp    cpsum=sum(CM%cptime_acc_ifc)
!tmp    if(cpsum>0.0) then
!tmp      do i=1,9
!tmp        write(fid_alog,100) proname(i),CM%cptime_acc_ifc(i),CM%cptime_acc_ifc(i)*100.0/cpsum
!tmp      end do
!tmp      write(fid_alog,*) " "
!tmp100     format(A20,":",ES13.4,"s,",F6.2,"%")
!tmp      CM%cptime_acc_ifc(1:20)=0.0
!tmp    end if
!tmp  end if

  ! end added for sheba

!!! stop


END SUBROUTINE ifc_cloud_micro
```

---

## 3. Substep dt setup — `class_Cloud_Micro.F90:380-390` (verbatim)

```fortran
    ! set time step
    time_step = dt_step

    CM%dt_dy = dt_step
    !  time step for collection process
    CM%dt_cl = dt_step/real(n_step_cl,PS_KIND)
    !  time step for vapor deposition process
    !CM%dt_vp = dt_step/real(n_step_vp)
    CM%dt_vp = CM%dt_cl/real(n_step_vp,PS_KIND)
    CM%n_step_cl=n_step_cl
    CM%n_step_vp=n_step_vp
```

Note: `dt_cl = dt_step / n_step_cl` (collision substep), and `dt_vp = dt_cl / n_step_vp` (vapor substep is a subdivision of the collision substep, not of the full dynamic step). So total vapor substeps per host step = `n_step_cl * n_step_vp`, each of length `dt_step/(n_step_cl*n_step_vp)`.

---

## 4a. `update_mesrc` — `class_Group.F90:11150-11220` (verbatim)

Diagnoses the hydrometeor phase flag `mes_rc(n)`: 0 = no water, 1 = vapor only, 2 = rain (warm), 3 = ice only, 4 = mixed rain+ice. Also sets `mark_cm` on each group.

```fortran
  subroutine update_mesrc(ag,gr,gs,ga,mes_rc,nacat)
    ! **********************************************************************
    ! diagnose hydrometeor flags
    ! **********************************************************************
    ! rain group, ice group
    type (Group), intent(inout) :: gr,gs
    type (group), dimension(*),intent(inout) :: ga
    ! thermo variable object
    type (AirGroup), intent(inout)  :: ag
    ! message from reality-check
    integer,dimension(*),intent(inout)   :: mes_rc
    integer,intent(in)  :: nacat
    ! 3D coordinates
    !integer :: ID(*),JD(*),KD(*)
    ! local vars
    integer :: n,i,ic
    real(PS),dimension(LMAX)         :: M_tot, M_v, M_tr, M_ts

    ! +++ initialization +++

    do n=1,ag%L
      mes_rc(n) = 1
      M_v(n) = ag%tv(n)%rv*ag%tv(n)%den
      M_tr(n)=0.0_PS
      M_ts(n)=0.0_PS
      gr%mark_cm(n)=0
      gs%mark_cm(n)=0
    enddo
    do ic=1,nacat
      do n=1,ag%L
        ga(ic)%mark_cm(n)=0
      enddo
    enddo

    do i=1,gr%n_bin
      do n=1,gr%L
        M_tr(n)=M_tr(n)+gr%MS(i,n)%mass(1)
      enddo
    enddo
    do i=1,gs%n_bin
      do n=1,gs%L
        M_ts(n)=M_ts(n)+gs%MS(i,n)%mass(1)
      enddo
    enddo
    do n=1,gs%L
      M_tot(n)=M_v(n)+M_tr(n)+M_ts(n)
    enddo

    do n=1,ag%L
      if( M_tot(n) <= 0.0 ) then
        mes_rc(n) = 0
      else
        if(M_tr(n)>0.0_PS) then
          if(M_ts(n)>0.0_PS) then
            mes_rc(n)=4
          else
            gs%mark_cm(n)=3
            mes_rc(n)=2
          end if
        else
          gr%mark_cm(n)=3
          if(M_ts(n)>0.0_PS) then
            mes_rc(n)=3
          else
            gs%mark_cm(n)=3
            mes_rc(n)=1
          end if
        end if
      endif
    enddo
  end subroutine update_mesrc
```

---

## 4b. `diag_t` — `mod_amps_core.F90:12449-12550` (verbatim)

Diagnoses T from theta_il and hydrometeor specific humidities `qr_0`, `qi_0` (accumulated per phase using `mes_rc` mask; rain mass minus aerosol mass `rmt-rmat`, ice mass minus aerosol mass `imt-imat`). Iterative inversion of the theta_il→T relation with a 253 K threshold branch.

```fortran
  subroutine diag_t(T_a_r,ag,gr,gs,thil,mes_rc,flagp_r,flagp_s)
    use mod_amps_utility, only: cal_growth_mode_inl_vec,random_genvar
    ! **********************************************************************
    ! diagnose temperature from theta il and hydrometeor specific humidity
    ! **********************************************************************
    ! rain group, ice group
    type (Group), intent(inout) :: gr,gs
    ! thermo variable object
    type (AirGroup), intent(inout)  :: ag
    ! message from reality-check
    integer,dimension(*)   :: mes_rc
    ! theta il
    real(MP_KIND),dimension(*)  :: thil
    ! 3D coordinates
    !integer :: ID(*),JD(*),KD(*)
    !
    integer,intent(in)  :: flagp_r,flagp_s
    !
    real(PS),dimension(*),intent(inout) :: T_a_r

    integer :: n,i
    integer,dimension(LMAX) :: ierror1
    real(PS),dimension(LMAX) :: qr_0,qi_0
    real(PS) :: T_a0,til

    do n=1,ag%L
      qr_0(n)=0.0_PS
      qi_0(n)=0.0_PS
    enddo

    if(flagp_s>0) then
      do i=1,gs%N_BIN
        do n=1,ag%L
          if(mes_rc(n)==3.or.mes_rc(n)==4) then
            qi_0(n)=qi_0(n)+&
              max(0.0_PS,gs%MS(i,n)%mass(imt)-gs%MS(i,n)%mass(imat))
          end if
        enddo
      enddo
    endif

    if(flagp_r>0) then
      do i=1,gr%N_BIN
        do n=1,ag%L
          if(mes_rc(n)==2.or.mes_rc(n)==4) then
            qr_0(n)=qr_0(n)+&
              max(0.0_PS,gr%MS(i,n)%mass(rmt)-gr%MS(i,n)%mass(rmat))
          endif
        enddo
      enddo
    endif
    do n=1,ag%L

      qr_0(n)=qr_0(n)/ag%TV(n)%den
      qi_0(n)=qi_0(n)/ag%TV(n)%den

      ! +++ calculate T_il from thil +++
      til=thil(n)*(ag%TV(n)%P/p00)**(Racp)

      T_a0=ag%TV(n)%T
      T_a_r(n)=Til*(1.0_PS+(L_e*qr_0(n)+L_s*qi_0(n))/(c_pa*253.0_PS))

      ierror1(n)=0
      if(T_a_r(n)>=253.0_PS) then
        T_a_r(n)=0.5*(Til+sqrt(Til**2+4.0_PS*Til/c_pa*(L_e*qr_0(n)+L_s*qi_0(n))))
        if(T_a_r(n)<253.0_PS) then
          ierror1(n)=1
        end if
      endif
    enddo
!dbg    if(any(ierror1(1:ag%L)>0)) then
!dbg      do n=1,ag%L
!dbg        if(ierror1(n)>0) then
!dbg          write(*,*) "diag_t>something wrong for T_a_r",T_a_r(n),ag%TV(n)%T &
!dbg                   ,thil(n),ag%TV(n)%P,qr_0(n),qi_0(n)
!dbg        endif
!dbg      enddo
!dbg!     stop
!dbg    endif

!tmp       ! re-initialize thermo variables
!tmp       ag%TV(n)=make_thermo_var3(real(ag%TV(n)%rv,MP),real(ag%TV(n)%den*1.0e+3_PS,MP) &
!tmp              ,real(ag%TV(n)%P*1.0e-1_PS,MP),real(T_a_r,MP),real(ag%TV(n)%W*1.0e-2_PS,MP) &
!tmp              ,ag%estbar,ag%esitbar)
!tmp
!tmp       ! re-analyze the field flag
!tmp       if(qr_0>0.0_PS) then
!tmp          if(qi_0>0.0_PS) then
!tmp             mes_rc(n)=4
!tmp          else
!tmp             mes_rc(n)=2
!tmp          endif
!tmp       else
!tmp          if(qi_0>0.0_PS) then
!tmp             mes_rc(n)=3
!tmp          else
!tmp             mes_rc(n)=1
!tmp          endif
!tmp       endif

!!c       write(*,'("diag_T>k,i,j,inmlt_t,T_new,T_old",4I5,10ES15.6)') KD(n),ID(n),JD(n),inmlt_t,T_a_r,T_a0
  end subroutine diag_t
```

Constants used: `Racp`, `p00`, `L_e` (latent heat evap), `L_s` (latent heat sublimation), `c_pa` (specific heat air). Note the linearized first guess divides by fixed 253 K, then if `T_a_r >= 253` it refines with the quadratic-root formula.

---

## 4c. `update_group_all` — `class_Group.F90:1684-3770` (2086 lines)

This routine is ~2086 lines of per-bin state advancement arithmetic. Below is the full **signature + declaration block** verbatim, the top-level control-flow skeleton with exact line anchors, and the orchestration-relevant blocks verbatim (init, the rain block in full, the vapor-field-update block, and the aerosol block header/end). The rain / ice / aerosol per-bin arithmetic bodies are structurally analogous; I mark their exact line ranges rather than reproduce all inner arithmetic.

**Key orchestration semantics:** `update_vapor` argument = **0 on the collision substep call** (line 1159) and **1 on the vapor substep call** (line 1353). The `update_vapor==1` branch (3349-3420) is what recomputes the air vapor mixing ratio `ag%tv(n)%rv` and `rv2(n)` from mass conservation (`qtp - (totmixr+totmixs)/den`). `dtd` is taken from `gr%dt`/`gs%dt`/`ga%dt`, which the caller sets to `dt_cl` (collision) or `dt_vp` (vapor) before calling.

### Signature + declarations (1684-1802, verbatim)

```fortran
  subroutine update_group_all(gr,gs,ga,ag,&
       level,nacat,update_vapor,mes_rc,&
       qtp,rv2,&
       iupdate_gr,iupdate_gs,iupdate_ga,&
       flagp_r,flagp_s,flagp_a,eps_ap0,ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    use class_Thermo_Var, only: &
       get_sat_vapor_pres_lk
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Note: time increment has to be specified for each group
    !       before this one is called.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! level of complexity
    integer, intent(in)  :: level
    integer,intent(in)  :: update_vapor,nacat
    integer,intent(in) :: ID(*),JD(*),KD(*)
    integer,dimension(*),intent(in) :: iupdate_gr,iupdate_gs
    integer,dimension(mxntend,*),intent(in) :: iupdate_ga
    type (group), intent(inout)  :: gr,gs
    type (group), dimension(*),intent(inout) :: ga
    type (airgroup), intent(inout) :: ag
    ! message from reality-check
    integer,dimension(*)   :: mes_rc
    real(MP_KIND),dimension(*),intent(in) :: qtp
    real(PS),dimension(*),intent(inout) :: rv2
    real(PS),dimension(*),intent(in)  :: eps_ap0
    integer,intent(in)  :: flagp_r,flagp_s,flagp_a
    ! new space
    real(8),dimension(mxnbin,LMAX)     :: ndxdt,ndxdt2
    !real(PS) :: ndxdt3
    real (ps)     :: var_n1,var_n2
    ! original mass concentration
    real(8),dimension(mxnbin,LMAX)  :: om
    ! original number concentration
    real(8),dimension(mxnbin,LMAX)  :: oc
    ! original total volume variables
    real(8),dimension(mxnbin,LMAX)  :: oQ
    real(8),dimension(mxnbin,LMAX)  :: oq_alen!,phi_dbg,psi_dbg
    real(8),dimension(mxnbin,LMAX)  :: oq_clen,ndxdt_alen,ndxdt_clen

    ! mass fraction of aerosol particles to total mass of hydrometeors and
    ! mass fraction of soluble aerosol particles to total aerosol mass
    ! those are used to maintain positive aerosol masses, and use fraction
    ! of smaller bin.
    real(ps),dimension(max(gr%n_bin,gs%n_bin),LMAX) :: fap
    !real(ps) :: fapt_r,faps_r,fapt_s,faps_s, fap3
    ! minimum possible fractions
    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=1.0e-5_ps&
!!c    real(ps),parameter :: min_fapt_r=1.0e-18_ps,min_faps_r=0.0_ps&
         ,min_fapt_s=1.0e-18_ps,min_faps_s=0.0_ps,sep_faps=1.0e-5_PS
    ! possible ratio of mean mass to bin boundaries
    real(ps),parameter :: brat1=1.001,brat2=0.999
    real(ps),parameter :: mlmt=1.0e-30,nlmt=1.0e-30,m_lmt_ap=1.0e-25
    !     the minimum possible background concentration for large particle
    !     is assumed to be 1.0e-5 cm^-3
    !real(ps),parameter :: n_lmt_ap=1.0e-5,n_max_ap=1.0e+4
    real(ps),parameter :: n_lmt_ap=0.0d0,n_max_ap=1.0e+5
!parcel model    real(ps),parameter :: n_lmt_ap=1.0e-15,n_max_ap=1.0e+4
    !     The minimum radius possible for accumulation particles
    real(ps),parameter :: r3_lmt=1.0e-18
    !     The minimum radius possible for nucleation particles
    real(PS),parameter :: r3_lmt_nuc=1.0e-21


!org    real(ps),parameter :: mx_aprat=0.99999
    real(ps),parameter :: mx_aprat=0.95
    real(ps) :: m_lmt
    ! realistic bounds for length predictions
    real(ps),parameter :: min_hexlen=1.0e-4,max_hexlen=3.0,&
                          max_roslen=3.0,max_irrlen=3.0,max_exice=1.0
    ! realistic bounds for crystal mass
    real(ps),parameter :: m_icmin=4.763209003e-12

    real(PS) :: tmass1,tmass2!,tmass3
    real(PS) :: phi,psi,pag,pcg,den_ip1
    real(PS) :: a_min,c_min,d_min,ag_min,cg_min

    ! minimum, and maximum volume of circumscribing sphere that corresponds to 1 um of
    ! hex ice crystal and sphere of 10cm radius
    real(PS),parameter :: V_csmin_hex=1.184768784e-11,V_csmax=4188.79020478639
!!c    real,parameter :: V_csmin_hex=1.184768784e-11,V_csmax=4188.79020478639e+10


    ! maximum and minimum aspect ratio
    real(PS),parameter  :: phi_max=2.0e+1,phi_min=5.0e-3

    ! accuracy limit to prevent dubious crystal features by numerical errors
    real(PS),parameter :: accuracy_lmt=1.0e-6

    ! vapor specific humidity cacluated with activation scheme
    real(PS) :: rv_n
    ! total mixing ratio of hydrometeors
    real(PS),dimension(LMAX) :: totmixr,totmixs
    real(PS),dimension(LMAX) :: s_n
    real(PS) :: e_n

    real(8) :: dtd,con8

!    integer,dimension(mxnbin*LMAX*mxntend) :: ierror
    integer,dimension(mxnbin,LMAX) :: icond1,icond2,icond3 !,ierror2,ierror3

    integer  :: i,j,k,n,ic,L!,ij,var_status
    !integer :: ie
    !integer :: ic1,ic2
    integer :: iq1

    !integer :: imark

    integer  :: maxice_kidx, maxice_bidx, maxice_kidx2, maxice_bidx2
    real(PS) :: maxice_s, maxice_s2, maxice_d(gs%n_tendpros), maxice_m(16)

    L=max(gr%L,gs%L,ga(1)%L)


    do n=1,L
      totmixs(n)=0.0_PS
      totmixr(n)=0.0_PS
    enddo
```

### Top-level control-flow skeleton (line anchors)

| Lines | Block | Purpose |
|---|---|---|
| 1799-1802 | init | zero `totmixr(n)`, `totmixs(n)` |
| **1807-2116** | **`if( flagp_r > 0 )`** | **`[WARM]` RAIN update — reproduced in full below.** advances `mass(rmt)`, `con`, aerosol masses `rmat`/`rmas`, accumulates `totmixr` (2091) |
| 2118-3346 | `if( flagp_s > 0 )` | `[ICE-ONLY skip]` ICE update — analogous per-bin arithmetic on ice mass variables (`imt`, `imr`, `ima`, `imc`, `imw`, `imf`, `imat`, `imas`), ice shape lengths, growth mode; accumulates `totmixs` (3342). Header + max-tracking at 2120-2203. |
| **3348-3420** | **`if(update_vapor==1)`** | **`[WARM-relevant]` vapor field update — reproduced in full below.** recomputes `ag%tv(n)%rv`, `rv2(n)`, saturation `s_n(n)`; clamps negatives / supersaturation |
| 3425-3769 | `if(flagp_a/=0)` do `ic=1,nacat` | AEROSOL update — per-category interstitial aerosol advance; header/cycle logic reproduced below (3425-3436) |
| 3770 | `end subroutine` | — |

### RAIN block — 1807-2116 (verbatim, `[WARM]`)

```fortran
    if( flagp_r > 0 ) then

      dtd=gr%dt

!      do in=1,gr%n_bin*L
!        n=(in-1)/gr%N_BIN+1
!        i=in-(n-1)*gr%N_BIN
      do n = 1, L
      do i = 1, gr%N_BIN
        ndxdt(i,n) = 0.0d+0
        ndxdt2(i,n) = 0.0d+0
      enddo
      enddo

      do j=1,gr%n_tendpros
        if(iupdate_gr(j)==0) cycle
        do n = 1, L
        do i = 1, gr%N_BIN
          ndxdt(i,n) = ndxdt(i,n)+gr%MS(i,n)%dmassdt(rmt,j)
          ndxdt2(i,n) = ndxdt2(i,n)+gr%MS(i,n)%dcondt(j)
        enddo
        enddo
      enddo

      if ( debug ) then
         do n = 1, L
         do i = 1, gr%n_bin
         do j = 1, gr%n_tendpros
            if(iupdate_gr(j)==0) cycle
            if(gr%MS(i,n)%dmassdt(rmt,j)>1.0e+03.or.&
               gr%MS(i,n)%dmassdt(rmt,j)<-1.0e+03) then
               write(*,201) n,i,j,gr%MS(i,n)%dmassdt(rmt,j)
201            format("Warning mt_rain>grid,bin,process,dmassdt",3i5,es15.6)
            endif
         enddo
         enddo
         enddo
         do n = 1, L
         do i = 1, gr%n_bin
         do j = 1, gr%n_tendpros
            if(iupdate_gr(j)==0) cycle
            if(gr%MS(i,n)%dcondt(j)<-1.0e+8 .or. &
               gr%MS(i,n)%dcondt(j)>1.0e+8 ) then
               write(*,*) "Warning con_r:bin,grid,process",i,n,j,gr%MS(i,n)%dcondt(j)
            endif
         enddo
         enddo
         enddo
      end if

      do n = 1, L
      do i = 1, gr%N_BIN
        om(i,n)=gr%MS(i,n)%mass(rmt)
        gr%MS(i,n)%mass(rmt)=gr%MS(i,n)%mass(rmt)+ndxdt(i,n)*dtd

        oc(i,n)=gr%MS(i,n)%con
        gr%MS(i,n)%con=gr%MS(i,n)%con+ndxdt2(i,n)*dtd
      enddo
      enddo
      !
      ! check
      !
      icond1(1:gr%n_bin,1:L)=0
      icond2(1:gr%n_bin,1:L)=0
      do n = 1, L
      do i = 1, gr%N_BIN
        if(gr%MS(i,n)%mass(rmt)<1.0e-30.or.gr%MS(i,n)%mass(rmt)<om(i,n)*1.0e-6) then
          ! case of zero out.  Due to single-presicion of mass variables.
          icond1(i,n)=1
        elseif(gr%MS(i,n)%mass(rmt)>1.0e-1) then
           icond1(i,n)=2
        end if

        if(gr%MS(i,n)%con<1.0e-30.or.gr%MS(i,n)%con<oc(i,n)*1.0e-6) then
          ! case of zero out.  Due to single-presicion of con variable.
          icond2(i,n)=1
        else
          icond2(i,n)=3
        endif
      enddo
      enddo

      do n = 1, L
      do i = 1, gr%N_BIN
        if(icond1(i,n)==1.or.icond2(i,n)==1) then
          gr%MS(i,n)%con=0.0_PS
          gr%MS(i,n)%mass(rmt)=0.0_PS
        endif
      end do
      end do
      do n = 1, L
      do i = 1, gr%N_BIN
        if(icond1(i,n)/=1.and.icond2(i,n)==3) then
          gr%MS(i,n)%con= &
                  gr%MS(i,n)%mass(rmt)&
                  /max(brat1*gr%binb(i),min(gr%MS(i,n)%mass(rmt)/gr%MS(i,n)%con&
                 ,brat2*gr%binb(i+1)))
        endif
      enddo
      enddo

      if ( debug ) then
         do n = 1, L
         do i = 1, gr%N_BIN
            if(icond1(i,n)==2) then
               write(*,251) i,n,gr%MS(i,n)%mass(rmt),gr%MS(i,n)%con
251            format("warning:qr is large at bin,grid:",2i5,2es15.6)
            endif
         end do
         end do
      end if

      if(level>=4) then
         do n = 1, L
         do i = 1, gr%N_BIN
          ndxdt(i,n) = 0.0d+0
          ndxdt2(i,n) = 0.0d+0
        enddo
        enddo
        do j=1,gr%n_tendpros
          if(iupdate_gr(j)==0) cycle
          do n = 1, L
          do i = 1, gr%N_BIN
            ! mass by total mass of aerosols (g/g)
            ndxdt(i,n)=ndxdt(i,n)+gr%MS(i,n)%dmassdt(rmat,j)

            ! mass by soluble mass of aerosols (g/g)
            ndxdt2(i,n)=ndxdt2(i,n)+gr%MS(i,n)%dmassdt(rmas,j)
          enddo
          enddo
        enddo

        if ( debug ) then
           do n = 1, L
           do i = 1, gr%n_bin
           do j = 1, gr%n_tendpros
              if(iupdate_gr(j)==0) cycle
              if(gr%MS(i,n)%dmassdt(rmat,j)>1.0e+08.or.&
                 gr%MS(i,n)%dmassdt(rmat,j)<-1.0e+08)then
                 write(*,236) i,n,j,gr%MS(i,n)%dmassdt(rmat,j)
236              format("Warning aptm in rain>bin,grid,process,dmassdt",3i5,es15.6)
              endif
           enddo
           enddo
           enddo
           do n = 1, L
           do i = 1, gr%n_bin
           do j = 1, gr%n_tendpros
              if(iupdate_gr(j)==0) cycle
              if(gr%MS(i,n)%dmassdt(rmas,j)>1.0e+08.or.&
                 gr%MS(i,n)%dmassdt(rmas,j)<-1.0e+08)then
                 write(*,237) i,n,j,gr%MS(i,n)%dmassdt(rmas,j)
237              format("Warning apsm in solid>bin,grid,process,dmassdt",3i5,es15.6)
              endif
           enddo
           enddo
           enddo
        end if

        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gr%MS(i,n)%mass(rmat)=gr%MS(i,n)%mass(rmat)+ndxdt(i,n)*dtd
            fap(i,n)=gr%MS(i,n)%mass(rmat)/gr%MS(i,n)%mass(rmt)
          endif
        enddo
        enddo
        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            if(fap(i,n)>1.0_PS) then
              gr%MS(i,n)%mass(rmt)=gr%MS(i,n)%mass(rmat)/mx_aprat
              gr%MS(i,n)%con= &
                  gr%MS(i,n)%mass(rmt)&
                  /max(brat1*gr%binb(i),min(gr%MS(i,n)%mass(rmt)/gr%MS(i,n)%con&
                 ,brat2*gr%binb(i+1)))

            elseif(fap(i,n)<=min_fapt_r) then

              gr%MS(i,n)%mass(rmat)=min(max(m_lmt_ap, &
                                   min_fapt_r*gr%MS(i,n)%mass(rmt))&
                                  ,mx_aprat*gr%MS(i,n)%mass(rmt))
            else
              gr%MS(i,n)%mass(rmat)=min(mx_aprat*gr%MS(i,n)%mass(rmt), &
                                    max(m_lmt_ap, &
                                    gr%MS(i,n)%mass(rmat)))

            endif

          endif
        enddo
        enddo
        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            gr%MS(i,n)%mass(rmas)=gr%MS(i,n)%mass(rmas)+ndxdt2(i,n)*dtd
            fap(i,n)=gr%MS(i,n)%mass(rmas)/gr%MS(i,n)%mass(rmat)
          endif
        enddo
        enddo
        do n = 1, L
        do i = 1, gr%N_BIN
          if(icond1(i,n)/=1.and.icond2(i,n)/=1) then
            if(fap(i,n)<min_faps_r) then
              gr%MS(i,n)%mass(rmas)=min(max(m_lmt_ap,&
                                 min_faps_r*gr%MS(i,n)%mass(rmat)) &
                                ,gr%MS(i,n)%mass(rmat))
            endif
          endif
        enddo
        enddo

      endif

      do i=1,gr%n_bin
        do n=1,L
          totmixr(n)=totmixr(n)+max(0.0_PS,gr%MS(i,n)%mass(rmt)-gr%MS(i,n)%mass(rmat))
        enddo
      end do

      if(debug) then
        do n=1,L
         if(sum(gr%MS(1:gr%n_bin,n)%con)>0.0) then
         write(*,'("update_group:bin,i,j,k,total rcon, others",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
                      ,sum(gr%MS(1:gr%n_bin,n)%con) &
                      ,gr%MS(1:gr%n_bin,n)%con

         write(*,'("update_group:bin,i,j,k,total act, total vap mass tendency",4i5,100es15.6)') n,ID(n),JD(n),KD(n) &
                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,3)) &  ! activation
                      ,sum(gr%MS(1:gr%n_bin,n)%dmassdt(rmt,1))    ! vapor
         write(*,*) "cal_model_tendency:con-bin,i,j,k,total rim, total melt, total mossop, total coal"
         write(*,'(4I5,10ES15.6)') n,ID(n),JD(n),KD(n) &
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(4)) &  ! riming
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(10)) &  ! melting_shedding
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(9)) &  ! mossop and hallet
                      ,sum(gr%MS(1:gr%n_bin,n)%dcondt(2))    ! collision-coal
         endif
       enddo
      endif

    endif
```

### Vapor-field update — 3346-3420 (verbatim, `[WARM-relevant]`, runs only when `update_vapor==1`)

```fortran
    endif

    ! +++ update vapor field +++
    if(update_vapor==1) then

      dtd=max(gr%dt,gs%dt,ga(1)%dt)
      do n=1,L
        ag%tv(n)%rv=ag%tv(n)%rv+ag%tv(n)%dmassdt/ag%tv(n)%den*dtd
        rv2(n)=qtp(n)-(totmixr(n)+totmixs(n))/ag%TV(n)%den
        e_n=ag%TV(n)%P*ag%tv(n)%rv/(Rdvchiarui+ag%tv(n)%rv)
        s_n(n)=e_n/get_sat_vapor_pres_lk(1,ag%TV(n)%T_m,ag%estbar,ag%esitbar)-1.0

!!c        write(*,*) "ck rv",kd(n),id(n),jd(n),ag%TV(n)%rv,ag%TV(n)%dmassdt, &
!!c                    rv2(n),qtp(n),totmixr(n),totmixs(n)

        if(debug .and. s_n(n)>0.10_PS) then
           write(*,*) "update_group:k,i,j,p,qv,qr,qi,t,tn,sv1,svn1,sn,rv2",KD(n),ID(n),JD(n),ag%TV(n)%P,&
               ag%tv(n)%rv,totmixr(n),totmixs(n),ag%TV(n)%T,&
               ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n),rv2(n)
        endif
        if(ag%tv(n)%rv<0.0_PS.or.s_n(n)>0.50) then
          if(ag%tv(n)%rv<0.0_PS) then
             ! assume vapor mixing ratio to 10% of ice saturation.
             rv_n=(-0.9_PS+1.0_PS)*ag%tv(n)%e_sat(2)*Rdvchiarui/&
                      (ag%tv(n)%P-(-0.9+1.0_PS)*ag%tv(n)%e_sat(2))

             if ( debug ) then
                write(*,*) "rv is negative"
                write(*,*) "update_group neg", KD(n),ID(n),JD(n),ag%TV(n)%P,&
                     qtp(n),ag%tv(n)%rv,rv_n,totmixr(n),totmixs(n),ag%TV(n)%T,&
                     ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n)

                write(*,*) "mes_rc,density of air",mes_rc(n),ag%TV(n)%den
                write(*,*) "rain"
                write(*,*) (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
                do i=1,gr%N_BIN
                   write(*,*) i,(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%N_tendpros)
                end do

                write(*,*) "ice"
                write(*,*) (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)
                do i=1,gs%N_BIN
                   write(*,*) i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gs%N_tendpros)
                end do
                write(*,*) ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),&
                     ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
             end if

            ag%tv(n)%rv=rv_n
          else
             if ( debug ) then
                write(*,*) "rv(n) is super saturated"
                write(*,*) "update_group su", KD(n),ID(n),JD(n),ag%TV(n)%P,&
                     qtp(n),ag%tv(n)%rv,totmixr(n),totmixs(n),ag%TV(n)%T,&
                     ag%TV(n)%T_n,ag%TV(n)%s_v(1),ag%TV(n)%s_v_n(1),s_n(n)
                write(*,*) "mes_rc,density of air",mes_rc(n),ag%TV(n)%den
                write(*,*) "rain"
                write(*,*) (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
                do i=1,gr%N_BIN
                   write(*,*) i,(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%N_tendpros)
                end do

                write(*,*) "ice"
                write(*,*) (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)
                do i=1,gs%N_BIN
                   write(*,*) i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gs%N_tendpros)
                end do
                write(*,*) ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),&
                     ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
             end if
          endif
        endif
      enddo

    endif
```

### Aerosol block header — 3422-3436 (verbatim); body 3437-3769 is per-category interstitial-aerosol arithmetic; routine closes at 3770

```fortran
    !
    ! Aerosol
    !
    if(flagp_a/=0) then

      do ic=1,nacat
        ! +++ for aerosols +++
        if((flagp_a==-1.or.flagp_a==-4).and.ic==2) then
          cycle
        elseif((flagp_a==-2.or.flagp_a==-5).and.ic/=2) then
          cycle
        elseif(flagp_a==-3.or.flagp_a==-6) then
          !write(fid_alog,*) "flaga", MAXVAL(ga(1)%MS(:,:)%con), MAXLOC(ga(1)%MS(:,:)%con), MINVAL(ga(1)%MS(:,:)%con), MINLOC(ga(1)%MS(:,:)%con)
          cycle
        end if
        ...  ! [3438-3769: per-bin aerosol mass/con advance, analogous to rain block]
      end do      ! ic
    endif         ! flagp_a
  end subroutine update_group_all
```

---

## Orchestration summary for the M2 host loop

Per host call `ifc_cloud_micro`, the microphysics state advance is:

```
ini_cloud_micro → reality_check → check_water_apmass("make_cloud",0)
→ cal_micro_tendency:
     for it_cl = 1 .. n_step_cl        (dt = dt_cl; masks 889-893)
       [refresh: update_mesrc→diag_t→update_airgroup→diag_pq]   (only it_cl>1, plus level>=6 mv_ice2liq block at it_cl==1)
       coalescence(rain,rain)          [WARM]
       coalescence(ice,ice)            [ICE-ONLY]
       coalescence(ice,rain)           [ICE-ONLY riming]
       melting_shedding                [MIXED]
       hydrodyn_breakup(rain)          [WARM]
       hydrodyn_breakup(ice)           [ICE-ONLY]
       Ice_Nucleation1                 [ICE-ONLY]
       repair(af_col, .false.,.true.)
       update_group_all(update_vapor=0)
       cal_dmtend_scale
       for it_vp = 1 .. n_step_vp      (dt = dt_vp = dt_cl/n_step_vp; masks 1183-1191)
         [refresh: update_mesrc→diag_t→update_airgroup→diag_pq]
         cal_aptact_var8_vec | cal_aptact_var8_kc04dep   [WARM activation]
         vapor_deposition(rain)        [WARM]
         vapor_deposition(ice)         [ICE-ONLY]
         Ice_Nucleation2               [ICE-ONLY]
         repair(af_vap, .true.,.false.)
         update_group_all(update_vapor=1)   ← recomputes rv/rv2/s_n
         cal_dmtend_scale
     [final refresh: update_mesrc→diag_t→update_airgroup→diag_pq]
→ return_output → check_water_apmass("return_model",1) → update_terminal_vel
```

Warm-only reduction: keep `coalescence(rain,rain)`, `hydrodyn_breakup(rain)`, warm activation (`cal_aptact_var8_*`), `vapor_deposition(rain)`, and the full refresh/`repair`/`update_group_all`/`cal_dmtend_scale` skeleton; drop the four `[ICE-ONLY]` calls (ice–ice aggregation, ice–rain riming, `hydrodyn_breakup(ice)`, `Ice_Nucleation1/2`, `vapor_deposition(ice)`). `melting_shedding` and the `mv_ice2liq` (level>=6) block are ice-sourced and only fire when `flagp_s>0`.

Relevant files: `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Cloud_Micro.F90`, `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_lib.F90`, `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_core.F90`, `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Group.F90`.