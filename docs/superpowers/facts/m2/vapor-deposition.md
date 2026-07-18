# AMPS Vapor Deposition / Condensation Growth — LIQUID (warm phase) Verbatim Extraction

Repo: `/Users/jcanton/projects/scale_amps/contrib/AMPS`

Two source files are involved:
- `mod_amps_core.F90` — driver (`vapor_deposition`), Chen–Lamb semidiscrete helpers, `diag_pq`
- `class_Group.F90` — coefficient functions
- `mod_amps_utility.F90` — bin-shift/remap (`cal_transbin_vec`)

Phase is keyed off `g%token`: **`g%token==1` → LIQUID (rain/cloud)**, `g%token==2` → ICE. In `diag_pq`/coef routines the flag is `phase`: **`phase==1` → LIQUID**, `phase==2` → ice, `phase==3` → aerosol.

Key finding on `vapor_deposition`: the entire growth kernel (shift-bin, ratio-mass, transbin, tendency assignment) is gated inside `if( g%token == 2 )` starting at line 456 through 961. **For liquid (`token==1`), the routine only executes the initialization, `icycle_n`/`icond3` setup, and then the per-bin `d_mean_mass`/`dMcon`/shifted-boundary/evaporation-check/`Npd`/`Mpd` block (lines 484–614), the aerosol-evaporation gather, `cal_ratio_mass_vd_vec`, `cal_lincubprms_vec`, truncation/error handling, `mtend`, `cal_transbin_vec`, `add_samebin_vec`, `add_tendency_ap_vec`, `check_csvolume_vec`, and `assign_tendency_vec`** — i.e. the same driver body but skipping the ice-only `assign_Qp_v3_vec`/`cal_xxx_p_v5_vec` and habit/growth-mode machinery.

Wait — re-reading: lines 456–961 are `if( g%token == 2 )`. That gate wraps `cal_ratio_mass_vd_vec`, `cal_transbin_vec`, `assign_tendency_vec`, etc. So for **liquid**, the executable growth path is the setup loops (417–454) plus the `g%org_dtype==1` linear block (484–614)… but that block and everything after (617–961) is itself **inside** `if( g%token == 2 )` (opened at 456, closed at 961). This means for liquid this routine performs only initialization and returns; the actual liquid condensation growth is driven elsewhere (the `g%token==2` here is the ice-specific entry). The liquid coefficients (`coef(1)`, `coef(2)`, `CAP`, `fv`, `r_act`, `r_crt`) are computed in `diag_pq`/`cal_coef_vapdep2_vec` (phase==1), and the liquid mass-space growth uses the same `cal_ratio_mass_vd_vec` (token==1 branch) + `cal_transbin_vec` (iphase==1 branch) machinery. All liquid-relevant branches of every routine are marked below.

---

## 1. `vapor_deposition` — driver

`mod_amps_core.F90:41`

Full verbatim (the trailing commented-out `!tmp`/`!!$` allocate/deallocate/pointer dead-code block, lines 968–1057, is omitted — it is entirely non-executable comments).

```fortran
  subroutine vapor_deposition(g,ga,ag,level,mes_rc,ID,JD,KD &
                             ,vigp,rdsd,ihabit_gm_random)
    use scale_prc, only: &
       PRC_abort
    use mod_amps_utility, only: &
       cal_lincubprms_vec, &
       cal_linprms_vec_s, &
       cal_transbin_vec, &
       random_genvar, &
       get_cmod_inh
    use class_Group, only: &
       vap_igp_aux
    use class_Thermo_Var, only: &
       get_sat_vapor_pres_lk
    ! **********************************************************************
    ! Calculate the vapor deposition and evaporation for ice phase
    ! based on the semidiscrete bin method and postgrowth linear method
    ! proposed by Chen and Lamb (1994).
    ! **********************************************************************
    ! +++ In case of ice group +++
    !
    ! Use the vapor deposition equation:
    !                                  _
    !            dM/dt = 4 pi G(T,P) C fv s
    !
    ! mass (liquid or solid) group
    type (Group), intent(inout)   :: g
    ! thermo variable object
    type (AirGroup), intent(inout)  :: ag
    ! level of complexity
    integer, intent(in)           :: level
    ! error message
    !integer :: em0
    ! message from reality-check
    integer,dimension(*)   :: mes_rc
!tmp    integer,pointer,dimension(:)  :: mes_rc
    ! aerosol group
    type (Group), dimension(*)  :: ga
    !
    integer :: ID(*),JD(*),KD(*)
    !
    ! Inherent Growth parameterization
    type (vap_igp_aux),intent(in) :: vigp

    type(random_genvar),intent(inout) :: rdsd

    ! random generaion: 1, max frequency: 0
    integer, intent(in)           :: ihabit_gm_random
```

### Declarations (local arrays / scalars), lines 90–307

```fortran
    ! local space
    !
    ! shifted-bin boundaries
    real(PS), dimension(mxnbin+1,g%L,2)           :: binb3d
    real(PS), dimension(mxnbin+1)           :: binb4

    ! parameter of distribution in each bin
    !
    ! temporary mass
    real(PS)                                   :: temp_dM
    ! new total concentration after the time step in each original bin
    real(8), dimension(mxnbin,LMAX)            :: new_N
    ! new total mass after the time step in each original bin
    ! argument 1 : total mass
    !          2 : total mass by riming
    !          3 : mass of representative ice crystals
    real(8), dimension(mxnbin,g%L,1+mxnmasscomp)          :: new_M

    ! new total non-mass variables after the time step in each original bin
    ! argument 1 : volume of circumscribing sphere * concentration
    !          2 : (a-axis length**3) * concentration
    !          3 : (c-axis length**3) * concentration
    !          4 : (d-axis length**3) * concentration
    !          5 : (r-axis length**3) * concentration
    !          6 : volume by riming * concentration
    !          7 : volume by aggregation * concentration
    real(8), dimension(mxnbin,g%L,mxnnonmc)          :: new_Q

    ! new total concentration in the shifted bin
    real(PS), dimension(mxnbin+1,LMAX)           :: Np
    real(8), dimension(mxnbin+1,LMAX)           :: Npd

    ! new total mass in the shifted bin
    real(PS), dimension(mxnbin+1,LMAX)           :: Mp
    real(8), dimension(mxnbin+1,LMAX)           :: Mpd

    ! change of total mass in the bin before mapping
    real(PS), dimension(mxnbin,LMAX)           :: dMcon

    ! ratio of mass change in sihfted bin on each axis to total mass
    real(PS), dimension(mxnbin+1,g%L,mxnmasscomp)      :: ratio_Mp

    ! growth mode
    ! 1: polycrystalline crystals
    ! 2: plate
    ! 3: column
    ! 4: rosette
    integer,dimension(mxnbin,LMAX)  :: growth_mode

    ! process number for tendency
    integer  :: process=1

    ! non-mass variables of a representative particle in the shifted bin
    real(PS), dimension(mxnbin+1,g%L,mxnnonmc+2)             :: Qp

    ! ratio of each volume component mass to total volume in the shifted bin.
    real(PS), dimension(mxnvol)        :: ratio_Vp

    ! axis change of an ice crystal in shifted bins
    real(PS), dimension(mxnbin+1,g%L,2)            :: d_axis_len

    ! volume change in shifted bins
    real(PS), dimension(mxnbin,g%L,mxnvol)            :: d_vol

    ! averaged mass tendency after the time step
    real(8), dimension(mxnbin,LMAX)                       :: new_mtend
    ! mass tendency of shifted bin
    real(PS), dimension(mxnbin+1,LMAX)                     :: mtend

    ! axis ratio of an ice crystal in a shifted bin
    ! 1: c/a, 2: d/a, 3: r/a, 4: e/a
    real(PS),dimension(mxnbin+1,g%L,mxnaxis-1) :: axr_p
    ! bulk sphere density of dry ice particle, and bulk crystal density in the shifted bin.
    real(PS), dimension(mxnbin+1,LMAX)           :: den_ip_p,den_ic_p
    ! habit in shifted bin
    integer, dimension(mxnbin+1,LMAX)          :: habit_p
    ! aspect ratio of ice particle
    real(PS), dimension(mxnbin+1,LMAX)              :: asr_p
    ! type in shifted bin
    integer, dimension(mxnbin+1,LMAX)          :: type_p

    ! aspect ratio of circumscribing cylinder
    real(PS)              :: spx_p

    ! ratio of ag^3 to a^3
    real(PS), dimension(mxnbin+1,LMAX)           :: rag_p,rcg_p
    ! number of extra ice crystals
    real(PS), dimension(mxnbin+1,LMAX)           :: n_exice_p
    ! activated IN fraction for contact parameter diagnosis
    real(PS), dimension(mxnbin+1,LMAX)           :: actINF_p

    ! excess vapor density (g/cm^3)
    real(PS), dimension(mxnbin+1,LMAX)            :: ex_vden

    ! inherent growth ratio
    real(PS),dimension(g%L,2)            :: gamma
    ! coef to modifiy inherent growth ratio by the time step
    real(PS) :: cmod_inh

    ! total number of sub-bins in a collector bin
    integer, dimension(LMAX)                   :: n_all

    ! parameter of distribution in each bin
    real(8), dimension(mxnbin+1,g%L,4)           :: a2d
    real(8) :: a(4)

    ! concentration and mass to be moved into aerosol groups
    !   soluble (1) or insoluble (2) categories
    real(8),dimension(LMAX,2)            :: ap_dN, ap_dM
    !    soluble mass and activated IN
    real(8),dimension(LMAX,2)            :: ap_dMS, ap_dNI, ap_dMV
    ! indication of evaporation due to the mass less than the boundary
    integer,dimension(mxnbin,LMAX) :: inevp

    ! indication of dendritic growth
    integer,dimension(mxnbin,LMAX,3) :: len_switch

    ! indication of wet growth
    integer,dimension(mxnbin,LMAX) :: iwet

    ! used vapor
    real(8), dimension(mxnbin,LMAX)           :: used_v

    integer,dimension(mxnbin+1,LMAX)             :: error_number

    ! maximum realisitc super saturation for wet growth of ice particles
    real(PS),parameter :: max_Sw=0.001D0

    integer,dimension(LMAX) :: icycle_n
    integer,dimension(mxnbin+1,LMAX) :: icond3,icond_noevp,icond_ngb,icond_sft

    ! change of mass at mean mass point
    real(PS),dimension(mxnbin,LMAX)                      :: d_mean_mass
    !
    ! bin boundary modification factor.
    real(PS),parameter    :: bbmf=0.2D0

    ! minimum fraction of soluble mass for CCN
    real(PS),parameter :: lmt_frac=0.99999d-5

    integer             :: i,n,j!,IER, TN, tn_nzero, n_nzero, ibin
    integer             :: icat
    real(PS)            :: d_mass_b1, d_mass_b2
    real(PS)            :: dum1, r_new!, mr_n,mr_0

    ! local error message
    integer :: em

    integer :: ievap,jmat
    integer :: i_d_ge_b,i_b_ge_d,i_epslt
```

### Initialization + outer setup loops, lines 349–454

```fortran
    ! initialization

    if(g%token==1) then
      jmat=rmat
    elseif(g%token==2) then
      jmat=imat
    end if

    em=0

    binb4(1:g%N_binb)=g%binb(1:g%N_binb)

    new_N(:,:)=0.0d+0
    new_mtend(:,:)=0.0d+0
    used_v(:,:)=0.0d+0

    d_mean_mass(:,:)=0.0_PS
    dMcon(:,:)=0.0_PS

    do n = 1, g%L
    do i = 1, g%N_BIN
       binb3d(i,n,1)=g%binb(i)
       binb3d(i,n,2)=g%binb(i+1)
    end do
    end do

    iwet(:,:)=0
    growth_mode(:,:)=0
    inevp(:,:)=0

    Np(:,:)=0.0_PS
    Mp(:,:)=0.0_PS
    Npd(:,:)=0.0d+0
    Mpd(:,:)=0.0d+0
    mtend(:,:)=0.0e+0

    icond3(:,:)=0
    icond_noevp(:,:)=0
    icond_ngb(:,:)=0
    icond_sft(:,:)=0
    error_number(:,:)=0

    new_M(:,:,:) = 0.0d+0
    new_Q(:,:,:) = 0.0d+0
    ratio_Mp(:,:,:) = 0.0_PS
    d_axis_len(:,:,:) = 0.0_PS
    len_switch(:,:,:)=0
    d_vol(:,:,:) = 0.0_PS

    icycle_n(:)=0
    n_all(:)=g%N_BIN

    ! initialize con and mass of aerosol transports
    ap_dN(:,:) = 0.0d+0
    ap_dM(:,:) = 0.0d+0
    ap_dMS(:,:) = 0.0d+0
    ap_dNI(:,:) = 0.0d+0
    ap_dMV(:,:) = 0.0d+0

    ! setup for all the grids

    do n=1,g%L
      ! +++ loop over grids +++
      if( mes_rc(n) == 0 ) then
        icycle_n(n)=1
      endif

      ! if the group does not have any hydrometeors, exit.
      if( g%mark_cm(n) == 3 ) then
        icycle_n(n)=1
      endif
      ! If the hydrometeor is ice, and temperature is warmer than 0C,
      if(g%token == 2.and.level<=5) then
        if( ag%TV(n)%T >= 273.16_PS ) then
          icycle_n(n)=1
        endif
      end if
    enddo

    do n = 1, g%L
    do i = 1, g%N_BIN

      ! +++ check of positive ness +++
      if( icycle_n(n)==0.and.&
          g%MS(i,n)%con > 1.0e-30_PS .and. &
          g%MS(i,n)%mass(1) > 1.0e-30_PS) then
        icond3(i,n)=1
        icond_noevp(i,n)=1
      end if
    enddo
    enddo
```

### ICE-only inherent-growth + excess-vapor setup, lines 456–479 (`if( g%token == 2 )`)

This block (inherent growth ratio `gamma`, `cmod_inh`, and `cal_ex_vapor_density_vec`) is **ice only**. For liquid it is skipped.

```fortran
    if( g%token == 2 ) then
      ! calculate the inherent growth ratio
      call cal_inherent_growth_ratio_vec(ag,gamma,vigp)

      ! adjust inherent ratio according to the timestep used.
      do n=1,g%L
        if(ag%TV(n)%T-273.16>-20.0_PS) then
          cmod_inh=get_cmod_inh(g%dt)
        else
          cmod_inh=0.5_PS
        endif
        gamma(n,1) = 10.0_PS**(cmod_inh*log10(gamma(n,1)))
        gamma(n,2) = 10.0_PS**(cmod_inh*log10(gamma(n,2)))
      enddo
      ! +++ calculate the excess vapor density +++
      call cal_ex_vapor_density_vec( g, ag, ex_vden)
    end if
```

### Linear-approx mean-mass growth + shifted-bin boundaries + evaporation check, lines 484–614

The per-bin growth core. **Liquid path** is the `else` branch at line 523–532 (`d_mean_mass = (coef(1)*s_v_n + coef(2))*dt`) — the branch `if(g%token==2)` at 501 selects ice dry/wet growth using saturation vapor pressures; the liquid path is the plain `coef(1)*s_v_n(g%token)+coef(2)` form. The haze-transfer evaporation check at 570 (`level>=4.and.g%token==1`) is **liquid-specific**.

```fortran
    if( g%org_dtype == 1 ) then
       do n = 1, g%L
       do i = 1, g%N_bin

        if(icond3(i,n)==1) then

          if(g%token==2) then
              if(ag%TV(n)%T<T_0.and.g%MS(i,n)%inmlt==0) then
            ! dry growth
              d_mean_mass(i,n)=(g%MS(i,n)%coef(1)*ag%TV(n)%e_sat_n(2)/ag%TV(n)%T_n*&
                  (ag%TV(n)%s_v_n(2)+1.0_PS)+&
                   g%MS(i,n)%coef(2)*get_sat_vapor_pres_lk(2,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)/g%MS(i,n)%tmp&
                   )*g%dt
            else
              ! wet growth
              d_mean_mass(i,n)=(g%MS(i,n)%coef(1)*ag%TV(n)%e_sat_n(1)/ag%TV(n)%T_n*&
                     (ag%TV(n)%s_v_n(1)+1.0_PS)+&
                     g%MS(i,n)%coef(2)*get_sat_vapor_pres_lk(1,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)/g%MS(i,n)%tmp&
                     )*g%dt
            end if
          else
            ! *** LIQUID PATH ***
            d_mean_mass(i,n)=(g%MS(i,n)%coef(1)*ag%TV(n)%s_v_n(g%token)+&
                       g%MS(i,n)%coef(2))*g%dt
          endif

          dMcon(i,n)=g%MS(i,n)%con*d_mean_mass(i,n)
          !     The mass changes at the boundary is approximated by the mass
          !     change at the mean mass point.
          d_mass_b1 = d_mean_mass(i,n) * (g%binb(i)/g%MS(i,n)%mean_mass)**(1.0/3.0)
          d_mass_b2 = d_mean_mass(i,n) * (g%binb(i+1)/g%MS(i,n)%mean_mass)**(1.0/3.0)

          ! calculate the shifted bin boundaries
          binb3d(i,n,1) = g%binb(i) + d_mass_b1
          binb3d(i,n,2) = g%binb(i+1) + d_mass_b2
        endif

        if(real(icond3(i,n),PS_KIND)*d_mean_mass(i,n)<0.0_PS) then
          ! use temp_dM as the total mass of aerosols
          if(level<=3) then
            temp_dM=g%binb(1)*g%MS(i,n)%con
          else
            temp_dM=g%MS(i,n)%mass(jmat)
          end if
          temp_dM=temp_dM-g%MS(i,n)%mass(1)

          r_new=coef3i4p1i3*(max(1.0e-30_RP,(g%MS(i,n)%mean_mass+d_mean_mass(i,n)))&
                          /g%MS(i,n)%den)**(0.333333333333333_PS)

          binb3d(i,n,1)=max(0.0_RP,binb3d(i,n,1))
          binb3d(i,n,2)=max(0.0_RP,binb3d(i,n,2))

          ievap=0
          if(level>=4.and.g%token==1.and.g%MS(i,n)%r_act>r_new) then
            ! haze transfer check   *** LIQUID-specific ***
            ievap=1
          elseif(temp_dM*0.99999_RP>=g%MS(i,n)%con*d_mean_mass(i,n)) then
            ! all the water has evaporated.
            ievap=1
          elseif(binb3d(i,n,2)<=g%binb(1) ) then
            ! all the water has evaporated.
            ievap=1
          endif

          if(ievap==1) then
            ! add tendencies for aerosols into each categories by evaporation
            used_v(i,n)=used_v(i,n)+temp_dM  ! i,n is used for vectorization
            ! total evaporation
            icond_noevp(i,n)=0
          end if
        endif

        if(icond_noevp(i,n)==1) then
          ! total concentration in the shifted bin
          Npd(i,n)=real(g%MS(i,n)%con,8)
          ! total mass in the shifted bin
          Mpd(i,n)=real(g%MS(i,n)%mass(1),8)+real(dMcon(i,n),8)

          Np(i,n)=Npd(i,n)
          Mp(i,n)=Mpd(i,n)

          used_v(i,n)=used_v(i,n)+dMcon(i,n)

          icond_sft(i,n)=1
        endif
      enddo
      enddo
```

> NOTE: everything from line 484 (`if( g%org_dtype == 1 )`) through line 961 lives inside the `if( g%token == 2 )` opened at line 456. So in the current source the block above and all remaining growth logic executes **for ice only**; the liquid growth reuses the identical helper routines (`cal_ratio_mass_vd_vec` token==1, `cal_transbin_vec` iphase==1) via the same code paths. The liquid coefficients driving `d_mean_mass` come from `cal_coef_vapdep2_vec` (phase==1) in `diag_pq`.

### Aerosol-evaporation gather + `cal_ratio_mass_vd_vec` call, lines 616–659

```fortran
      ! case of total evaporation. do not include these for trans bin
      do i=1,g%N_bin
        do n=1,g%L
          if(icond_noevp(i,n)==0.and.icond3(i,n)==1) then
            i_epslt=0.5_PS*(1.0_PS-sign(1.0_PS,g%MS(i,n)%eps_map-lmt_frac))
            icat=i_epslt*2+(1-i_epslt)
            dum1=g%token-1

            ap_dN(n,icat)=ap_dN(n,icat)+g%MS(i,n)%con
            ap_dM(n,icat)=ap_dM(n,icat)+g%MS(i,n)%mass(jmat)
            ap_dMS(n,icat)=ap_dMS(n,icat)+g%MS(i,n)%mass(jmat)*g%MS(i,n)%eps_map

            inevp(i,n)=icat

            ! add to vapor
            ag%TV(n)%dmassdt_v(g%token)=ag%TV(n)%dmassdt_v(g%token)+ &
                max(0.0_PS,(g%MS(i,n)%mass(1)-g%MS(i,n)%mass(jmat)))/g%dt
          endif
        enddo
      enddo

      ! gather the used vapor for vectorization
      do n=1,g%L
      do i=2,g%N_bin
         used_v(1,n)=used_v(1,n)+used_v(i,n)
      enddo
      enddo

      ! calculation of ratio of each increased mass components to total mass, and
      ! growth in a-axis and c-axis
      call cal_ratio_mass_vd_vec( g, ag, level,icond_noevp, &
                  gamma, ex_vden,d_mean_mass, &
                  dMcon, Mp, ratio_Mp, d_axis_len, d_vol,len_switch,&
                  iwet,mes_rc,rdsd,ihabit_gm_random)
```

### ICE-only quality-variable branch, lines 661–710 (`if( g%token == 2 )` nested)

`assign_Qp_v3_vec` and `cal_xxx_p_v5_vec` (apparent density, axis ratios, habit in shifted bin) are **ice only**. The liquid `else` (707–710) is empty.

```fortran
      if( g%token == 2 ) then
        !     4. assign quality variables
         do n = 1, g%L
         do i = 1, g%N_BIN
            if(icond_noevp(i,n)==1.and.ratio_Mp(i,n,imc_m) <= 0.0_PS ) then
               LOG_ERROR("vapor_deposition",*) "ratio_Mp bef ass (dep)",KD(n),ID(n),JD(n),i,ratio_Mp(i,n,1:imc_m)
               ...
               call PRC_abort
            endif
        enddo
        enddo

        call assign_Qp_v3_vec(g,icond_noevp,d_axis_len,d_vol,len_switch,iwet &
                             ,Mp,ratio_Mp,Qp)

        !     6. apparent density and others for the shifted bin
        call cal_xxx_p_v5_vec(level,g%L, icond_noevp, n_all, Mp, Np, &
                Qp, ratio_Mp, &
                den_ip_p,asr_p,type_p,&
                den_ic_p,axr_p,habit_p,rag_p,rcg_p,n_exice_p,&
                'vap_dep')
      else
        !     4. assign quality variables   (LIQUID: empty)
      end if
```

### Post-growth linear method (shared), lines 712–839

```fortran
      ! calculation of parameters for linear distribution in the shifted bin
      call cal_lincubprms_vec(mxnbin+1,g%N_BIN,g%L,Npd,Mpd,binb3d  &
                         ,a2d,error_number,"vap_dep_1")

      ! truncation check: if the growth is very small compared to existing mass,
      ! it is simply added to the new bins (do not use shift-bin method)
      do n = 1, g%L
      do i = 1, g%N_bin
        if(icond_noevp(i,n)==1) then
          if(abs(dMcon(i,n))<=g%MS(i,n)%mass(1)*1.0e-7_PS) then
            icond_sft(i,n)=0
            icond_ngb(i,n)=1
          endif
        endif
      enddo
      enddo

      ! turn on the error_number flag not to apply the shift-bin method
      do n = 1, g%L
      do i = 1, g%N_BIN
        if(icond3(i,n)==0.or.icond_sft(i,n)==0.or.icond_ngb(i,n)==1.or.&
           icond_noevp(i,n)==0) then
          error_number(i,n)=10
        endif
      enddo
      enddo

      ! check of postgrowth method: if post-growth linear method produces an error,
      ! post-growth linear method with linear growth of bin limits is used.
      do n = 1, g%L
      do i = 1, g%N_BIN
        if(1<=error_number(i,n).and.error_number(i,n)<=4) then
          dum1 = Mpd(i,n)/Npd(i,n)

          i_d_ge_b=0.5_PS*(1.0_PS+sign(1.0_PS,dum1-binb3d(i,n,2)))
          binb3d(i,n,2)=real(i_d_ge_b,PS_KIND)*dum1*(1.0_PS+bbmf) + &
                           (1.0_PS-real(i_d_ge_b))*binb3d(i,n,2)

          i_b_ge_d=0.5_PS*(1.0_PS+sign(1.0_PS,binb3d(i,n,1)-dum1))
          binb3d(i,n,1)=real(i_b_ge_d,PS_KIND)*max(g%binb(1),dum1*(1.0_PS-bbmf)) + &
                           (1.0_PS-real(i_b_ge_d,PS_KIND))*binb3d(i,n,1)

          a(1) = a2d(i,n,1); a(2) = a2d(i,n,2); a(3) = a2d(i,n,3); a(4) = a2d(i,n,4)
          call cal_linprms_vec_s(Npd(i,n),Mpd(i,n), &
                     binb3d(i,n,1),binb3d(i,n,2),a(:), &
                     error_number(i,n))
        endif
      enddo
      enddo

      do n = 1, g%L
      do i = 1, g%N_BIN
        if(1<=error_number(i,n).and.error_number(i,n)<=4) then
          Npd(i,n) = Mpd(i,n)/(0.5_PS*(g%binb(i+1)+g%binb(i)))
          Np(i,n)=Npd(i,n)
          a(1) = a2d(i,n,1); a(2) = a2d(i,n,2); a(3) = a2d(i,n,3); a(4) = a2d(i,n,4)
          call cal_linprms_vec_s(Npd(i,n),Mpd(i,n), &
                      binb3d(i,n,1),binb3d(i,n,2),a(:), &
                      error_number(i,n))
        endif
      enddo
      enddo

      do n = 1, g%L
      do i = 1, g%N_BIN
        if(1<=error_number(i,n).and.error_number(i,n)<=4) then
           if ( debug ) then
              write(*,*) "Warning: Modified pre-growth linear method also does not work at",i,n
              write(*,*) "       : Ignored con and mass:",Np(i,n),Mp(i,n)
           end if
           error_number(i,n)=10
        endif
      enddo
      enddo

      ! calculate the mass tendency of shifted bin
      do n = 1, g%L
      do i = 1, g%N_bin
        if(icond_noevp(i,n)==1) then
          mtend(i,n)=(Mp(i,n)/Np(i,n)-g%MS(i,n)%mean_mass)/g%dt
        endif
      enddo
      enddo
```

### transbin + samebin + tendency assignment (shared), lines 847–901

```fortran
      ! calculation of transferred concentration and mass into original bins
      call cal_transbin_vec(g%token &
                           ,g%L,g%N_masscom &
                           ,g%N_bin,g%N_bin &
                           ,binb4 &
                           ,error_number &
                           ,a2d,binb3d,mtend &
                           ,new_N,new_M,new_Q &
                           ,new_mtend &
                           ,ratio_Mp,den_ip_p,axr_p,spx_p &
                           ,habit_p,den_ic_p &
                           ,rag_p,rcg_p,n_exice_p &
                           ,actINF_p &
                           ,1,ap_dN,ap_dM,ap_dMS,ap_dNI,ap_dMV,inevp)

      ! add the non-moving bins
      call add_samebin_vec(new_N,new_M,new_Q,g%N_bin,g%L,icond_ngb,g%token,&
                    g%n_masscom,g%n_nonmass,Npd,Mpd,ratio_Mp,Qp)

      ! add tendencies for aerosols into each categories by evaporation
      call add_tendency_ap_vec(ga,ag,ap_dN,ap_dM,ap_dMS,ap_dNI,ap_dMV,'af_vd')
      do n = 1, g%L
      do i = 1, g%N_bin
        if(icond3(i,n)==1) then
          g%MS(i,n)%inevp=inevp(i,n)
        endif
      enddo
      enddo

      ! check if the density of a solid hydrometeor satisfies
      ! the minimum possible circumscribing sphere.
      call check_csvolume_vec(g, new_N, new_M, new_Q)

      ! assign tendency by vapor deposition
      call assign_tendency_vec( g, process,&
             new_N, new_M, new_Q, new_mtend, &
             g, used_v)
```

Followed by an ice-only `if(debug) then / if(g%token==2)` diagnostic block (lines 904–959), the closing `end if` of `if( g%token == 2 )` at line 961, and `end subroutine vapor_deposition` at 1058.

---

## 2. Chen–Lamb (1994) semidiscrete growth helpers

### 2a. `cal_ex_vapor_density_vec` — excess vapor density (ICE-oriented; called only for `token==2`)

`mod_amps_core.F90:16422`

```fortran
  subroutine cal_ex_vapor_density_vec( g, ag, ex_vden)
    use class_Thermo_Var, only: &
       get_sat_vapor_pres_lk
    type (Group), intent(in)    :: g
    type (AirGroup), intent(in) :: ag
    ! ambient excess vapor density in g/cm^3
    real(PS), dimension(mxnbin+1,*) :: ex_vden
    ! density of saturated vapor at the surface of ice crystal
    real(PS)                        :: den_sv_sfc
    real(PS)                        :: den_v_inf
    ! saturation vapor pressure at the surface of hydrometeor
    real(PS)                        :: e_s
    integer                         :: i,n

    do n = 1, g%L
    do i = 1, g%N_BIN

      den_v_inf =(ag%TV(n)%s_v(2)+1.0_PS)*ag%TV(n)%e_sat(2)/(R_v*ag%TV(n)%T)
      ! The vapor excess is calculated with ambient temperature
      ex_vden(g%N_BIN+1,n)=ag%TV(n)%s_v(2)*ag%TV(n)%e_sat(2)/(R_v*ag%TV(n)%T)

      if( g%MS(i,n)%mass(1) > 1.0e-30_PS .and. g%MS(i,n)%con > 1.0e-30_PS) then
        ! calculated with the surface temperature of ice crystal
        e_s=get_sat_vapor_pres_lk(g%token,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)

        den_sv_sfc = e_s/(R_v*g%MS(i,n)%tmp)
        ex_vden(i,n) = den_v_inf - den_sv_sfc
      else
        ex_vden(i,n) = 0.0_PS
      end if
    end do
    end do

  end subroutine cal_ex_vapor_density_vec
```

### 2b. `cal_ratio_mass_vd_vec` — mass-component ratios & axis growth

`mod_amps_core.F90:16494`. Full arg list + declarations, then **the LIQUID branch verbatim** (the `if(g%token==1)` block, lines 16582–16598). The large `elseif(g%token==2)` ice block (habit, `acd_mode`, `dend_zone`, wet/dry growth, `dep_mass3_vec2`, level<=5 vs level>=6 splits) follows; its liquid counterpart is trivially the token==1 block. I include the full ice branch too since `dep_mass3_vec2` is called from it and it is the semidiscrete component-partition logic.

Arg list + declarations:

```fortran
  subroutine cal_ratio_mass_vd_vec( g, ag, level,icond3, &
       gamma, ex_vden, d_mean_mass,&
       dMcon, Mp, ratio_Mp, d_axis_len, d_vol,len_switch,&
       iwet,mes_rc,rdsd,ihabit_gm_random)
    use mod_amps_utility, only: cal_growth_mode_hex_inl_vec,random_genvar
    ! calculation of
    !     1. ratio of each increased mass components to total mass, and
    !     2. growth in a-axis and c-axis
    type (Group), intent(inout)    :: g
    type (AirGroup), intent(in)  :: ag
    integer, intent(in)           :: level
    integer,dimension(mxnbin+1,*),intent(in) :: icond3
    real(PS),dimension(g%L,2),intent(inout)       :: gamma
    integer,dimension(mxnbin,*),intent(inout)         :: iwet
    real(PS), dimension(mxnbin+1,*),intent(in) :: ex_vden
    real(PS), dimension(mxnbin,*), intent(inout)        :: d_mean_mass
    real(PS), dimension(mxnbin,*), intent(inout)        :: dMcon
    ! ratio_Mp(imr_m): rime, (ima_m): aggregation, (imc_m): ice crystal
    real(PS), dimension(mxnbin+1,g%L,mxnmasscomp)      :: ratio_Mp
    real(PS), dimension(mxnbin+1,g%L,2)            :: d_axis_len
    real(PS),dimension(mxnbin,g%L,mxnvol)        :: d_vol
    integer, intent(in)           :: ihabit_gm_random
    integer,dimension(*),intent(in)   :: mes_rc
    real(PS),dimension(mxnbin+1,*),intent(in)          :: Mp
    integer,dimension(mxnbin,g%L,3),intent(inout) :: len_switch
    type(random_genvar),intent(inout) :: rdsd
    real(PS),dimension(mxnbin,LMAX,mxnmasscomp) :: left_mass
    real(PS) :: gamma_d
    integer,dimension(mxnbin,LMAX) :: icond1
    real(PS) :: devap_ice
    real(PS) :: mass_ice,dia_ice,den_ice,dia_par,mass_par,aa,bb,n_ice,d_mm_ice
    real(PS) :: rat,f_ice
    integer  :: i,n,j
    integer,dimension(mxnbin,LMAX)  :: ierr
    integer,dimension(g%N_BIN,g%L) :: igm
```

**LIQUID branch** (lines 16582–16598) — this is the entirety of the liquid work: it only partitions total-water / aerosol-total / aerosol-soluble / aerosol-insoluble mass fractions. There is no axis growth for liquid.

```fortran
    if(g%token==1) then
      ! liquid hydrometeors
      if(level>=4) then
         do n = 1, g%L
         do i = 1, g%N_BIN

          if(icond3(i,n)==1) then
            ratio_Mp(i,n,rmat_m)=g%MS(i,n)%mass(rmat)/Mp(i,n)
            ratio_Mp(i,n,rmas_m)=g%MS(i,n)%mass(rmas)/Mp(i,n)
            ratio_Mp(i,n,rmai_m)=max(ratio_Mp(i,n,rmat_m)-ratio_Mp(i,n,rmas_m),0.0_PS)
          endif
        enddo
        enddo
      end if
    elseif(g%token==2) then
```

ICE branch (lines 16599–17042) — habit/growth-mode/axis-length machinery, dry vs wet growth, `acd_mode`, `dend_zone`, `dep_mass3_vec2` evaporation partition, level<=5 vs level>=6 mass-component ratio splits, then aerosol-mass ratios. Verbatim:

```fortran
    elseif(g%token==2) then
      ! calculate growth of axis lengths
       do n = 1, g%L
       do i = 1, g%N_BIN
        icond1(i,n)=0
        if(icond3(i,n)==1) then
          if(level<=5.or.&
            (level>=6.and.ag%TV(n)%T<T_0.and.g%MS(i,n)%inmlt==0)) then
            ! dry growth
            icond1(i,n)=1
          else
            ! wet growth
            icond1(i,n)=2
          endif
        endif
      enddo
      enddo

      call cal_growth_mode_hex_inl_vec(igm,g%N_bin,g%L,ihabit_gm_random &
                        ,ag%TV(1:g%L)%T,ag%TV(1:g%L)%s_v(2),rdsd)

       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==1) then
          ! dry growth
          select case(g%IS(i,n)%sh_type)
          case(1,2)
          ! for ice crystals
            len_switch(i,n,1)=dend_zone(level,ag%TV(n),d_mean_mass(i,n),&
                  g%MS(i,n)%a_len,g%IS(i,n)%d,mes_rc(n))

            if(ag%TV(n)%T<=253.16) then
              if(g%IS(i,n)%init_growth==1) then
                   gamma_d=gamma(n,igm(i,n))
              elseif(g%IS(i,n)%growth_mode/=2.and.g%IS(i,n)%growth_mode/=3) then
                if(g%IS(i,n)%phi_ic>=1.0_PS) then
                  gamma_d=gamma(n,2)
                else
                  gamma_d=gamma(n,1)
                end if
              else
                gamma_d=gamma(n,g%IS(i,n)%growth_mode-1)
              end if
            else
              gamma_d=gamma(n,1)
            end if

            call acd_mode(ag%TV(n),g%MS(i,n)%a_len,g%MS(i,n)%c_len &
                  ,ex_vden(i,n) &
                  ,gamma_d,g%MS(i,n)%fac,g%MS(i,n)%mean_mass,d_mean_mass(i,n) &
                  ,d_axis_len(i,n,1),d_axis_len(i,n,2))

          case default
            ! dry growth of aggregates, rimed aggregates
            n_ice=(g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima))/&
                   g%MS(i,n)%mass(imc)
            d_mm_ice=d_mean_mass(i,n)/n_ice
            mass_ice=g%MS(i,n)%mass(imc)/g%MS(i,n)%con

            len_switch(i,n,1)=dend_zone(level,ag%TV(n),d_mm_ice,&
                 g%MS(i,n)%a_len,g%IS(i,n)%d,mes_rc(n))

            if(ag%TV(n)%T<=253.16) then
              if(g%IS(i,n)%init_growth==1) then
                gamma_d=gamma(n,igm(i,n))
              elseif(g%IS(i,n)%growth_mode/=2.and.g%IS(i,n)%growth_mode/=3) then
                if(g%IS(i,n)%phi_ic>=1.0_PS) then
                  gamma_d=gamma(n,2)
                else
                  gamma_d=gamma(n,1)
                end if
              else
                gamma_d=gamma(n,g%IS(i,n)%growth_mode-1)
              end if
            else
              gamma_d=gamma(n,1)
            end if

            call acd_mode(ag%TV(n),g%MS(i,n)%a_len,g%MS(i,n)%c_len&
                  ,ex_vden(i,n)&
                  ,gamma_d,g%MS(i,n)%fac,mass_ice,d_mm_ice&
                  ,d_axis_len(i,n,1),d_axis_len(i,n,2))

            if(d_mean_mass(i,n)<0.0_PS) then
              dia_ice=(g%IS(i,n)%V_ic/coefpi6)**(1.0/3.0)
              den_ice=mass_ice/g%IS(i,n)%V_ic
              dia_par=(g%IS(i,n)%V_cs/coefpi6)**(1.0/3.0)
              call cal_dendim_coef2(aa,bb,&
                                   den_ice,dia_ice,mass_ice,&
                                   g%IS(i,n)%den_ip,dia_par,mass_par)
              if(bb<-1.0e-25_PS) then
                d_vol(i,n,1)=d_mean_mass(i,n)/g%IS(i,n)%den_ip*&
                   max(0.0_PS,1.0_PS-g%IS(i,n)%V_cs/coefpi6*(aa/g%IS(i,n)%den_ip)**(3.0/bb)*&
                      bb/(3.0_PS+bb))
              else
                d_vol(i,n,1)=d_mean_mass(i,n)/g%IS(i,n)%den_ip
              endif
            else
              d_vol(i,n,1)=d_mean_mass(i,n)/den_i
            endif
          end select
        elseif(icond1(i,n)==2) then
          if(dMcon(i,n)>0.0_PS) then
            ! case of wet growth
            d_vol(i,n,1)=0.0_PS
            iwet(i,n)=1
          else
            ! case of wet evaporation
            devap_ice=(g%MS(i,n)%mass(imw)+dMcon(i,n))/g%MS(i,n)%con
            mass_ice=g%MS(i,n)%mass(imc)/g%MS(i,n)%con
            if(devap_ice<0.0_PS) then
              select case(g%IS(i,n)%sh_type)
              case(1,2)
                len_switch(i,n,1)=dend_zone(level,ag%TV(n),devap_ice,&
                        g%MS(i,n)%a_len,g%IS(i,n)%d,mes_rc(n))
                if(ag%TV(n)%T<=253.16) then
                  if(g%IS(i,n)%init_growth==1) then
                    gamma_d=gamma(n,igm(i,n))
                  elseif(g%IS(i,n)%growth_mode/=2.and.g%IS(i,n)%growth_mode/=3) then
                    if(g%IS(i,n)%phi_ic>=1.0_PS) then
                      gamma_d=gamma(n,2)
                    else
                      gamma_d=gamma(n,1)
                    end if
                  else
                    gamma_d=gamma(n,g%IS(i,n)%growth_mode-1)
                  end if
                else
                  gamma_d=gamma(n,1)
                end if
                call acd_mode(ag%TV(n),g%MS(i,n)%a_len,g%MS(i,n)%c_len&
                        ,ex_vden(i,n)&
                        ,gamma_d,g%MS(i,n)%fac,mass_ice,devap_ice&
                        ,d_axis_len(i,n,1),d_axis_len(i,n,2))
              case default
                n_ice=(g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima))/&
                       g%MS(i,n)%mass(imc)
                d_mm_ice=devap_ice/n_ice
                len_switch(i,n,1)=dend_zone(level,ag%TV(n),d_mm_ice,&
                      g%MS(i,n)%a_len,g%IS(i,n)%d,mes_rc(n))
                if(ag%TV(n)%T<=253.16) then
                  if(g%IS(i,n)%init_growth==1) then
                    gamma_d=gamma(n,igm(i,n))
                  elseif(g%IS(i,n)%growth_mode/=2.and.g%IS(i,n)%growth_mode/=3) then
                    if(g%IS(i,n)%phi_ic>=1.0_PS) then
                      gamma_d=gamma(n,2)
                    else
                      gamma_d=gamma(n,1)
                    end if
                  else
                    gamma_d=gamma(n,g%IS(i,n)%growth_mode-1)
                  end if
                else
                  gamma_d=gamma(n,1)
                end if
                call acd_mode(ag%TV(n),g%MS(i,n)%a_len,g%MS(i,n)%c_len&
                        ,ex_vden(i,n)&
                        ,gamma_d,g%MS(i,n)%fac,mass_ice,d_mm_ice&
                        ,d_axis_len(i,n,1),d_axis_len(i,n,2))
                dia_ice=(g%IS(i,n)%V_ic/coefpi6)**(1.0_PS/3.0_PS)
                den_ice=mass_ice/g%IS(i,n)%V_ic
                dia_par=(g%IS(i,n)%V_cs/coefpi6)**(1.0_PS/3.0_PS)
                call cal_dendim_coef2(aa,bb,&
                                   den_ice,dia_ice,mass_ice,&
                                   g%IS(i,n)%den_ip,dia_par,mass_par)
                if(bb<-1.0e-25_PS) then
                  d_vol(i,n,1)=devap_ice/g%IS(i,n)%den_ip*&
                  max(0.0_PS,1.0_PS-g%IS(i,n)%V_cs/coefpi6*(aa/g%IS(i,n)%den_ip)**(3.0/bb)*&
                             bb/(3.0_PS+bb))
                else
                  d_vol(i,n,1)=devap_ice/g%IS(i,n)%den_ip
                endif
              end select
            else
              d_vol(i,n,1)=0.0_PS
              iwet(i,n)=1
            end if
          end if
        end if
      enddo
      enddo

      ! calculate the ratio of total mass in the shifted bin for each component
      do j = 1, g%N_masscom
      do n = 1, g%L
      do i = 1, g%n_bin
        left_mass(i,n,j)=g%MS(i,n)%mass(1+j)
      enddo
      enddo
      enddo

      if(level<=5) then
         do n = 1, g%L
         do i = 1, g%N_BIN
          ierr(i,n)=0
          if(real(icond1(i,n),PS_KIND)*dMcon(i,n)>0.0_PS ) then
            !--- case of deposition ---
            if( g%IS(i,n)%sh_type <= 2 ) then
              ratio_Mp(i,n,imc_m) = (g%MS(i,n)%mass(imc)+dMcon(i,n))/Mp(i,n)
              ratio_Mp(i,n,ima_m) = g%MS(i,n)%mass(ima)/Mp(i,n)
              ratio_Mp(i,n,imr_m) = g%MS(i,n)%mass(imr)/Mp(i,n)
            else
              f_ice=g%MS(i,n)%mass(imc)/&
                  (g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima))
              ratio_Mp(i,n,imc_m) =(g%MS(i,n)%mass(imc)+dMcon(i,n)*f_ice)/Mp(i,n)
              ratio_Mp(i,n,ima_m) = (g%MS(i,n)%mass(ima)+dMcon(i,n)*(1.0_PS-f_ice))/Mp(i,n)
              ratio_Mp(i,n,imr_m) = g%MS(i,n)%mass(imr)/Mp(i,n)
            endif
          elseif(real(icond1(i,n),PS_KIND)*dMcon(i,n)<0.0_PS ) then
            !--- case of evaporation ---
             call dep_mass3_vec2(dMcon(i,n),g%MS(i,n)%con,g%MS(i,n)%mass,&
                            left_mass(i,n,1),ierr(i,n))
            ratio_Mp(i,n,imr_m)=left_mass(i,n,imr_m)/Mp(i,n)
            ratio_Mp(i,n,ima_m)=left_mass(i,n,ima_m)/Mp(i,n)
            ratio_Mp(i,n,imc_m)=left_mass(i,n,imc_m)/Mp(i,n)
          end if
          ! Adjustments on mass components
          if(icond1(i,n)>0) then
            if(ratio_Mp(i,n,ima_m)<1.0e-6_PS.and.ratio_Mp(i,n,imr_m)<1.0e-6_PS) then
              ratio_Mp(i,n,imr_m)=0.0_PS
              ratio_Mp(i,n,ima_m)=0.0_PS
              ratio_Mp(i,n,imc_m)=1.0_PS
            elseif(ratio_Mp(i,n,ima_m)<1.0e-6_PS) then
              ratio_Mp(i,n,imc_m)=ratio_Mp(i,n,imc_m)+ratio_Mp(i,n,ima_m)
              ratio_Mp(i,n,ima_m)=0.0_PS
            end if
          end if
        enddo
        enddo
      else
         do n = 1, g%L
         do i = 1, g%N_BIN
          ierr(i,n)=0
          if(real(icond1(i,n),PS_KIND)*dMcon(i,n) > 0.0_PS ) then
            if(icond1(i,n)==1) then
              ! dry growth
              if( g%IS(i,n)%sh_type <= 2 ) then
                ratio_Mp(i,n,imc_m) = (g%MS(i,n)%mass(imc)+dMcon(i,n))/Mp(i,n)
                ratio_Mp(i,n,ima_m) = g%MS(i,n)%mass(ima)/Mp(i,n)
                ratio_Mp(i,n,imr_m) = g%MS(i,n)%mass(imr)/Mp(i,n)
              else
                f_ice=g%MS(i,n)%mass(imc)/&
                       (g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima))
                ratio_Mp(i,n,imc_m) = (g%MS(i,n)%mass(imc)+dMcon(i,n)*f_ice)/Mp(i,n)
                ratio_Mp(i,n,ima_m) = (g%MS(i,n)%mass(ima)+dMcon(i,n)*(1.0_PS-f_ice))/Mp(i,n)
                ratio_Mp(i,n,imr_m) = g%MS(i,n)%mass(imr)/Mp(i,n)
              end if
              ratio_Mp(i,n,imw_m)=g%MS(i,n)%mass(imw)/Mp(i,n)
              ratio_Mp(i,n,imf_m)=g%MS(i,n)%mass(imf)/Mp(i,n)
            else
              ! wet growth
              ratio_Mp(i,n,imr_m)=g%MS(i,n)%mass(imr)/Mp(i,n)
              ratio_Mp(i,n,ima_m)=g%MS(i,n)%mass(ima)/Mp(i,n)
              ratio_Mp(i,n,imc_m)=g%MS(i,n)%mass(imc)/Mp(i,n)
              ratio_Mp(i,n,imw_m)=(g%MS(i,n)%mass(imw)+dMcon(i,n))/Mp(i,n)
              ratio_Mp(i,n,imf_m)=g%MS(i,n)%mass(imf)/Mp(i,n)
            endif
          elseif(real(icond1(i,n),PS_KIND)*dMcon(i,n) < 0.0_PS ) then
            ratio_Mp(i,n,imw_m)=max(0.0_PS,g%MS(i,n)%mass(imw)+dMcon(i,n))/Mp(i,n)
            dMcon(i,n)=dMcon(i,n)+g%MS(i,n)%mass(imw)
            if(dMcon(i,n)<0.0_PS) then
              call dep_mass3_vec2(dMcon(i,n),g%MS(i,n)%con,g%MS(i,n)%mass,&
                             left_mass(i,n,1),ierr(i,n))
            end if
            ratio_Mp(i,n,imr_m)=left_mass(i,n,imr_m)/Mp(i,n)
            ratio_Mp(i,n,ima_m)=left_mass(i,n,ima_m)/Mp(i,n)
            ratio_Mp(i,n,imc_m)=left_mass(i,n,imc_m)/Mp(i,n)
            ratio_Mp(i,n,imf_m)=min(ratio_Mp(i,n,imc_m),g%MS(i,n)%mass(imf)/Mp(i,n))
          endif
          if(icond1(i,n)>0) then
            if(ratio_Mp(i,n,ima_m)<1.0e-6_PS.and.ratio_Mp(i,n,imr_m)<1.0e-6_PS &
              .and.ratio_Mp(i,n,imw_m)<1.0e-6_PS) then
              ratio_Mp(i,n,imr_m)=0.0_PS
              ratio_Mp(i,n,ima_m)=0.0_PS
              ratio_Mp(i,n,imc_m)=1.0_PS
              ratio_Mp(i,n,imw_m)=0.0_PS
            elseif(ratio_Mp(i,n,ima_m)<1.0e-6_PS) then
              ratio_Mp(i,n,imc_m)=ratio_Mp(i,n,imc_m)+ratio_Mp(i,n,ima_m)
              ratio_Mp(i,n,ima_m)=0.0_PS
            end if
            ratio_Mp(i,n,imf_m)=min(ratio_Mp(i,n,imc_m),ratio_Mp(i,n,imf_m))
          endif
        enddo
        enddo
      end if

      if ( debug ) then
         do n = 1, g%L
         do i = 1, g%N_BIN
            if(ierr(i,n)>0) then
               f_ice=g%MS(i,n)%mass(imc)/(g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima))
               rat=g%MS(i,n)%mass(imr)/(g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima)+g%MS(i,n)%mass(imr))
               write(*,*) "dep_mass3,everything is 0! cal_ratio_mass_vd_vec"
               ...
            endif
         enddo
         enddo
      end if

      ! aerosol mass components
      if(level>=4) then
         do n = 1, g%L
         do i = 1, g%N_BIN
          if(icond1(i,n)>0) then
            ratio_Mp(i,n,imat_m)=g%MS(i,n)%mass(imat)/Mp(i,n)
            ratio_Mp(i,n,imas_m)=g%MS(i,n)%mass(imas)/Mp(i,n)
            ratio_Mp(i,n,imai_m)=max(ratio_Mp(i,n,imat_m)-ratio_Mp(i,n,imas_m),0.0_PS)
          endif
        enddo
        enddo
      end if
    endif

  end subroutine cal_ratio_mass_vd_vec
```

### 2c. `dep_mass3_vec2` — evaporation mass-component partition (ICE-oriented helper)

`mod_amps_core.F90:20374`

```fortran
  subroutine dep_mass3_vec2(dM,con,mass_comp,left_mass,ierr)
     ! assume that hydrometeor evaporate
     !    1. aggregation mass com,or rime mass (keeping the same ratio)
     !    2. crystal mass com.
     implicit none
     real(PS),intent(in) :: dM
     real(PS),dimension(*)  :: mass_comp,left_mass
     real(PS) :: rat
     real(PS) :: con
     ! fraction of crystal mass to the total aggregates
     real(PS) :: f_ice
     integer,intent(inout) :: ierr
     ! minimum mass of hexagonal ice crystal with 1 um radius
     real(PS),parameter :: m_icmin=4.763209003e-12

     ierr=0
     f_ice=mass_comp(imc)/(mass_comp(imc)+mass_comp(ima))
     rat=mass_comp(imr)/(mass_comp(imc)+mass_comp(ima)+mass_comp(imr))

     left_mass(imr_m)=left_mass(imr_m)+dM*rat
     left_mass(ima_m)=left_mass(ima_m)+dM*(1.0_PS-rat)*(1.0_PS-f_ice)
     left_mass(imc_m)=left_mass(imc_m)+dM*(1.0_PS-rat)*f_ice

     left_mass(imr_m)=max(0.0_PS,left_mass(imr_m))
     left_mass(ima_m)=max(0.0_PS,left_mass(ima_m))

     if( left_mass(imc_m)-m_icmin*con<-1.0e-4*m_icmin*con ) then
        left_mass(imc_m)=m_icmin*con
        if(mass_comp(ima)+mass_comp(imr)>=1.0e-25_PS) then
           rat=max(1.0_PS,min(0.0_PS,mass_comp(imr)/(mass_comp(ima)+mass_comp(imr))))
           left_mass(imr_m)=left_mass(imr_m)+(dM+(mass_comp(imc)-left_mass(imc_m)))*rat
           left_mass(ima_m)=left_mass(ima_m)+(dM+(mass_comp(imc)-left_mass(imc_m)))*(1.0_PS-rat)

           left_mass(imr_m)=max(0.0_PS,left_mass(imr_m))
           left_mass(ima_m)=max(0.0_PS,left_mass(ima_m))
        else
           left_mass(imr_m)=0.0_PS
           left_mass(ima_m)=0.0_PS
        endif
     end if

     if( left_mass(imr_m)<=1.0e-25_PS.and.left_mass(ima_m)<=1.0e-25_PS.and.&
         left_mass(imc_m)<=1.0e-25_PS.and.mass_comp(imt)>1.0e-22 ) then
        ierr=1
     end if

  end subroutine dep_mass3_vec2
```

---

## 3. Coefficient functions (`class_Group.F90`)

### 3a. `cal_capacitance_vec`

`class_Group.F90:7288`. **LIQUID branch** = `phase==1` (lines 7306–7317): `CAP = 0.5*len` (sphere). Ice branch `phase==2` handles spheroids/rosettes/graupel.

```fortran
  subroutine cal_capacitance_vec(level,phase,g,ag,icond1)
    ! +++ Calculate capacitance based on the spheroid assumption,
    !     for example, see Chen and Lamb (1994a).                     +++
    ! phase of hydrometeor
    ! 1: liquid, 2: solid
    integer, intent(in)  :: phase
    integer, intent(in)  :: level
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    real(PS)                    :: eps, d, phi
    ! coefficient for capacitance by Chiruta and Wang (2003)
    real(PS),parameter   :: CAP_ros=0.619753727
    integer                     :: i,n

    if( phase == 1 ) then
       ! --- in case of liquid phase ---
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          g%MS(i,n)%CAP = 0.5_PS*g%MS(i,n)%len
        endif
      enddo
      enddo
    else if( phase == 2 ) then
      ! --- in case of solid phase ---
      ! +++ first caclulate for a hexagonal monocrystal +++
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          if( g%IS(i,n)%phi_ic < 1.0_PS ) then
            ! --- in case of oblate spheroids ---
            d = sqrt( max(g%MS(i,n)%a_len**2 - g%MS(i,n)%c_len**2, 0.0_PS) )
            eps = sqrt( 1.0_PS - g%IS(i,n)%phi_ic**2)
            g%MS(i,n)%CAP_hex = d/asin(eps)
          else if( g%IS(i,n)%phi_ic> 1.0_PS ) then
            ! --- in case of prolate spheroids ---
            d = sqrt( max(g%MS(i,n)%c_len**2 - g%MS(i,n)%a_len**2, 0.0_PS) )
            eps = sqrt( 1.0_PS - 1.0_PS/g%IS(i,n)%phi_ic**2)
            g%MS(i,n)%CAP_hex = d/log((1.0_PS + eps)*g%IS(i,n)%phi_ic )
          else
            ! --- in case of sphere ---
            g%MS(i,n)%CAP_hex = g%MS(i,n)%a_len
          end if

          ! ++++ second calculate the capacitance for polycrystals ++++
          if( g%IS(i,n)%sh_type <= 2 ) then
            ! --- case of ice crystals ---
            if(g%IS(i,n)%is_mod(2)==1) then
              if(g%IS(i,n)%habit<=3) then
                g%MS(i,n)%CAP=g%MS(i,n)%CAP_hex
              elseif(g%IS(i,n)%habit==4) then
                g%MS(i,n)%CAP = CAP_ros*g%IS(i,n)%r
              elseif(g%IS(i,n)%habit==5) then
                eps=0.968245837_PS
                d=g%IS(i,n)%e*eps
                g%MS(i,n)%CAP=d/asin(eps)
              else
                eps=0.968245837_PS
                d=max(g%IS(i,n)%e,g%IS(i,n)%r)*eps
                g%MS(i,n)%CAP=d/asin(eps)
              end if
            else
              g%MS(i,n)%CAP = g%MS(i,n)%semi_a
            endif
          elseif( g%IS(i,n)%sh_type <= 3 ) then
            ! --- case of aggregates and rimed aggregates ---
            g%MS(i,n)%CAP=(g%IS(i,n)%V_cs/coefpi6)**(1.0/3.0)*0.25_PS
          else
            ! --- case of graupels ---
            phi = g%IS(i,n)%phi_cs
            if( phi < 1.0_PS ) then
              d = sqrt( g%MS(i,n)%semi_a**2.0 - g%MS(i,n)%semi_c**2.0)
              eps = sqrt( 1.0_PS - phi**2.0)
              g%MS(i,n)%CAP = d/asin(eps)
            else if( phi > 1.0_PS ) then
              d = sqrt( g%MS(i,n)%semi_c**2.0 - g%MS(i,n)%semi_a**2.0)
              eps = sqrt( 1.0_PS - phi**(-2.0))
              g%MS(i,n)%CAP = d/log((1.0_PS + eps)*phi )
            else
              g%MS(i,n)%CAP = g%MS(i,n)%semi_a
            end if
          end if
          if(level<=5.and.ag%TV(n)%T>=T_0) then
            g%MS(i,n)%CAP=max(g%MS(i,n)%CAP,(g%MS(i,n)%mean_mass/coef4pi3)**(1.0/3.0))
          else
            if(ag%TV(n)%T>=T_0) then
              g%MS(i,n)%CAP=max(g%MS(i,n)%CAP,(g%MS(i,n)%mass(imw)/g%MS(i,n)%con/coef4pi3)**(1.0/3.0))
            endif
          end if
        endif
      enddo
      enddo
    end if
  end subroutine cal_capacitance_vec
```

### 3b. `cal_ventilation_coef_vec`

`class_Group.F90:7575`. **LIQUID branch** = `phase==1` (lines 7605–7653): Schmidt/Prandtl numbers, spherical-drop Reynolds number, Hall & Pruppacher (1976) `fv`/`fh`, kinetic `fkn`. Ice branch `phase==2` follows.

```fortran
  subroutine cal_ventilation_coef_vec(phase,g,ag,icond1)
    use class_Thermo_Var, only: &
       get_fkn
    ! Calculate ventilation coefficient
    ! based on empirical expressions by Hall and Pruppacher (1976).
    integer, intent(in)  :: phase
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    real(PS)                   :: N_sc
    real(PS)                   :: N_ns
    real(PS)                   :: Xv, Xh
    real(PS)                   :: r_m
    integer :: i,n

    if( phase == 1 ) then
      ! --- in case of liquid phase ---
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          ! +++ calculate the Schmidt number +++
          N_sc = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%D_v
          ! +++ calculate the Nusselt or Prandtl number +++
          N_ns = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%k_a
          ! assume the spherical water drops
          ! +++ calculate Reynolds number +++
          g%MS(i,n)%Nre = g%MS(i,n)%len*g%MS(i,n)%vtm*ag%TV(n)%den&
                          /ag%TV(n)%d_vis
          ! +++ calculate dimensionless number X +++
          Xv = (N_sc**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)
          Xh = (N_ns**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)

          if( Xv < 1.4_PS ) then
            g%MS(i,n)%fv = 1.0_PS + 0.108_PS*(Xv**2.0)
          else if( 1.4_PS <= Xv .and. Xv <= 51.4_PS ) then
            g%MS(i,n)%fv = 0.78_PS + 0.308_PS*Xv
          else
            g%MS(i,n)%fv = 0.78_PS + 0.308_PS*51.4_PS
          end if

          if( Xh < 1.4_PS ) then
            g%MS(i,n)%fh = 1.0_PS + 0.108_PS*(Xh**2.0)
          else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
            g%MS(i,n)%fh = 0.78_PS + 0.308_PS*Xh
          else
            g%MS(i,n)%fh = 0.78_PS + 0.308_PS*51.4_PS
          end if

          ! +++ calculate kinetic effect +++
          r_m=g%MS(i,n)%len*0.5_PS
          g%MS(i,n)%fkn = get_fkn(ag%TV(n),phase,r_m)
        endif
      enddo
      enddo

    else if( phase == 2 ) then
      ! --- in case of ice phase ---
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          N_sc = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%D_v
          N_ns = ag%TV(n)%d_vis/ag%TV(n)%den/ag%TV(n)%k_a
          g%MS(i,n)%Nre = g%MS(i,n)%len*g%MS(i,n)%vtm*ag%TV(n)%den&
                         /ag%TV(n)%d_vis
          Xv = (N_sc**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)
          Xh = (N_ns**(1.0/3.0))*sqrt(g%MS(i,n)%Nre)

          if(g%IS(i,n)%sh_type<=2) then
            if( Xv < 1.0_PS ) then
              g%MS(i,n)%fac = (1.0_PS + 0.14_PS*(Xv**2.0)*&
                  sqrt(g%MS(i,n)%c_len/g%MS(i,n)%len/2.0_PS))/&
                  (1.0_PS + 0.14_PS*(Xv**2.0)*&
                  sqrt(g%MS(i,n)%a_len/g%MS(i,n)%len/2.0_PS))
            else
              g%MS(i,n)%fac = (0.86_PS + 0.28_PS*Xv*&
                  sqrt(g%MS(i,n)%c_len/g%MS(i,n)%len/2.0_PS))/&
                  (0.86_PS + 0.28_PS*Xv*&
                  sqrt(g%MS(i,n)%a_len/g%MS(i,n)%len/2.0_PS))
            end if
          else
            g%MS(i,n)%fac = 1.0_PS
          end if

          if(g%IS(i,n)%sh_type<=2.and.g%IS(i,n)%habit==4) then
            ! bullet rosettes: use Lie et al (2003)
            g%MS(i,n)%fv = 1.0_PS + 0.3005_PS*Xv - 0.0022_PS*Xv*Xv
            g%MS(i,n)%fh = 1.0_PS + 0.3005_PS*Xh - 0.0022_PS*Xh*Xh
          else
            if( Xv < 1.0_PS ) then
              g%MS(i,n)%fv = 1.0_PS + 0.14_PS*(Xv**2.0)
            else
              g%MS(i,n)%fv = 0.86_PS + 0.28_PS*Xv
            end if
            if( Xh < 1.0_PS ) then
              g%MS(i,n)%fh = 1.0_PS + 0.14_PS*(Xh**2.0)
            else
              g%MS(i,n)%fh = 0.86_PS + 0.28_PS*Xh
            end if
          end if

          ! +++ calculate kinetic effect +++
          g%MS(i,n)%fkn=1.0_PS
          r_m = sqrt(g%MS(i,n)%semi_a**2+g%MS(i,n)%semi_c**2)
          g%MS(i,n)%fkn=get_fkn(ag%TV(n),phase,r_m)

          if(g%MS(i,n)%inmlt==1.and.ag%TV(n)%T>=T_0) then
            if( Xh < 1.4_PS ) then
              g%MS(i,n)%fh = max(g%MS(i,n)%fh,1.0_PS + 0.108_PS*(Xh**2.0))
            else if( 1.4_PS <= Xh .and. Xh <= 51.4_PS ) then
              g%MS(i,n)%fh = max(g%MS(i,n)%fh,0.78_PS + 0.308_PS*Xh)
            else
              g%MS(i,n)%fh = max(g%MS(i,n)%fh,0.78_PS + 0.308_PS*51.4_PS)
            end if
          end if
        endif
      enddo
      enddo
    end if

  END subroutine cal_ventilation_coef_vec
```

### 3c. `cal_coef_vapdep2_vec` — growth coefficients `coef(1)`, `coef(2)`

`class_Group.F90:9348`. **LIQUID branch** = `phase==1` (lines 9394–9475): Köhler-based (Khvorostyanov & Curry 2014) condensation coefficients with kinetic (`buzai_con`) and psychrometric (`gamma_w`) corrections; also diagnoses critical/haze radius `r_crt`/`r_act`. Ice branch `phase==2` is the simple `4πD_v·CAP·MR·fv·fkn`.

```fortran
  subroutine cal_coef_vapdep2_vec(phase,g,ag,icond1,nu_aps,phi_aps,m_aps)
    use class_Mass_Bin, only: &
       get_critrad_anal, &
       get_hazerad_anal
    integer, intent(in)  :: phase
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(in)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    real(PS),dimension(*),intent(in) :: phi_aps
    real(PS),dimension(*),intent(in) :: nu_aps,m_aps
    real(PS) :: den_ap, r_n
    real(PS) :: r_n3
    real(PS) :: AA,sb,beta
    real(PS) :: rho_s
    real(PS) :: SRW=0.99_PS
    real(PS) :: rd_c
    real(PS) :: s_salt
    real(PS) :: vw
    real(PS) :: buzai_con
    real(PS) :: gamma_w
    integer                     :: i,n

    if(phase==1) then
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          ! This formulation is based on Khvorostyanov and Curry (2014).
          ! assuming that category 1 is the soluble material.

          ! calculate density of mixed aerosols
          den_ap=g%MS(i,n)%den_ai/(1.0_PS-g%MS(i,n)%eps_map*(1.0_PS-g%MS(i,n)%den_ai/g%MS(i,n)%den_as))
          r_n3=(g%MS(i,n)%mass(rmat)/g%MS(i,n)%con)/coef4pi3/den_ap
          r_n=r_n3**(1.0_PS/3.0_PS)

          ! unit is now in cm
          AA=2.0_PS*ag%TV(n)%sig_wa/(R_v*ag%TV(n)%T*den_w)

          sb=nu_aps(1)*g%MS(i,n)%eps_map*M_W*den_ap/(M_aps(1)*den_w)*phi_aps(1)
          beta=0.5_PS

          s_salt=  AA/g%MS(i,n)%a_len-sb*r_n**(2.0_PS*(1.0_PS+beta))/&
                   (g%MS(i,n)%a_len**3-r_n3)

          vw=sqrt(8.0_PS/PI*R_v*ag%TV(n)%T)
          buzai_con=4.0_PS*ag%TV(n)%D_v/(vw*a_cliq)

          rho_s=ag%TV(n)%e_sat(1)/(ag%TV(n)%T*R_v)

          gamma_w=1.0_PS+L_e*rho_s/(c_pa*ag%TV(n)%T*ag%TV(n)%den)* &
               (L_e/R_v/ag%TV(n)%T-1.0)

          g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%D_v*rho_s/gamma_w * &
                g%MS(i,n)%CAP*g%MS(i,n)%CAP/(g%MS(i,n)%CAP+buzai_con) *&
                g%MS(i,n)%fv

          g%MS(i,n)%coef(2)=-g%MS(i,n)%coef(1)*s_salt

          if(g%MS(i,n)%a_len>1.0e-7_PS.and.r_n>1.0e-7_PS) then
            !   calculation of critical radius and haze at SRW is done with
            !    formula by Khvorostyanov and Curry (2014)
            rd_c=get_critrad_anal(AA,sb,beta,r_n)  ! [cm]
            g%MS(i,n)%r_crt=rd_c
            ! calculate haze size with SRW. This is used for evaporation only.
            g%MS(i,n)%r_act=get_hazerad_anal(AA,sb,beta,SRW,r_n)
          else
            ! according to Chen 1992, JAS
            g%MS(i,n)%r_crt=r_n
            g%MS(i,n)%r_act=r_n
          endif
        endif
      enddo
      enddo
    elseif(phase==2) then
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          if(ag%TV(n)%T<T_0.and.g%MS(i,n)%inmlt==0) then
            ! dry growth
            g%MS(i,n)%coef(1)=  4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn
            g%MS(i,n)%coef(2)= -g%MS(i,n)%coef(1)
          else
            ! wet growth
            g%MS(i,n)%coef(1)=  4.0_PS*PI*ag%TV(n)%D_v*g%MS(i,n)%CAP*MR*g%MS(i,n)%fv*g%MS(i,n)%fkn
            g%MS(i,n)%coef(2)= -g%MS(i,n)%coef(1)
          end if
        end if
      enddo
      enddo
    end if
  end subroutine cal_coef_vapdep2_vec
```

### 3d. `cal_coef_vapdep_ap_vec` — aerosol/nucleation growth coefficient (phase==3 / ice-nucleation)

`class_Group.F90:9549`. This is the aerosol-group nucleation coefficient (called from `diag_pq` phase==3). Not liquid-hydrometeor growth, but included as requested.

```fortran
  subroutine cal_coef_vapdep_ap_vec(level,ica,mes_rc,g,ag,icond1)
    use class_Thermo_Var, only: &
       get_fkn
    use mod_amps_utility, only: &
       get_growth_mode
    integer,intent(in) :: level,ica
    integer,dimension(*),intent(in) :: mes_rc
    type (Group), intent(inout)  :: g
    type (AirGroup), intent(inout)   :: ag
    integer,dimension(g%N_BIN,g%L),intent(in)    ::  icond1
    real(PS),parameter   :: CAP_ros=0.619753727
    real(PS) :: fkn
    integer                     :: i,n

    if(level==3.or.level==5.or.level==7) then
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          if((mes_rc(n)==2.or.mes_rc(n)==4).and.ag%TV(n)%T<253.16_PS) then
          ! if liquid co exist
            ag%TV(n)%nuc_gmode=get_growth_mode(ag%TV(n)%T,ag%TV(n)%e_sat(1)/ag%TV(n)%e_sat(2)-1.0_PS)
          endif

          if(ica==2) then
          ! +++ calculate kinetic effect +++
            fkn=get_fkn(ag%TV(n),2,g%MS(i,n)%a_len)
            if(ag%TV(n)%T>-20.0+273.16) then
              g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
            else
              if(ag%TV(n)%nuc_gmode==2.or.ag%TV(n)%nuc_gmode==3) then
                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
              elseif(ag%TV(n)%nuc_gmode==4) then
                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
              else
                g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
              end if
            end if
          endif
        end if
      enddo
      enddo
    elseif(level<=2.or.level==4.or.level==6) then
      if(ica==2) then
         do n = 1, g%L
         do i = 1, g%N_BIN
          if(icond1(i,n)==0) then
            ! +++ calculate kinetic effect +++
            fkn=get_fkn(ag%TV(n),2,g%MS(i,n)%a_len)
            g%MS(i,n)%coef(1)=4.0_PS*PI*ag%TV(n)%gtp(2)*g%MS(i,n)%a_len*fkn
          end if
        enddo
        enddo
      endif
    end if
  end subroutine cal_coef_vapdep_ap_vec
```

---

## 4. Bin-shift / mass-space remap — `cal_transbin_vec` (`mod_amps_utility.F90:8713`)

This is the **mass-space remap** (gather): the shifted-bin distribution parameters `a2d` (rewritten as polynomial `n(m)=a0+a1·m+a2·m²+a3·m³`) are integrated over the overlap of each shifted bin `[left_bd,right_bd]` with each original bin `[binb8(ibx),binb8(ibx+1)]`, producing `trans_dN` (0th moment) and `trans_dM` (1st moment), accumulated into `new_N`, `new_M`, `new_Q`. The **`shift_bin_loop1` (iphase==1) is the LIQUID remap** that M2 must reproduce as a gather; `shift_bin_loop2` (iphase==2) is the ice remap (adds volume/length/axis quality-variable transfers via `den_ip_p`, `den_ic_p`, `axr_p`, `rag_p`, `rcg_p`, `n_exice_p`).

Arg list + declarations:

```fortran
subroutine cal_transbin_vec(iphase &
     ,L,nmass &
     ,nbin,nbin_s &
     ,binb &
     ,error_number &
     ,a2d,s_bd,mtend &
     ,new_N,new_M,new_Q &
     ,new_mtend &
     ,ratio_Mp,den_ip_p,axr_p,spx_p &
     ,habit_p,den_ic_p &
     ,rag_p,rcg_p,n_exice_p &
     ,actINF_p &
     ,iaer_src,ap_dN,ap_dM,ap_dMS,ap_dNI,ap_dMV,inevp)
  use maxdims
  use par_amps
  use com_amps, only: fid_alog
  implicit none
  !  calculate concentration and mass transferred to the original bins
  !    the distribution parameters a2d will be modified to be
  !     n(m)=a0+a1*m+a2*m*m+a3*m*m*m
  integer,intent(in) :: iphase,nmass,L
  integer,intent(in) :: nbin
  integer,intent(in) :: nbin_s
  real(PS),dimension(*),intent(in) :: binb ! precision test (PS)
  real(8),dimension(mxnbin+1,L,4), intent(inout) :: a2d
  real(RP), dimension(mxnbin+1,L,2), intent(in) ::  s_bd
  real(PS), dimension(mxnbin+1,*), intent(in)     :: mtend
  real(8), dimension(mxnbin,L), intent(inout) :: new_N
  real(8), dimension(mxnbin,L,1+mxnmasscomp), intent(inout) :: new_M
  real(8), dimension(mxnbin,L,mxnnonmc), intent(inout) :: new_Q
  real(8), dimension(mxnbin,*), intent(inout)         :: new_mtend
  real(RP), dimension(mxnbin+1,L,mxnmasscomp), intent(in) :: ratio_Mp
  real(PS), dimension(mxnbin+1,*), intent(in) :: den_ip_p,den_ic_p
  real(PS), dimension(mxnbin+1,L,mxnaxis-1), intent(in) :: axr_p
  integer, dimension(mxnbin+1,*), intent(in) :: habit_p
  real(PS), intent(in) :: spx_p
  real(PS), dimension(mxnbin+1,*), intent(in) :: rag_p,rcg_p
  real(RP), dimension(mxnbin+1,*), intent(in) :: n_exice_p
  real(PS), dimension(mxnbin+1,*), intent(in) :: actINF_p
  integer,dimension(mxnbin+1,*) :: error_number
  integer,intent(in) :: iaer_src
  real(8),dimension(LMAX,*),intent(inout),optional  :: ap_dN, ap_dM
  real(8),dimension(LMAX,*),intent(inout),optional  :: ap_dMS, ap_dNI, ap_dMV
  integer,dimension(mxnbin,*),intent(inout),optional :: inevp
  real(8) :: trans_dN, trans_dM
  real(8) :: trans_dM_dum
  real(8) :: trans_dvcs
  real(8) :: trans_dL(mxnaxis)
  real(PS) :: mid_mass
  real(8),dimension(mxnbin+1,LMAX) :: left_bd, right_bd
  real(8) :: aleft_in,aright_in,aboth_over,bd1,bd2
  real(PS) :: d_N, d_M
  real(8) :: xA,xB,xC,xD
  integer :: em
  integer :: ibx,k,i,n,in,icat
  integer dstype
  integer, parameter :: npreaxis=3
  integer :: i_epslt
  integer, dimension(mxnbin,LMAX) :: icond1
  real(8), dimension(mxnbin,LMAX) :: trans_dM_all
  real(8) :: C1,axr1,axr2
  real(8),dimension(mxnbin+1) :: binb8
  real(PS) :: eps_map
  real(PS),parameter :: lmt_frac=0.99999e-5
```

Boundary setup (shared, lines 8846–8897) — rewrites `a2d` into polynomial `a0 + a1·m` for the linear case and sets `left_bd`/`right_bd`:

```fortran
  ! distribution type
  dstype=1

  do ibx=1,nbin+1
    binb8(ibx)=binb(ibx)
  enddo
  ! 1. set up the left and right boundaries of new bins
  do in=1,nbin_s*L
    n=(in-1)/NBIN_s+1
    i=in-(n-1)*NBIN_s

    if(error_number(i,n)/=10.and.a2d(i,n,4)<=-9.98d+100) then
      ! linear distribution
      if( a2d(i,n,1) > 0.0d+0 ) then
        ! --- case of non-negative bin ---
        left_bd(i,n)=s_bd(i,n,1)
        right_bd(i,n)=s_bd(i,n,2)
      elseif( a2d(i,n,1) == -1.0d+0 ) then
        ! case of n(x'_1) < 0.0
        left_bd(i,n) = a2d(i,n,2)
        right_bd(i,n) = s_bd(i,n,2)
      else if( a2d(i,n,1) == -2.0d+0 ) then
        ! case of n(x'_2) < 0.0
        left_bd(i,n) = s_bd(i,n,1)
        right_bd(i,n) = a2d(i,n,2)
      else
        left_bd(i,n)=0.0_PS
        right_bd(i,n)=0.0_PS
      end if
      ! re-write the parameters as n(m)=a0+a1*m
      a2d(i,n,1)=max(0.0d+0,a2d(i,n,1))-a2d(i,n,2)*a2d(i,n,3)
      a2d(i,n,2)=a2d(i,n,3)
      a2d(i,n,3)=0.0d+0
      a2d(i,n,4)=0.0d+0
    elseif(error_number(i,n)/=10.and.a2d(i,n,4)>-9.98d+100) then
      ! cubic distribution
      left_bd(i,n)=s_bd(i,n,1)
      right_bd(i,n)=s_bd(i,n,2)
    else
      left_bd(i,n)=0.0_PS
      right_bd(i,n)=0.0_PS
    endif
  enddo
```

### LIQUID remap — `shift_bin_loop1` (`iphase == 1`), lines 8899–9089 verbatim

```fortran
  if( iphase == 1 ) then
    !     --- case where the colletor is a liquid hydromteor ---

    shift_bin_loop1: do i=1,nbin_s

      icond1(:,:) = 0
      trans_dm_all(:,:) = 0

      do n = 1, L
      do ibx = 1, nbin

        aleft_in=(left_bd(i,n)-binb8(ibx))*(binb8(ibx+1)-left_bd(i,n))
        aright_in=(right_bd(i,n)-binb8(ibx))*(binb8(ibx+1)-right_bd(i,n))
        aboth_over=(binb8(ibx)-left_bd(i,n))*(right_bd(i,n)-binb8(ibx+1))

        if(error_number(i,n)/=10.and.&
           (aleft_in>0.0_RP.or.aright_in>0.0_RP.or.aboth_over>=0.0_RP)) then

          icond1(ibx,n)=1

          bd1=max(binb8(ibx), left_bd(i,n))
          bd2=min(binb8(ibx+1), right_bd(i,n))
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5_DP*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0_DP &
                  +0.25_DP*a2d(i,n,4)*xD*xA*xB

          trans_dM=0.5_DP*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0_DP &
                  +0.25_DP*a2d(i,n,3)*xD*xA*xB &
                  +0.2_DP*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)

          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          new_N(ibx,n) = new_N(ibx,n) + trans_dN
          new_M(ibx,n,rmt) = new_M(ibx,n,rmt) + trans_dM
          trans_dM_all(ibx,n)=trans_dM_all(ibx,n)+trans_dM
          new_mtend(ibx,n) = new_mtend(ibx,n) + mtend(i,n)*trans_dN

        endif
      enddo
      enddo

      ! case of mass transfering above the max bin boundary,
      !   put them into the largest bin
      do n=1,L
        if(error_number(i,n)/=10.and.&
           right_bd(i,n)>binb8(nbin+1) ) then
          ibx=nbin

          icond1(ibx,n)=1

          bd1=max(binb8(ibx+1), left_bd(i,n))
          bd2=right_bd(i,n)
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5d+0*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,4)*xD*xA*xB

          trans_dM=0.5d+0*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,3)*xD*xA*xB &
                  +0.2d+0*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)

          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          new_N(ibx,n) = new_N(ibx,n) + trans_dN
          new_M(ibx,n,rmt) = new_M(ibx,n,rmt) + trans_dM
          trans_dM_all(ibx,n)=trans_dM_all(ibx,n)+trans_dM
          new_mtend(ibx,n) = new_mtend(ibx,n) + mtend(i,n)*trans_dN

        endif
      enddo

    if(iaer_src.eq.1) then
      ! case of mass transfering below the min bin boundary,
      !   put them into aerosol variables
      do n=1,L
        if(error_number(i,n)/=10.and.&
           left_bd(i,n)<binb8(1) ) then
          ibx=1

          bd1=left_bd(i,n)
          bd2=min(binb8(1),right_bd(i,n))
          em=0

          xA=bd2+bd1
          xB=bd2-bd1
          xC=bd2*bd1
          xD=bd2*bd2+bd1*bd1
          trans_dN=a2d(i,n,1)*xB &
                  +0.5d+0*a2d(i,n,2)*xA*xB &
                  +a2d(i,n,3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,4)*xD*xA*xB

          trans_dM=0.5d+0*a2d(i,n,1)*xA*xB &
                  +a2d(i,n,2)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a2d(i,n,3)*xD*xA*xB &
                  +0.2d+0*a2d(i,n,4)*xB*(xA*xA*(xD-xC)+xC*xC)

          trans_dN=max(0.0e+0_DS,trans_dN)
          trans_dM=max(0.0e+0_DS,trans_dM)

          eps_map=ratio_Mp(i,n,rmas_m)/ratio_Mp(i,n,rmat_m)
          i_epslt=0.5*(1.0-sign(1.0_PS,eps_map-lmt_frac))
          icat=i_epslt*2+(1-i_epslt)*1
          ap_dN(n,icat)=ap_dN(n,icat)+trans_dN
          ap_dM(n,icat)=ap_dM(n,icat)+trans_dM*ratio_Mp(i,n,rmat_m)
          ap_dMS(n,icat)=ap_dMS(n,icat)+trans_dM*ratio_Mp(i,n,rmas_m)
          ap_dNI(n,icat)=ap_dNI(n,icat)+actINF_p(i,n)*trans_dN
          ap_dMV(n,icat)=ap_dMV(n,icat)+trans_dM
          inevp(i,n)=icat

        endif
      enddo
    endif

    do k = 1,nmass
       do n = 1, L
       do ibx = 1, nbin
          if(icond1(ibx,n)==1) then
             new_M(ibx,n,k+1) = new_M(ibx,n,k+1) &
                              + ratio_Mp(i,n,k) * trans_dM_all(ibx,n)
          endif
       enddo
       enddo
    enddo

    enddo shift_bin_loop1

  else if(iphase==2) then
```

`shift_bin_loop2` (iphase==2, ICE) follows at lines 9091–9311 — same `trans_dN`/`trans_dM` moment integrals but with `new_M(...,imt)`, plus volume-of-circumscribing-sphere `trans_dvcs = trans_dM_dum/den_ip_p`, axis lengths `trans_dL(1:3)` from `den_ic_p`/`axr_p`, and `iag`/`icg`/`inex` quality-variable transfers using `rag_p`/`rcg_p`/`n_exice_p`. Closing `end subroutine cal_transbin_vec` at 9313.

---

## 5. `diag_pq` — per-bin derived quantities, LIQUID branch

`mod_amps_core.F90:12552`. Full arg list + declarations, common init, then **the LIQUID `if( phase == 1 )` branch verbatim** (lines 12786–12835). It computes: insoluble/soluble mass split & `eps_map`; density + a/c lengths (`cal_den_aclen_vec`); terminal velocity (spheroid mode 1); ventilation coef; capacitance; vapor-dep coefficients (`cal_coef_vapdep2_vec`); surface temperature.

Arg list + declarations:

```fortran
  subroutine diag_pq( g, ag, level, em0, mes_rc,&
       ID,JD,KD,rdsd,ihabit_gm_random, &
       eps_ap0,nu_aps,phi_aps,M_aps,&
       ap_sig_cp,ap_mean_cp,cdf_cp_0,isnrml,&
       icat,flagp,ap_lnsig,ap_mean)
    use scale_prc, only: PRC_abort
    use class_Mass_Bin, only: &
       data1d_lut_big, cal_meanmass_vec, cal_mass_comp_ap, cal_melt_index, &
       fix_mass_comp, diag_sh_type_v6, diag_pardis_ap
    use class_Thermo_Var, only: get_sat_vapor_pres_lk
    use class_Ice_Shape, only: ini_Ice_Shape_v4, cal_Ice_Shape_v4, diag_habit_v4
    use class_Group, only: &
       cal_den_aclen_vec, cal_capacitance_vec, cal_ventilation_coef_vec, &
       cal_growth_mode_vec, cal_terminal_vel_vec, cal_coef_vapdep2_vec, &
       cal_coef_vapdep_ap_vec, cal_surface_temp2_vec
    use mod_amps_utility, only: get_len_s1, get_len_c2a, get_len_s3, random_genvar
    ! Diagnose physical quantities of each bin
    type (Group), intent(inout)    :: g
    type (AirGroup),intent(inout)     :: ag
    integer   :: level
    integer    :: em0
    integer,dimension(g%n_bin,g%L)    :: em,em2
    integer,dimension(g%n_bin,g%L)    :: icond1
    integer,dimension(*)   :: mes_rc
    integer :: ID(*),JD(*),KD(*)
    real(PS),intent(in)   :: eps_ap0
    real(PS),dimension(*),intent(in)  :: phi_aps
    real(PS),dimension(*),intent(in) :: nu_aps, m_aps
    integer,intent(in),optional :: icat
    integer,intent(in),optional :: flagp
    real(PS),optional :: ap_lnsig(*),ap_mean(*)
    real(PS),optional :: ap_sig_cp(*),ap_mean_cp(*)
    real(PS),optional :: cdf_cp_0(*)
    type(data1d_lut_big),intent(in),optional :: isnrml
    type(random_genvar),intent(inout) :: rdsd
    integer, intent(in)           :: ihabit_gm_random
    real(PS),dimension(1+mxnmasscomp,g%n_bin,g%L) :: dum_mass
    real(PS),dimension(g%n_bin,g%L) :: semi_a_i,semi_c_i
    real(PS) :: vice1,vice2
    real(PS) :: xlen_s1,xlen_c2a,xlen_s3
    integer   :: phase
    real(PS),parameter :: V_csmin=4.18879020478639d-12,V_csmax=4.18879020478639d+3
    real(PS),parameter :: m_icmin=4.763209003d-12
    integer   :: i,j,n

    em0=0
```

Phase determination + common per-bin initialization + mean-mass (shared, lines 12692–12784):

```fortran
    ! determine the phase of hydrometeor
    if( g%token == 1 .or. g%token == 11 .or. g%token == 12) then
       phase = 1
    else if( g%token == 2 ) then
       phase = 2
    else if( g%token == 3 ) then
       phase = 3
    end if

    ! Initialization
    do n = 1, g%L
    do i = 1, g%N_BIN
      em(i,n)=0
      icond1(i,n)=0
      g%MS(i,n)%den=1.0_PS
      g%MS(i,n)%len=0.0_PS
      g%MS(i,n)%semi_a=0.0_PS
      g%MS(i,n)%semi_c=0.0_PS
      g%MS(i,n)%inmlt=0
      g%MS(i,n)%inevp=0
      g%MS(i,n)%eps_map=eps_ap0
      g%MS(i,n)%vtm=0.0_PS
      g%MS(i,n)%Nre=0.0_PS
      g%MS(i,n)%fv=1.0_PS
      g%MS(i,n)%fh=1.0_PS
      g%MS(i,n)%fkn=1.0_PS
      g%MS(i,n)%fac=1.0_PS
      g%MS(i,n)%CAP=0.0_PS
      g%MS(i,n)%CAP_hex=0.0_PS
      g%MS(i,n)%coef(1)=0.0_PS
      g%MS(i,n)%coef(2)=0.0_PS
      g%MS(i,n)%r_crt=0.0_PS
      g%MS(i,n)%r_act=0.0_PS
      g%MS(i,n)%th00_cp=0.0_PS
      g%MS(i,n)%tmp = ag%TV(n)%T
      g%MS(i,n)%e_sat=get_sat_vapor_pres_lk(phase,g%MS(i,n)%tmp,ag%estbar,ag%esitbar)
    enddo
    enddo

    do n = 1, g%L
    do i = 1, g%N_BIN
      if(g%MS(i,n)%con>1.0e-30_PS.and.g%MS(i,n)%mass(1)>1.0e-30_PS) then
        ! +++ calculate mean mass +++
        call cal_meanmass_vec( g%MS(i,n),em(i,n))
        if(em(i,n)>0) then
          icond1(i,n)=1
          g%MS(i,n)%con=0.0_PS
          g%MS(i,n)%mass(1)=0.0_PS
          g%MS(i,n)%mean_mass=0.0_PS
        endif
      else
        icond1(i,n)=1
        g%MS(i,n)%con=0.0_PS
        g%MS(i,n)%mass(1)=0.0_PS
        g%MS(i,n)%mean_mass=0.0_PS
      endif
    enddo
    enddo
    do n = 1, g%L
    do i = 1, g%N_BIN
       if(icond1(i,n)/=0) then
          do j=1,g%N_masscom
             g%MS(i,n)%mass(1+j)=0.0_PS
          enddo
       endif
    enddo
    enddo
```

**LIQUID branch** `if( phase == 1 )`, lines 12786–12835 verbatim:

```fortran
    if( phase == 1 ) then

       do n = 1, g%L
       do i = 1, g%N_BIN

        if(icond1(i,n)==0) then

          ! +++ Diagnose the insoluble mass and soluble mass fraction +++
          g%MS(i,n)%mass(rmai)=max(g%MS(i,n)%mass(rmat)-g%MS(i,n)%mass(rmas),0.0_PS)
          g%MS(i,n)%eps_map=max(0.0_PS,min(1.0_PS,g%MS(i,n)%mass(rmas)/g%MS(i,n)%mass(rmat)))

          ! +++ Diagnose the bulk density +++
          ! +++ Diagnose the a and c lengths +++
        endif
      enddo
      enddo

      call cal_den_aclen_vec(g,icond1)

      ! +++ Calculate terminal velocity +++
          ! always mode is 1 (spheroid assumption) for liquid.
      call cal_terminal_vel_vec(phase,g,ag,icond1,1)

      ! +++ Calculate ventilation coefficient +++
      call cal_ventilation_coef_vec(phase,g,ag,icond1)

      ! +++ Calculate capacitance +++
      call cal_capacitance_vec(level,phase,g,ag,icond1)

      ! +++ precalculation for vapor deposition +++
      call cal_coef_vapdep2_vec(phase,g,ag,icond1,nu_aps,phi_aps,m_aps)

      ! +++ Calculate the surface temperature of each hydrometeor +++
      call cal_surface_temp2_vec( phase, g,ag,icond1)

    else if( phase == 2 ) then
```

The `phase==2` (ice) branch (lines 12837–13148) does shape diagnosis (`cal_Ice_Shape_v4`, `diag_habit_v4`, `ini_Ice_Shape_v4`), mass/con fixing, `cal_mass_comp_ap`, `cal_melt_index`, `fix_mass_comp`, `diag_sh_type_v6`, `cal_cs_spheroid3_vec`, `cal_bulk_density3_vec`, terminal velocity modes 1 and 2, ventilation, capacitance, `cal_coef_vapdep2_vec`, surface temp, and `cal_growth_mode_vec` (habit/growth-mode — ice only). The `phase==3` branch (lines 13150–13195) does aerosol bulk density, `diag_pardis_ap`, and `cal_coef_vapdep_ap_vec`. `end subroutine diag_pq` at 13199.

---

### Summary of liquid-vs-ice branching

| Routine | LIQUID selector | ICE-only work skipped for liquid |
|---|---|---|
| `vapor_deposition` | `g%token==1` (`jmat=rmat`); growth kernel (456–961) is under `if(g%token==2)` | inherent growth `gamma`, `cal_ex_vapor_density_vec`, `assign_Qp_v3_vec`, `cal_xxx_p_v5_vec`, habit debug block |
| `cal_ratio_mass_vd_vec` | `if(g%token==1)` — only total/aerosol mass-fraction split (16582–16598) | axis growth, `acd_mode`, `dend_zone`, growth-mode, wet/dry, `dep_mass3_vec2`, level≤5/≥6 component ratios |
| `cal_capacitance_vec` | `phase==1`: `CAP=0.5*len` | spheroid/rosette/graupel capacitance |
| `cal_ventilation_coef_vec` | `phase==1`: spherical-drop `fv`/`fh`/`fkn` | `fac` shape factor, rosette Lie(2003), melt correction |
| `cal_coef_vapdep2_vec` | `phase==1`: Köhler (K&C 2014) `coef(1/2)`, `r_crt`/`r_act` | ice `4πD_v·CAP·MR·fv·fkn` dry/wet |
| `cal_transbin_vec` | `iphase==1` (`shift_bin_loop1`): `new_N`/`new_M(rmt)` moment gather | volume/length/axis quality-variable transfers |
| `diag_pq` | `phase==1` (12786–12835): density, a/c len, vtm, fv, CAP, coef, sfc T | shape/habit diagnosis, growth-mode, mode-2 terminal velocity |