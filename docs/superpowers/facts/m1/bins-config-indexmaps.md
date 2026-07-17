# AMPS Bin-Grid + Configuration — Verbatim Fortran Extraction

## 1. `binmicrosetup_scale` — bin boundaries, haze split, sizing
**File:** `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_lib.F90`, lines 328–618 (lookup-table reads begin at line 621, excluded)

```fortran
subroutine binmicrosetup_scale(npr,nbr,ncr,npi,nbi,nci, &
                               npa,nba,nca,model_k, model_i, model_j,nbhzcl, &
                               estbar,esitbar)
  use scale_prc, only: &
     PRC_abort
  use maxdims
  use com_amps
  use mod_amps_utility, only: &
       read_AMPSTASK, &
       RDCETB, &
       RDAPTB, &
       RDSTTB, &
       init_inherent_growth_par, &
       init_osmo_par, &
       init_normal_lut, &
       init_inv_normal_lut, &
       read_seed
  use par_amps, only: &
       isplit_bin_liq, &
       isplit_bin_ice
  implicit none
!  arguments
      integer, intent(in)  :: npr, nbr, ncr, npi, nbi, nci, npa, nba, nca, model_k, model_i, model_j
      integer, intent(out) :: nbhzcl
      real(DS), DIMENSION(*) ::  ESTBAR
      real(DS), DIMENSION(*) ::  ESITBAR

      integer :: i, j

      integer :: ierr

!
!     IXCNFL ...... level of complexity for X hydrometeor
!     default total number of bins
      integer :: idbnum(10)
      DATA idbnum/0, 3, 5, 10, 15, 20, 30, 40, 60, 80/
!
!     setting for terminal velocity
      real(PS) :: VRR(11)
      DATA VRR/8.545e-1,8.59895E-01,8.65272E-01,8.75851E-01 &
           ,8.92279E-01  &
           ,9.11267E-01,9.28026E-01,9.40583E-01,9.49129E-01 &
           ,9.55522E-01,9.58816E-01/
      real(PS) :: VSR(5)
      DATA VSR/6.12185E-01,6.30215E-01,7.35575E-01,8.51963E-01, &
           8.95877E-01/

!     pi/6.0, 4.0*pi/3.0, 2.0*pi
      real(PS), parameter :: coef2=0.523598776, coef3=4.18879020478639, coef4=6.28318530717959
      real(PS) :: m0, m1

!     the minimum cloud droplet mass
      real(PS) :: c_mnm, c_max
!     1 um
!      parameter(c_mnm=4.188790205e-12)
!     0.1 um
      parameter(c_mnm=4.188790205e-15,c_max=6.54498e-8)
!     0.01 um
!     c_max -> rad=25 micron
!      parameter(c_mnm=4.188790205e-18,c_max=6.54498e-8)

      real(PS) :: dsrat, maxmass_r, dsrat_h
!     corresponds to 1 cm diameter of a rain drop
      parameter(maxmass_r=5.2359870E-01)

!     the mass boundaries to define pristince crystal, aggregates, hail
      real(PS) :: mg1, mg2, mg3, mg4
!ccc      parameter(mg1=4.18879020478639E-12,mg2=1.0e-7,mg3=1.0e-1
!ccc     *     ,mg4=1.0e+1)
      parameter(mg1=4.18879020478639E-12,mg2=1.0e-6,mg3=1.0e-2 &
           ,mg4=1.0e+1)
!ccc      parameter(mg1=4.18879020478639E-12,mg2=1.0e-6,mg3=1.0e+0
!ccc     *     ,mg4=1.0e+1)


      integer :: NBIN20(3), NBIN10(3), NBIN40(3)
!ccc      NBIN20=(/4,14,2/)
!ccc      NBIN10=(/2,7,1/)
!ccc      NBIN40=(/8,28,4/)
      NBIN20=(/4,10,6/)
      NBIN10=(/2,5,3/)
      NBIN40=(/8,20,12/)

!
!   Create the file id for log
!
      call set_ampslog

      !open(fid_alog,file=fname_ampslog)
      !write(*,*) "fid_alog",fid_alog

      ! set grid maximum limits
      n1mx = model_k
      n2mx = model_i
      n3mx = model_j
      nxpmax = n2mx
      nypmax = n3mx
      nzpmax = n1mx
      nzgmax = n4mx
      mxln = max(n2mx,n1mx,n3mx)
      mxlh = max(n2mx,n3mx)
      mxlv = n1mx

      LMAX = n1mx + 1

      ! set 6D microphysical variable maximum limits
      if (nbr /= 40 .and. nbr /= 80) then
         LOG_ERROR("binmicrosetup_scale",*) "# of bins for liquid is not 40 or 80", nbr
         call PRC_abort
      endif
      if (nbi /= 10 .and. nbi /= 20 .and. nbi /= 40) then
         LOG_ERROR("binmicrosetup_scale",*) "# of bins for ice is not 10 or 20 or 40", nbi
         call PRC_abort
      endif
      if (nbr == 40) then
         allocate(bu_fd(2,17835))
         allocate(bu_tmass(435))
      else if (nbr == 80) then
         allocate(bu_fd(2,62400))
         allocate(bu_tmass(780))
      endif
      bu_fd = 0.0_PS
      bu_tmass = 0.0_PS

      nprmx = npr
      npimx = npi
      npamx = npa
      nbrmx = nbr
      nbimx = nbi
      nbamx = nba
      ncrmx = ncr
      ncimx = nci
      ncamx = nca

      mxnbinr = nbrmx
      mxnbini = nbimx
      mxnbina = nbamx

      mxnbin = max(nbrmx,nbimx)
      mxnbinb = mxnbin + 1
!
!  Read the AMPSTASK.F
!
         call read_AMPSTASK

!     add haze and cloud droplets category
!tmp tempei 2012/10         nbr = nbr+nbin_h
         nbhzcl=nbin_h
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "in binmic",nbr,nbin_h,nbhzcl
         end if
!
!     calculation of bin boundaries
         if(size(binbr).lt.nbr) then
            LOG_ERROR("binmicrosetup_scale",*) "# of bins for liquid is not enough in binbr", size(binbr), nbr
            call PRC_abort
         endif
         if(size(binbi).lt.nbi) then
            LOG_ERROR("binmicrosetup_scale",*) "# of bins for ice is not enough in binbi", size(binbi), nbi
            call PRC_abort
         endif

!     for rain category, including cloud cat.
         binbr(1) = c_mnm
!sheba         dsrat_h=(minmass_r/c_mnm)**(1.0/(nbin_h))
         dsrat_h=(c_max/c_mnm)**(1.0/(nbin_h))
         do i=2,nbin_h
            binbr(i)=binbr(i-1)*dsrat_h
         enddo
!sheba         binbr(nbin_h+1) = minmass_r
         binbr(nbin_h+1) = c_max

!sheba         dsrat=(maxmass_r/ &
!sheba              minmass_r*sth_r**(0.5* &
!sheba              real((nbr-nbin_h-1)*(nbr-nbin_h-0)))) &
!sheba              **(1.0/real(nbr-nbin_h-0))

         dsrat=(maxmass_r/c_max)**(1.0/(nbr-nbin_h))

         binbr(nbr+1)=maxmass_r

         if ( debug ) then
            write(fid_alog,*) "bin boundaries of rain"
            write(fid_alog,*) "c_mnm: ",c_mnm, c_max
!sheba         write(fid_alog,*) "minmass_r: ",minmass_r
            write(fid_alog,*) "maxmass_r: ",maxmass_r
            write(fid_alog,*) "nbin_h, dsrat_h: ",nbin_h,dsrat_h
            write(fid_alog,*) "theta: ",sth_r
            write(fid_alog,*) "bin#,binbr,srat,rad"
            do i=1,nbin_h+1
               write(fid_alog,'(I5,3ES15.6)') i,binbr(i),dsrat_h,(binbr(i)/coef3)**(1.0_RP/3.0_RP)
            end do
         end if
         do i=nbin_h+2,nbr+1
            binbr(i)=binbr(i-1)*dsrat
            if(debug) write(fid_alog,'(I5,3ES15.6)') i,binbr(i),dsrat,(binbr(i)/coef3)**(1.0_RP/3.0_RP)
!sheba            dsrat=dsrat/sth_r
         end do


!cc         do i=2,nbr+1
!cc            write(fid_alog,'(I5,2ES15.6)') i,binbr(i),binbr(i)/binbr(i-1)
!cc         end do
!ccc         do i=3,nbr+1
!ccc            binbr(i)=srat_r*binbr(i-1)+sadd_r
!ccc         end do


!     for ice category
         if(nbi==40) then
            srat_s=(mg2/mg1)**(1.0/(real(NBIN40(1),PS_KIND)))
            binbi(1) = minmass_s
            do i=2,NBIN40(1)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg3/mg2)**(1.0/(real(NBIN40(2),PS_KIND)))
            do i=NBIN40(1)+2,NBIN40(1)+NBIN40(2)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg4/mg3)**(1.0/(real(NBIN40(3),PS_KIND)))
            do i=NBIN40(1)+NBIN40(2)+2,NBIN40(1)+NBIN40(2)+NBIN40(3)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
         elseif(nbi==20) then
            srat_s=(mg2/mg1)**(1.0/(real(NBIN20(1),PS_KIND)))
            binbi(1) = minmass_s
            do i=2,NBIN20(1)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg3/mg2)**(1.0/(real(NBIN20(2),PS_KIND)))
            do i=NBIN20(1)+2,NBIN20(1)+NBIN20(2)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg4/mg3)**(1.0/(real(NBIN20(3),PS_KIND)))
            do i=NBIN20(1)+NBIN20(2)+2,NBIN20(1)+NBIN20(2)+NBIN20(3)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
         elseif(nbi==10) then
            srat_s=(mg2/mg1)**(1.0/(real(NBIN10(1),PS_KIND)))
            binbi(1) = minmass_s
            do i=2,NBIN10(1)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg3/mg2)**(1.0/(real(NBIN10(2),PS_KIND)))
            do i=NBIN10(1)+2,NBIN10(1)+NBIN10(2)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
            srat_s=(mg4/mg3)**(1.0/(real(NBIN10(3),PS_KIND)))
            do i=NBIN10(1)+NBIN10(2)+2,NBIN10(1)+NBIN10(2)+NBIN10(3)+1
               binbi(i)=srat_s*binbi(i-1)
            end do
         else
            LOG_ERROR("binmicrosetup_scale",*) "# of bins for ice is not 40, 20, 10", nbi
            call PRC_abort
            binbi(1) = minmass_s
            do i=2,nbi+1
               binbi(i)=srat_s*binbi(i-1)+sadd_s
            end do
         end if
!     write out the bin boundaries
         !if ( IsMaster ) then
         if ( debug ) then
           write(fid_alog,*) "bin boundaries of ice"
           write(fid_alog,*) "minmass_s: ",minmass_s
           write(fid_alog,*) "bin#,binbi,srat"
           write(fid_alog,'(I5,2ES15.6)') 1,binbi(1)
           do i=2,nbi+1
              write(fid_alog,'(I5,2ES15.6)') i,binbi(i),binbi(i)/binbi(i-1)
           end do
        end if

!
!     define split bin
!
         do i=1,nbr+1
!!!           if((binbr(i)/coef3)**(1.0/3.0).gt.25.0e-4) then
           if((binbr(i)/coef3)**(1.0/3.0).gt.20.0e-4) then
!           if((binbr(i)/coef3)**(1.0/3.0).gt.14.0e-4) then
             isplit_bin_liq=i
             exit
           endif
         end do
         do i=1,nbi+1
            if( binbi(i) > binbr(isplit_bin_liq) ) then
               isplit_bin_ice = i-1
               exit
            end if
         end do
         if ( IsMaster .and. debug ) then
           write(fid_alog,*) "Split bin for liquid and ice: ",isplit_bin_liq, isplit_bin_ice
         end if
```

**Supplementary — `com_amps.F90` lines 33–49** (declarations of `nbin_h`, tokens, `binbr`/`binbi` fixed size 101; `nbin_h` has no compile-time default — it is set only via namelist `PARAM_ATMOS_PHY_MP_AMPS_bin`):

```fortran
  integer :: level_comp, debug_level, coll_level, out_type, T_print_period, &
       token_c, token_r, token_s, token_a,&
       dtype_c, dtype_r, dtype_s, dtype_a, hbreak_c, hbreak_r, hbreak_s, &
       flagp_c, flagp_r, flagp_s, flagp_a, act_type,&
       nrmic, nrtime,nbin_h
  common/BMINOP/level_comp, debug_level, out_type, T_print_period, &
       token_c, token_r, token_s, token_a, &
       dtype_c, dtype_r, dtype_s, dtype_a(4), hbreak_c, hbreak_r, hbreak_s, &
       flagp_c, flagp_r, flagp_s, flagp_a, act_type,&
       nrmic, nrtime(17),nbin_h

  real(PS) :: srat_c, srat_r, srat_s, srat_a, sadd_c, sadd_r, sadd_s, sadd_a, sth_r,&
       fcon_c, minmass_c, minmass_r, minmass_s, minmass_a, binbi,binbr
  
  common/BMREOP/srat_c, srat_r, srat_s, srat_a, sadd_c, sadd_r, sadd_s, sadd_a, sth_r,&
       fcon_c,minmass_c,minmass_r,minmass_s, minmass_a, binbi(101),binbr(101)
```

## 2. `par_amps.F90` — FULL file
**File:** `/Users/jcanton/projects/scale_amps/contrib/AMPS/par_amps.F90`, lines 1–102

```fortran
Module par_amps
  implicit none
  public
  save
  ! parameters
  !---------------------------------------------
  ! ice particles
  ! array specfication for qipv
  ! NOTE:
  !      - aerosol mass components have to be in consecutive order.
  !
  !      - volume, lengths components have to be in consecutive order
  !        following the concentration.
  !
  !    These are redefined at the initialization based on the configuration.
  !
  !    total mass of ice particle
  integer :: imt_q=1
  !    rime, aggregate, crystal, melt water, frozen water (nucleation)
  integer :: imr_q=10,ima_q=15,imc_q=11,imw_q=12,imf_q=16
  !    total aerosol mass, soluble aerosol mass
  integer :: imat_q=13,imas_q=14
  !    number concentration
  integer :: icon_q=2
  !    circiumcsribing volume (dry), a axis length, c axis length, dendritic length
  integer :: ivcs_q=3,iacr_q=4,iccr_q=5,idcr_q=6
  !    center of gravity coordinates: a and c, and extra crystalline structure
  integer :: iag_q=7,icg_q=8,inex_q=9

  !
  ! array specification for g%MS%mass,g%MS%dmassdt,new_M
  integer :: imt=1
  integer :: imr=2,ima=3,imc=4,imw=8,imf=9, imat=5,imas=6,imai=7
  ! array specification for ratio_M
  integer :: imr_m=1,ima_m=2,imc_m=3,imw_m=7,imf_m=8, imat_m=4,imas_m=5,imai_m=6
  ! array specification for new_Q
  integer :: ivcs=1,iacr=2,iccr=3,idcr=4,iag=5,icg=6,inex=7

  !---------------------------------------------
  ! array specfication for qrpv
  integer :: rmt_q=1
  integer :: rmat_q=3,rmas_q=4
  integer :: rcon_q=2
  ! array specification for g%MS%mass,g%MS%dmassdt,new_M
  integer :: rmt=1
  integer :: rmat=2,rmas=3,rmai=4
  ! array specification for ratio_M
  integer :: rmat_m=1,rmas_m=2,rmai_m=3

  !---------------------------------------------
  ! array specfication for qapv
  integer :: amt_q=1
  integer :: acon_q=2
  integer :: ams_q=3
  ! array specification for g%MS%mass,g%MS%dmassdt,new_M
  integer :: amt=1
  integer :: ams=2,ami=3

  ! number of mass component, starding and ending indexes
  integer :: nvar_mcp_ice, i1_mcp_ice, i2_mcp_ice
  ! number of volumne component
  integer :: nvar_vcp_ice, i1_vcp_ice, i2_vcp_ice
  ! number of axis component, starting and ending indexes
  integer :: nvar_acp_ice, i1_acp_ice, i2_acp_ice
  ! number of concentration component, starting and ending indexes
  integer :: nvar_ccp_ice, i1_ccp_ice, i2_ccp_ice
  ! number of non-mass component
  integer :: nvar_nonmcp_ice
  ! number of property variables
  integer :: nvar_pv_ice

  ! number of mass component, starding and ending indexes
  integer :: nvar_mcp_liq, i1_mcp_liq, i2_mcp_liq
  ! number of concentration component, starting and ending indexes
  integer :: nvar_ccp_liq, i1_ccp_liq, i2_ccp_liq
  ! number of property variables
  integer :: nvar_pv_liq
  ! number of non-mass component
  integer :: nvar_nonmcp_liq

  ! number of mass component, starding and ending indexes
  integer :: nvar_mcp_aer, i1_mcp_aer, i2_mcp_aer
  ! number of property variables
  integer :: nvar_pv_aer
  ! number of concentration component, starting and ending indexes
  integer :: nvar_ccp_aer, i1_ccp_aer, i2_ccp_aer


  !---------------------------------------------------
  ! largest bin number that splits cloud and precipitating particles
  integer :: isplit_bin_liq=10
  integer :: isplit_bin_ice=0


  !---------------------------------------------------
  ! debug timing
  integer :: bug_time

  ! debug layer
  integer :: bug_layer

end Module par_amps
```

## 3. `maxdims.F90` — FULL file
**File:** `/Users/jcanton/projects/scale_amps/contrib/AMPS/maxdims.F90`, lines 1–89

```fortran
Module maxdims
  use acc_amps
  implicit none
  public
  ! these parameters are set in init_AMPS subroutine
  !
  !
  ! max # of grid boxes in z, x, and y for the atmospheric part, 
  ! and z for the soil model.
  integer :: n1mx=125,n2mx=50,n3mx=2,n4mx=1
  integer :: nxpmax,nypmax,nzpmax,nzgmax

  integer :: mxln
  integer :: mxlh, mxlv
  !
  ! max # of grid systems
  integer,parameter :: maxgrds=2
  integer,parameter :: ntrgrds=maxgrds-1
  !
  ! # of total grib boxes that are passed to shared memory for 
  ! microphysical calculation. The best value depends on the machine.
  !integer,parameter :: lbnmx=200 ! for CASE_1
  !integer,parameter :: lensmx=lbnmx*330+50
  !
  ! # of soil types
  integer,parameter :: nstyp=12
  !
  ! 6D microphysical variables
  !    see the explanation at MASLWIKI.
  ! ***** This if for bulk micro par *****
  ! max # of parameters, bins, and categories for liquid 6d variable (qrp)
!  integer,parameter :: nprmx=4,nbrmx=1,ncrmx=2
  ! max # of parameters, bins, and categories for ice 6d variable (qip)
!  integer,parameter :: npimx=4,nbimx=1,ncimx=5
  ! max # of parameters, bins, and categories for aerosol 6d variable (qap)
!  integer,parameter :: npamx=4,nbamx=1,ncamx=2
  ! ***** This if for CLR *****
  ! max # of parameters, bins, and categories for liquid 6d variable (qrp)
!  integer,parameter :: nprmx=4,nbrmx=1,ncrmx=2
  ! max # of parameters, bins, and categories for ice 6d variable (qip)
!  integer,parameter :: npimx=4,nbimx=1,ncimx=3
  ! max # of parameters, bins, and categories for aerosol 6d variable (qap)
!  integer,parameter :: npamx=3,nbamx=1,ncamx=5
  ! ***** This if for AMPS *****
  ! max # of parameters, bins, and categories for liquid 6d variable (qrp)
  integer :: nprmx=6,nbrmx=40,ncrmx=1
  ! max # of parameters, bins, and categories for ice 6d variable (qip)
  integer :: npimx=18,nbimx=20,ncimx=1
  ! max # of parameters, bins, and categories for aerosol 6d variable (qap)
  integer :: npamx=5,nbamx=1,ncamx=4
  !
  ! max # of variable initialization files
  integer,parameter :: maxhfils=50
  !
  ! # of variable used in KUO cumulus parameterization
  !integer,parameter :: nkp=2*n1mx
  !
  ! critical (or minimum possible) distance between two grid interfaces in z direction.
  real(RP),parameter :: dzcrit=1.0_RP
  ! ---------------------------------------------------------------------
      
  ! ---------------------------------------------------------------------
  ! AMPS parameters for memory allocation.
  !   Other parameters are defined in par_micro.f90
  !
  ! number of AMPS objects 
  ! if running bulk microphysics parameterization, set LMAX=1.
  ! Otherwise, LMAX should be equal to lbin.
!  integer,parameter :: LMAX=1
  integer :: LMAX
  !
  ! max number of bins for each categories
  integer :: mxnbinr,mxnbini,mxnbina
  !
  ! max number of bins and bin boundaries for liq and ice hydrometeors
  integer :: mxnbin,mxnbinb
  !
  ! max number of mass components and volume components for ice particles
  integer,parameter :: mxnmasscomp=8,mxnvol=2
  !
  ! max number of microphysical tendencies, axis lengths, and non-mass component
  ! variables for ice particles.
  integer,parameter :: mxntend=12,mxnaxis=5,mxnnonmc=7
  !
  ! max number of mass components for liquid particles
  integer,parameter :: mxnmasscomp_r=3
  !
  ! ---------------------------------------------------------------------
end Module maxdims
```

## 4. `read_AMPSTASK` — FULL subroutine
**File:** `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_utility.F90`, lines 1323–1423

```fortran
  subroutine read_AMPSTASK
    use scale_prc, only: &
       PRC_abort
    use com_amps
    implicit none
      !integer,parameter :: log_fid=101,ctl_fid=102
      integer :: LOG_FID, CTL_FID
      integer :: ierr

!     +++ input variables +++
      namelist /AMPS_param/&
         level_comp, &
         debug, &
         debug_level, &
         coll_level, &
         out_type, &
         T_print_period, &
         output_format, &
         token_a, &
         token_c, &
         token_r, &
         token_s, &
         dtype_c, &
         dtype_r, &
         dtype_s, &
         dtype_a, &
         hbreak_c, &
         hbreak_r, &
         hbreak_s, &
         flagp_c, &
         flagp_r, &
         flagp_s, &
         flagp_a, &
         act_type, &
         srat_c, &
         srat_r, &
         sth_r, &
         srat_s, &
         srat_a, &
         fcon_c, &
         sadd_c, &
         sadd_r, &
         sadd_s, &
         sadd_a, &
         minmass_c, &
         minmass_r, &
         minmass_s, &
         minmass_a, &
         M_aps, &
         M_api, &
         den_aps, &
         den_api, &
         nu_aps, &
         phi_aps, &
         ap_lnsig, &
         ap_mean, &
         N_ap_ini, &
         eps_ap, &
         ap_mean_cp, &
         ap_sig_cp, &
         APSNAME, &
         DRCETB, &
         DRAPTB, &
         DRSTTB, &
         micexfg, &
         n_step_cl, &
         n_step_vp, &
         ihabit_gm_random, &
         CCNMAX, &
         frac_dust, &
         nucleation_halflife, &
         CRIC_RN_IMM

        CTL_FID = IO_get_available_fid()
        open(CTL_FID, &
          file='AMPSTASK.F', &
          form='formatted',   &
          status='old',       &
          iostat=ierr)
        if(ierr/=0) then
           LOG_ERROR("read_AMPSTASK",*) 'Cannot open AMPS PARAMETER file, please check.'
           call PRC_abort
        end if
        rewind(CTL_FID)
        read(CTL_FID,nml=AMPS_param)
        close(CTL_FID)

        if ( IsMaster .and. debug ) then
          LOG_FID = IO_get_available_fid()
          open(LOG_FID, &
            file='AMPSTASK'//'_'//string_fid_amps//'.log', &
            form='formatted',   &
            status='replace',       &
            iostat=ierr)

          write(LOG_FID,nml=AMPS_param)
          write(LOG_FID,*) "THIS IS STANDARD CODE"
          close(LOG_FID)
        end if

      end subroutine read_AMPSTASK
```

## 5. `AMPSTASK.F` — FULL file
**File:** `/Users/jcanton/projects/scale_amps/scale-rm/test/case/cloudlab/AMPSTASK.F`, lines 1–293

```fortran
&AMPS_param
! *********************************************************************************
! integer variable of options
! *********************************************************************************
! nbin = 0 means no construction of the group.
!
! 1. level of complexity
! level_com:1
!          :2 a, c axis lengths and dendric length  predicted
!          :3 rosette and irregular length added to 2.
!          :4 aerosol mass components are predicted to 2.
!          :5 aerosol mass components are predicted to 3.
!
! 1. level of complexity
level_comp = 7
! 2. debugging level
debug_level = 1
! 3. what model needs
!    1 -> only tendency is updated.
!    2 -> tendency and variable are updated.
out_type = 2
!
! 4. time period for printing (s)
!T_print_period = 1800
!T_print_period = 1200
T_print_period = 3600000

! 5. output file format
!      defaut : text format, 'binary' : binary format (stream)
output_format = 'binary',
!- Cloud_drop group -
! 5. token (type of hydrometeor)
token_c = 11
! 7. dis_type
dtype_c = 0
! 8. h_break_method
hbreak_c = 0
! 9. flag for prediction (concentration fixed?)
!    1 -> concentration is fixed
!    2 -> concentration is predicted.
flagp_c = 2
!
! - rain group -
! 10. token (type of hydrometeor)
token_r = 1
! 12. dis_type
dtype_r = 1
! 13. type of hydrodynamic breakup
hbreak_r = 2
! 14. flag for prediction (concentration fixed?)
flagp_r = 2
!
! - Solid hydrometeor group -
! 15. token (type of hydrometeor)
token_s = 2
! 17. dis_type
dtype_s = 1
! 18. type of hydrodynamic breakup
hbreak_s = 2
! 19. flag for prediction (concentration fixed?)
flagp_s = 2
!
! - aerosol group -
! 20. token (type of hydrometeor)
token_a = 3
! 21. dis_type for category 1, and 2
dtype_a = 3,4,3,3
! 23. flag for prediction (concentration fixed?)
! 0: no prediction in micro
! 1: predict concentration and diagnose mass content (standard deviation fixed)
! 2: predict concentration and mass content (standard deviation fixed)
! 3: predict concentration and mass content (mean radius fixed)
! -1: initialize aerosol groups, but predict only CCN category and pass it to dynamic model (standard deviation fixed)
! -2: initialize aerosol groups, but predict only IN category and pass it to dynamic model (standard deviation fixed)
! -3: initialize aerosol groups, but do not pass any change of aerosols to dynamic model (standard deviation fixed)
! -4: initialize aerosol groups, but predict only CCN category and pass it to dynamic model (mean radius fixed)
! -5: initialize aerosol groups, but predict only IN category and pass it to dynamic model (mean radius fixed)
! -6: initialize aerosol groups, but do not pass any change of aerosols to dynamic model (mean radius fixed)
!flagp_a = 1
flagp_a = 2
!flagp_a = -3
!
! 24. habit growth determination
!     0: max frequency, 1: random number method
ihabit_gm_random = 1
!
! *********************************************************************************
! real variable of options
! *********************************************************************************
! - Cloud_drop group -
! 1. size_ratio
srat_c = 0.0
! 2. size_addition
sadd_c = 0.0
! 3. minimum_mass
minmass_c = 0.0
! 4. fixed concentration
fcon_c = 25.0
!
! - Rain group -
! 5. size_ratio
!    for n=20
srat_r = 2.540068909
! 6. size_addition
sadd_r = 0.0
! 7. minimum_mass
minmass_r = 4.188790E-09
!    for n=20
sth_r = 1.02
!
! - Solid hydrometeor group -
! 8. size_ratio
! this is for 83 bins from 4.19e-12 to 10.0 g
!!!srat_s = 1.4142135623731
! this is for 20 bins from 4.19e-12 to 10.0 g
srat_s = 4.1581061E+00
! 9. size_addition
sadd_s = 0.0
! 10. minimum_mass
minmass_s = 4.18879020478639E-12
!
! - aerosol group -
! 11. size_ratio
srat_a = 1.0e-6
! 12. size_addition
sadd_a = 0.0
! 13. minimum_mass
minmass_a = 1.0E-21
!
! for physical property of CCN, and IN
!
! 14. molecular weights of soluble part
! NaCl
! M_aps = 58.45,58.45,58.45,58.45
! (NH4)2SO4
! M_aps = 132.14052,132.14052,132.14052,132.14052
! NH4HSO4
M_aps = 115.11,115.11,115.11,115.11,
! 15. molecular weights of insoluble part
! silver iodide (AgI)
M_api = 234.77,234.77,234.77,234.77
! Kaolinite (Al2Si2O5(OH)4)
! M_api = 258.16,258.16,258.16,258.16
! 16. bulk density of soluble material for CCN
! NaCl
! den_aps = 2.16,2.16,2.16,2.16
! (NH4)2SO4
! den_aps = 1.77,1.77,1.77,1.77
! NH4HSO4
den_aps = 1.79,1.79,1.79,1.79,
! 17. bulk density of insoluble material for CCN
! silver iodide (AgI)
den_api = 5.683,5.683,5.683,5.683
! Kaolinite (Al2Si2O5(OH)4)
!den_api = 2.6,2.6,2.6,2.6
! 18. number of ions in dessociated solutes
! NaCl
! nu_aps = 2.0, 2.0, 2.0, 2.0,
! (NH4)2SO4
! nu_aps = 3.0,3.0,3.0,3.0,
! NH4HSO4
nu_aps = 2.0,2.0,2.0,2.0,
! 19. molal coefficient (average)
phi_aps = 0.75, 0.75, 0.75, 0.75,
! 20. mass fraction of soluble material to total mass
! eps_ap = 1.0,0.0,0.2,0.7
eps_ap = 1.0,0.0,0.05,1.0
! 21. log of standard deviation (in natural log of radius in micro m)
! from Jaenicke (1993), marine type, mode II,I,III: log(sig)=0.210, 0.657, 0.396 in diameter
!    log(sig) in radius is exactly the same as those in diameter.
! ap_lnsig = 0.210, 1.0, 0.657, 0.396,
! MPACE (sig=2.04, 2.5)
!ap_lnsig = 0.712949807856125, 1.0, 0.916290731874155, 1.0,
! for W2 in Shipway and Hill (2012)
!ap_lnsig = 0.405465108108164, 1.0, 0.405465108108164, 0.405465108108164,
! for SHEBA
!ap_lnsig = 2.04, 1.0, 2.5, 1.0
ap_lnsig = 0.712949807856125,0.0,0.916290731874155,0.916290731874155
! 22. geometrical mean radius (in cm)
!  from Jaenicke marine type
! ap_mean = 0.133e-4,1.0e-4,0.004e-4,0.29e-4
! MPACE
! ap_mean = 0.052e-4, 1.0e-4, 1.3e-4, 0.001e-4,
! SHEBA
ap_mean = 0.052e-4, 1.0e-4, 1.3e-4, 1.3e-4
! for W2 in Shipway and Hill (2012)
!ap_mean = 0.1e-4, 1.0e-4, 0.1e-4, 0.1e-4,
!
! 23. geometrical mean contact parameter (in degree), normal distriubtion
ap_mean_cp = 132.0, 15.5, 132.0, 132.0,
!
! 24. standard deviation of contact parameter (in degree), normal distribution
!
ap_sig_cp = 20.0, 1.4, 20.0, 20.0,
!
! 25. Initial number concentraion of aerosol particles
!N_ap_ini = 350.0, 0.00015, 1.8, 1.8,
N_ap_ini = 317.0, 0.0, 0.0, 0.0,
!
! * character variable of options
! 1. DRCETB
! - directory where the collision coefficiency files are stored.
DRCETB ='/cluster/scratch/congchia/scale_amps/scale-rm/test/case/cloudlab/AMPS_DATA/collision_data'
! 2. IFMGTB
! - parameter or database file for mass-geometric relation of ice crystals.
! DRMGTB = 'prmRHI.v1.dat',
! 3. DRAPTB
! - directory where the aerosol activation file is stored.
DRAPTB = '/cluster/scratch/congchia/scale_amps/scale-rm/test/case/cloudlab/AMPS_DATA/apact'
! 4. DRSTTB
! - directory where the statistics lookup table is stored.
DRSTTB = '/cluster/scratch/congchia/scale_amps/scale-rm/test/case/cloudlab/AMPS_DATA/statpack'
!
! 5. Name of chemicals for aerosols
! APSNAME = 'NACL2','NACL2','NACL2'
! APSNAME = 'NH42SO4','NH42SO4','NH42SO4','NH42SO4'
APSNAME = 'NH4HSO4','NH4HSO4','NH4HSO4','NH4HSO4'
!
! 6. model type of CCN activation
! 1: implicit prediction of vapor (default)
! 2: water saturation adjustment
act_type = 1
!
! sensitivity test
!
coll_level = 1
!
! max CCN
!
CCNMAX = 400.0D0
!
! dust fraction
!
frac_dust = 0.01D0
!
! ice nucleation half life parameter
!
nucleation_halflife = 0.0001D0
!
! minimum dust radius (cm)
!
CRIC_RN_IMM = 0.25D-4
!
! simulation under fixed Relative Humidity and Temperature
!
! relative humidity over water (%)
!FRHW = 100.0
! relative humidity over ice (ratio)
!FRHI = -999.9
! temperature in Celsius
!FTMP = -5.0
!
! number of time steps relative to one dynamic step
!  collision processes
n_step_cl = 1
! number of time steps relative to collision processes
!  vapor deposition processes
n_step_vp = 10
!
! execution flag for microphysical processes
! on -> 1, off -> 0
!   some processes, other schemes may be chosen with this flag.
!
! index
! 1: printing
! 2: liq-liq collsion-coalescence process
! 3: ice-ice collsion-coalescence process
! 4: liq-ice collsion-coalescence process
! 5: update surface temperature
! 6: vapor deposition on liq
! 7: vapor deposition on ice
! 8: melting shedding process
! 9: hydrodynamic breakup of ice particles
! 10: ice nucleation process
! 11: hydrodynamic breakup of liquid hydrometeors
! 12: autoconversion of cloud droplet bin
! 13: ice nucleation; depositional nucleation
! 14: ice nucleation; contact freezing
! 15: ice nucleation; splinter nulcreation
! 16: ice nucleation; homogenious freezing
! 17: ice nucleation; immersion freezing (heterogenious freezing)
!         0                 1
!         1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
!micexfg = 1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1
!micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,1
micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,1,0,0,0,0,1,0,1

!micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,1,1,0,1,1,1,1,1
!micexfg = 1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,1
!micexfg = 1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,2,0,1,1,1
!micexfg = 1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1
!
/
```
(Note: lines 277–281 of the file label indices 13–17 as shown above — the file's own comment lists "16: ice nucleation; immersion freezing" and "17: ... homogenious" in the order quoted; quoted exactly as in the file.)

## 6a. `scale_atmos_phy_mp_amps.F90` — module-level config defaults
**File:** `/Users/jcanton/projects/scale_amps/scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90`, lines 121–179

```fortran
  integer, parameter :: k1eta = 2

  integer, public :: QA, I_QV, I_QL, I_QW, I_QI, I_QPPVL, I_QPPVI, I_QPPVA
  integer :: ivis = 1
  !integer, dimension(max_nmoments_liq)  :: I_scl2ship_l
  !integer, dimension(max_nmoments_ice)  :: I_scl2ship_i
  !integer, dimension(max_nmoments_aero) :: I_scl2ship_a

  ! internal time step
  integer :: TIME_AMPS

  logical  :: amps_debug         = .false.
  logical  :: amps_ignore        = .false.
  logical  :: l_restart          = .false.
  logical  :: l_fix_aerosols     = .true.
  logical  :: l_sediment         = .true.
  logical  :: l_no_ice_heat      = .false.
  logical  :: l_fill_aerosols    = .false.
  logical  :: l_bin_shift        = .false.
  logical  :: l_axis_limit       = .true.
  integer  :: ini_aerosol_prf    = 3
  integer  :: l_gaxis_version    = 1
  integer  :: l_aadv_version     = 2
  integer  :: l_reff_version     = 2
  integer  :: iadvv              = 1
  logical  :: fix_aerosol_type(4) = (/.true., .true., .true., .true./)

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

  integer  :: nx, ny, nz, nzh
  integer, parameter :: max_nmoments_liq=4, max_nmoments_ice=16, max_nmoments_aero=3
  integer, parameter :: nspecies=2, max_char_len=200
  integer, parameter :: naerosol=4 ! number of aerosol species
  integer, parameter  :: num_aero_bins(naerosol) = 1
  integer, parameter  :: num_aero_moments(naerosol) = 3
  integer, parameter  :: num_h_moments(nspecies) = (/  &
              4 & ! liquid
              ,16 & ! ice
              /)
  integer  :: num_h_bins(nspecies)= (/ &
  ! for fine
  !         80 & ! liquid
  !        ,83 &  ! ice
  ! for midium
  !         60 & ! liquid
  !        ,40 &  ! ice
  ! for coarse (default)
            40 & ! liquid
           ,20 &  ! ice
           /)
```

## 6b. Moment-name blocks
**Same file, lines 181–292**

```fortran
  character(len=9) :: mom_names_ice(max_nmoments_ice)= &
       (/ 'ice_imass'   &                        ! ice total mass
       ,  'rim_imass'   &                        ! rime mass
       ,  'agg_imass'   &                        ! aggragate mass
       ,  'cry_imass'   &                        ! crystal mass
       ,  'mtw_imass'   &                        ! melt water mass
       ,  'frw_imass'   &                        ! frozen water mass (nucleation)
       ,  'tae_imass'   &                        ! total aerosol mass
       ,  'sat_imass'   &                        ! soluble aerosol mass
       ,  'con_imass'   &                        ! number concentration
       ,  'vol_imass'   &                        ! circumscribing volume
       ,  'a_axis   '   &                        ! a-axis length
       ,  'c_axis   '   &                        ! c-axis length
       ,  'd_axis   '   &                        ! d-axis length
       ,  'ag_axis  '   &                        ! center of gravity a (polycrystals)
       ,  'cg_axis  '   &                        ! center of gravity c (polycrystals)
       ,  'ex_cry   '   &                        ! extra crystalline structure
       /)

  character(len=34) :: mom_descriptions_ice(max_nmoments_ice)= &
       (/ 'ice total mass                    '    &
       ,  'rime mass                         '    &
       ,  'aggragate mass                    '    &
       ,  'crystal mass                      '    &
       ,  'melt water mass                   '    &
       ,  'frozen water mass (nucleation)    '    &
       ,  'total aerosol mass                '    &
       ,  'soluble aerosol mass              '    &
       ,  'number concentration              '    &
       ,  'circumscribing volume             '    &
       ,  'a-axis length                     '    &
       ,  'c-axis length                     '    &
       ,  'd-axis length                     '    &
       ,  'center of gravity a (polycrystals)'    &
       ,  'center of gravity c (polycrystals)'    &
       ,  'extra crystalline structure       '    &
       /)

  character(len=6) :: mom_units_ice(max_nmoments_ice)= &
       (/ 'kg/kg '   &                        ! ice total mass
       ,  'kg/kg '   &                        ! rime mass
       ,  'kg/kg '   &                        ! aggregate mass
       ,  'kg/kg '   &                        ! crystal mass
       ,  'kg/kg '   &                        ! melt water mass
       ,  'kg/kg '   &                        ! frozen water mass (nucleation)
       ,  'kg/kg '   &                        ! total aerosol mass
       ,  'kg/kg '   &                        ! soluble aerosol mass
       ,  '1/kg  '   &                        ! number concentration
       ,  'cm3/kg'   &                        ! circumscribing volume
       ,  '1/kg  '   &                        ! a-axis length
       ,  '1/kg  '   &                        ! c-axis length
       ,  '1/kg  '   &                        ! d-axis length
       ,  '1/kg  '   &                        ! center of gravity a (polycrystals)
       ,  '1/kg  '   &                        ! center of gravity c (polycrystals)
       ,  '1/kg  '   &                        ! extra crystalline structure
       /)

  character(len=9) :: mom_names_liq(max_nmoments_liq)= &
       (/ 'liq_lmass'   &                        ! liquid total mass
       ,  'tae_lmass'   &                        ! total aerosol mass
       ,  'sae_lmass'   &                        ! soluble aerosol mass
       ,  'con_lmass'   &                        ! number concentration
       /)

  character(len=20) :: mom_descriptions_liq(max_nmoments_liq)= &
       (/ 'liquid total mass   '   &
       ,  'total aerosol mass  '   &
       ,  'soluble aerosol mass'   &
       ,  'number concentration'   &
       /)

  character(len=5) :: mom_units_liq(max_nmoments_liq)= &
       (/ 'kg/kg'   &                        ! liquid total mass, rho_r/rho_total
       ,  'kg/kg'   &                        ! total aerosol mass
       ,  'kg/kg'   &                        ! soluble aerosol mass
       ,  '1/kg '   &                        ! number concentration
       /)

  character(len=4) :: aero_names(naerosol)= &
       (/ 'ccn1'   &
       ,  'in1 '   &
       ,  'ccn2'   &
       ,  'ccn3'   &
       /)

  character(len=9) :: mom_names_aero(max_nmoments_aero)= &
       (/ 'ae_amass '   &
       ,  'con_amass'   &
       ,  'sae_amass'   &
       /)

  character(len=20) :: mom_descriptions_aero(max_nmoments_aero)= &
       (/ 'total aerosol mass  '  &
       ,  'number concentration'  &
       ,  'soluble aerosol mass'  &
       /)

  character(len=5) :: mom_units_aero(max_nmoments_aero)= &
      (/  'kg/kg'   &
       ,  '1/kg '   &
       ,  'kg/kg'   &
       /)

  character(len=32),dimension(20) :: proname_ifc,proname_mtd

  ! for liquid category
  integer :: split_bins = 21
  ! for ice category
  integer :: ibin_is
  integer :: ibin_sg
  integer :: ibin_gh
```

**Supplementary — tracer_setup namelist that sets `nbin_h` (same file, lines 357–381; `nbin_h` lives in `com_amps` and has no compile-time default):**

```fortran
    namelist / PARAM_ATMOS_PHY_MP_AMPS_bin / &
       num_h_bins,         & ! 40 or 80 bins for liquid, ONLY 20 bins for ice
       nbin_h,             & ! number of bins for haze particles, e.x. 20 for 40 liq. bins
       iadvv,              & ! 1 for Euler, 2 for PPM sedimentation scheme
       ini_aerosol_prf,    & ! 1 for SHEBA, 2 for MPACE, 3 for general use (1 and 2 OVERRIDE aerosol settings in AMPSTASK.F)
       l_restart,          & ! whether restart to fill aerosols, this should be replaced by SCALE restart namelist!!!!!!!
       l_fix_aerosols,     & ! invariant aerosols throughout integration
       l_sediment,         & ! sediment on or off
       l_no_ice_heat,      & ! ignore latent heat release due to ice processes
       l_fill_aerosols,    & ! fill aerosols in cloud-free region or not
       l_bin_shift,        & ! whether bin shift is performed after advection, default is false for testing phase, it is still under checking
       l_gaxis_version,    & ! axis definition, 1 (a, c, d, ag, cg), 2 (a, c/a, d/a, ag, cg/ag), 3 (a, c, d, ag/a, cg/a), default is 1
       l_aadv_version,     & ! 1 for all PPVs advected in the same fasion, 2 for modified advection for non-mass PPVs
       l_axis_limit,       & ! whether center of gravity axis (ag and cg) limit=1 is applied, default is true
       l_reff_version,     & ! using (1) maximum dimension or (2) equivalent spherical radius for calculating radiation
       fix_aerosol_type,   & ! types of aerosols that are fixed if l_fix_aerosols is true
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

## 6c. `make_Cloud_Micro_cnfg` — DUMMY-argument declaration section
**File:** `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Cloud_Micro.F90`, lines 229–322 (signature + config-scalar declarations with comments; remaining args, lines 323–358, are LUT array passthroughs)

```fortran
  subroutine make_Cloud_Micro_cnfg( &
       CM &
      ,dt_step, n_step_cl, n_step_vp, c_time  &
      ,NRTYPE,NRBIN, NRCAT  &
      ,NSTYPE,NSBIN, NSCAT  &
      ,NATYPE,NABIN, NACAT  &
      ,estbar,esitbar  &
      ,level_comp,debug_level,coll_level_IN,out_type,T_print_period,output_format &
      ,token_r,token_s,token_a,dtype_c,dtype_r,dtype_s,dtype_a &
      ,hbreak_c,hbreak_r,hbreak_s,flagp_c,flagp_r,flagp_s,flagp_a &
      ,act_type,micexfg &
      ,fcon_c,coef_ap,eps_ap &
      ,srat_c,srat_r,srat_s,srat_a &
      ,sadd_c,sadd_r,sadd_s,sadd_a &
      ,minmass_r,minmass_s,minmass_a &
      ,den_apt,den_aps,den_api &
      ,binbr,binbi &
      ,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_fd,bu_tmass &
      ,APSNAME,nu_aps,phi_aps,m_aps,ap_lnsig,ap_mean &
      ,ap_sig_cp,ap_mean_cp,cdf_cp_0,cdf_cp_180m0 &
      ,nr_drpdrp,nc_drpdrp,xs_drpdrp,dx_drpdrp,ys_drpdrp,dy_drpdrp,drpdrp &
      ,nr_hexdrp,nc_hexdrp,xs_hexdrp,dx_hexdrp,ys_hexdrp,dy_hexdrp,hexdrp &
      ,nr_bbcdrp,nc_bbcdrp,xs_bbcdrp,dx_bbcdrp,ys_bbcdrp,dy_bbcdrp,bbcdrp &
      ,nr_coldrp,nc_coldrp,xs_coldrp,dx_coldrp,ys_coldrp,dy_coldrp,coldrp &
      ,nr_gp1drp,nc_gp1drp,xs_gp1drp,dx_gp1drp,ys_gp1drp,dy_gp1drp,gp1drp &
      ,nr_gp4drp,nc_gp4drp,xs_gp4drp,dx_gp4drp,ys_gp4drp,dy_gp4drp,gp4drp &
      ,nr_gp8drp,nc_gp8drp,xs_gp8drp,dx_gp8drp,ys_gp8drp,dy_gp8drp,gp8drp &
      ,nok_igp,x_igp,a_igp,b_igp &
      ,n_osm_nh42so4,xs_osm_nh42so4,dx_osm_nh42so4,y_osm_nh42so4 &
      ,n_osm_sodchl,xs_osm_sodchl,dx_osm_sodchl,y_osm_sodchl &
      ,n_snrml,xs_snrml,dx_snrml,y_snrml &
      ,n_isnrml,xs_isnrml,dx_isnrml,y_isnrml &
      ,ihabit_gm_random &
      ,CCNMAX,CRIC_RN_IMM,frac_dust,nucleation_halflife &
      ,lbin)
!!c       lbin) result (CM)
!!c       binbr,binbi) result (CM)
    use class_AirGroup, only: &
       make_AirGroup_alloc
    use class_Group, only: &
       make_Group, &
       make_col_lut

    type (Cloud_Micro),intent(inout)    :: CM
!tmp    type (Cloud_Micro)    :: CM
!tmp    COMMON /ICMPRIV/CM
!tmp!$omp threadprivate(/ICMPRIV/)

    ! model time step in second
    real(DS), intent(in) :: dt_step
    real(PS) :: time_step
    ! number of time steps relative to one dynamic time step
    integer, intent(in) :: n_step_cl, n_step_vp
    ! current model time in second
    real(DS), intent(in) :: c_time
    !
    integer, intent(in) :: NRBIN, NSBIN, NABIN
    integer, intent(in) :: NRTYPE, NSTYPE, NATYPE
    integer, intent(in) :: NRCAT, NSCAT, NACAT
    integer, intent(in) :: lbin

    ! random generaion: 1, max frequency: 0
    integer, intent(in) :: ihabit_gm_random

    real(DS), intent(in) :: estbar(150), esitbar(111)

    integer, intent(in) :: coll_level_IN

    integer :: level_comp,debug_level,out_type,T_print_period,&
       token_r,token_s,token_a,dtype_c,dtype_r,dtype_s,dtype_a(*),&
       hbreak_c,hbreak_r,hbreak_s,flagp_c,flagp_r,flagp_s,flagp_a,&
       act_type,micexfg(*)
    real(PS) :: srat_c,srat_r,srat_s,srat_a,sadd_c,sadd_r,sadd_s,sadd_a,&
       fcon_c,minmass_r,minmass_s,minmass_a,&
       den_apt(*),den_aps(*),den_api(*),&
       binbr(*),binbi(*),&
       coef_ap(*),eps_ap(*)
    real(PS), intent(in) :: CCNMAX, CRIC_RN_IMM, frac_dust, nucleation_halflife

    ! output format
    character(len=16), intent(in) :: output_format

    ! bin boundaries for aerosols have to be introduced here.
    real(PS) :: binba(82)

    ! collisional-breakup variables
    integer, intent(in) :: imin_bk,imax_bk,jmin_bk,jmax_bk
    real(PS) :: bu_fd(2,*),bu_tmass(*)
    ! chemical notation for aerosols
    character (len=16) :: APSNAME(*)
    ! aerosol variables
    real(PS),dimension(*) :: nu_aps,phi_aps,m_aps,ap_lnsig,ap_mean &
           ,ap_sig_cp,ap_mean_cp,cdf_cp_0,cdf_cp_180m0
```

## 7. `make_Group` — component counts per token (incl. a_b/b_b)
**File:** `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Group.F90`, lines 201–339 (subroutine ends line 577; token 12 rain-category block continues at line 359)

```fortran
  subroutine make_Group(g, dL, ddt, tok, n_bin, dab, dbb, dmbinb, dtype,&
       dden,dden_as,dden_ai,deps_map,dbinb)
    use class_Mass_Bin, only: &
       make_mass_bin

    integer, intent(in) :: dL, tok, n_bin, dtype
!tmp    integer, optional,intent(in) :: dNtp
    real(PS),intent(in)    :: ddt, dab, dbb, dmbinb,dden,dden_as,dden_ai,deps_map
!tmp    real(PS),optional,dimension(n_bin+1) :: dbinb
    real(PS),dimension(*) :: dbinb
    type (Group)  :: g

    ! The boundary of the bin binb_i is given by
    ! binb_i = a_b * binb_{i-1} + b_b
    !real(PS)            :: a_b, b_b
    integer :: i,j


    allocate( g%MS(mxnbin,LMAX) )
    allocate( g%IS(mxnbin,LMAX) )
    allocate( g%mark_cm(LMAX) )
    allocate( g%mark_er(LMAX) )
    allocate( g%binb(mxnbinb) )

    ! definition of total number of grids point with microphysics in it.
    g%L = dL

    ! definition of time step
    g%dt = ddt

    ! definition of token
    g%token = tok

    ! definition of parameters defining bin property
    g%N_BIN = n_bin
    g%N_binb = n_bin + 1

    g%a_b = dab
    g%b_b = dbb
    g%mbinb = dmbinb


    ! +++ definition of concentration on the boundary of bins +++
    !     allocate memory to type of distribution +++
!!c    allocate( g%con_onbinb(g%N_binb), stat = var_Status)
!!c    if(var_Status /= 0 ) stop "Memory not available for con_onbinb &
!!c         in class_group"
!!c    if( present(conb)) then; g%con_onbinb = conb
!!c    else;                  g%con_onbinb   = 0.0_PS
!!c    end if


    ! +++ definition of distribution
    g%org_dtype = dtype

    ! +++ memory allocation +++
!!c    call allocate_group(g)


!!c    write(*,*) "make_group",dL,ddt,tok,n_bin,dab,dbb,dmbinb,dtype
!!c    write(*,*) "makegrou2",dden,dden_as,dden_ai,deps_map
!!c    write(*,*) "makegroup3",dbinb

    ! +++ calculation of bin boundaries +++
       do i = 1, g%N_binb
          g%binb(i) = dbinb(i)
!!c          write(*,*) "gbinb",i,g%binb(i)
       end do
!tmp    if( present(dbinb)) then
!tmp       do i = 1, g%N_binb
!tmp          g%binb(i) = dbinb(i)
!tmp!!c          write(*,*) "gbinb",i,g%binb(i)
!tmp       end do
!tmp    else
!tmp       g%binb(1) = g%mbinb
!tmp       do i = 2, g%N_binb
!tmp          g%binb(i) = g%a_b*g%binb(i-1) + g%b_b
!tmp       end do
!tmp    end if
    ! +++ definition of number of process for tendency +++
       g%N_tendpros = 12
!tmp    if( present(dNtp)) then
!tmp       g%N_tendpros = dNtp
!tmp    else
!tmp       g%N_tendpros = 12
!tmp    end if

    if( g%token == 1 ) then
      ! +++ construct the liquid objects +++
      g%N_axis = 0
      g%N_rimemass = 0
      g%N_aggmass = 0
      g%N_apmass=3
      g%N_meltmass=0
      g%N_masscom=g%N_apmass
      g%N_vol = 0
      g%N_nonmass = g%N_vol

      call initialize

    else if( g%token == 2 ) then
      ! +++ construct the ice objects +++
      g%N_axis = 5
      g%N_rimemass = 1
      g%N_aggmass = 1
      g%N_apmass=3
      g%N_meltmass=1
      g%N_frozenmass=1
      g%N_masscom = g%N_rimemass+g%N_aggmass+1+g%N_apmass+g%N_meltmass+g%N_frozenmass
      g%N_vol = 0
      g%N_nonmass = 1 + g%N_axis + 1 + g%N_vol

      call initialize

      call initialize_shape

    else if( g%token == 3 ) then
      ! +++ construct the aerosol objects +++
      g%N_axis = 0
      g%N_rimemass = 0
      g%N_aggmass = 0
      g%N_apmass=2
      g%N_meltmass=0
      g%N_masscom=g%N_apmass
      g%N_vol = 0
      g%N_nonmass = g%N_vol

      call initialize

    else if( g%token == 11 ) then
       ! +++ construct cloud droplets category +++
       g%N_axis = 0
       g%N_rimemass = 0
       g%N_aggmass = 0
       g%N_apmass=2
       g%N_meltmass=0
       g%N_masscom=g%N_apmass
       g%N_vol = 0
       g%N_nonmass = 0
```

(For completeness: token 12 "rain category" at lines 359–368 sets the same counts as token 11: `N_axis=0, N_rimemass=0, N_aggmass=0, N_apmass=2, N_meltmass=0, N_masscom=N_apmass, N_vol=0, N_nonmass=0`, with a duplicated `g%N_apmass=2` at line 369. Group type declaration `real(PS) :: a_b,b_b` is at class_Group.F90:103.)