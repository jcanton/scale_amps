I have everything needed. Here is the port spec.

---

# Port Spec — `cal_breakfragment` Cluster-1 (breakup math), standalone extraction for M2b

**Repo:** `/Users/jcanton/projects/scale_amps` · **Source files:** `contrib/AMPS/mod_amps_lib.F90`, `contrib/AMPS/mod_amps_core.F90`, `contrib/AMPS/mod_amps_utility.F90`, `contrib/AMPS/class_Group.F90`, `contrib/AMPS/acc_amps.F90`

Scope of this doc: the **pure breakup math (Cluster 1)** identified in `docs/superpowers/facts/m2/breakup-tables-and-icon4py-m2.md` PART A, extracted VERBATIM, plus the **minimal `%len`/`%vtm`** substitute that lets M2b skip the monolithic `diag_pq` (Cluster 2). The `cal_breakfragment` driver itself is already in the fact file (§A.1) — re-confirmed below against source.

---

## 0. Re-confirmation of the driver and call tree

`cal_breakfragment` at `contrib/AMPS/mod_amps_lib.F90:1831` — the verbatim body in fact-file §A.1 (lines 1831–2017, 187 lines) **matches source exactly**. The `use mod_amps_core, only: cal_Coalescence_Efficiency, cal_breakup_dis_LL, diag_pq` import and the two inner calls

```fortran
call cal_Coalescence_Efficiency(liquid,i,liquid,j,1,steady%TV(1),E_coal, D_L,D_S,S_T,S_C,DS_S,CKE)   ! mod_amps_lib.F90:169-170
if(CKE<=1.0e-20) cycle                                                                                  ! :172
call cal_breakup_dis_LL(liquid,i,j,imin_bk,jmin_bk,bu_tmass,bu_fd, D_L,D_S,S_T,S_C,CKE)                 ! :175-176
```

are confirmed. `diag_pq` (line 128-129) is the **only** Cluster-2 dependency in the driver, and it exists solely to populate `liquid%MS(:,1)%len` and `%vtm` before the pair loop. Section 8 below replaces it.

**Cluster-1 routines to port (all VERBATIM below):**

| # | Routine | file:line | lines |
|---|---|---|---|
| 1 | `cal_breakfragment` (driver) | mod_amps_lib.F90:1831–2017 | 187 (in fact file) |
| 2 | `cal_breakup_dis_LL` | mod_amps_core.F90:12019–12447 | 429 |
| 3 | `cal_Coalescence_Efficiency` (token==1 branch only) | mod_amps_core.F90:11796–12017 | 222 |
| 4 | `cal_sig_sf` | mod_amps_core.F90:27254–27302 | 49 |
| 4 | `cal_Hmusig` | mod_amps_core.F90:27304–27354 | 51 |
| 5 | `zbrent` (+ internal `func`) | mod_amps_utility.F90:12778–12907 | 130 |
| 6 | `getznorm2` | mod_amps_utility.F90:369–375 | 7 |
| 6 | `cdfnor` | mod_amps_utility.F90:10033–10230 | 198 |
| 6 | `cumnor` | mod_amps_utility.F90:10231–10435 | 205 |
| 6 | `d_swap` | mod_amps_utility.F90:10436–10466 | 31 |
| 6 | `dinvnr` | mod_amps_utility.F90:10467–10548 | 82 |
| 6 | `stvaln` | mod_amps_utility.F90:11698–11764 | 67 |
| 6 | `eval_pol` | mod_amps_utility.F90:11765–11808 | 44 |

**Correction to the fact-file A.2 inventory:** the "`devlpl`" named in the fact file is not present in this tree; the polynomial evaluator behind `stvaln` is **`eval_pol`** (`mod_amps_utility.F90:11765`). Also note `getznorm2` calls `cdfnor(which=1,...)`, which reaches **only `cumnor`** — `dinvnr`/`stvaln`/`eval_pol` are the `which=2/3/4` (inverse-CDF) path and are **not on the breakup hot path**. They are quoted in full as requested but M2b can defer them (they are needed only if the port ever needs the inverse normal CDF; `zbrent` does the root-finding here instead, using only the forward `getznorm2`).

---

## 1. Constants (`contrib/AMPS/acc_amps.F90`)

```fortran
real(MP_KIND), parameter :: PI = 3.141592653589793238462643_PS                 ! acc_amps.F90:23
real(MP_KIND), parameter :: coefpi6=0.523598776_MP_KIND, &
                            coefsq2p=2.506628274631_MP_KIND, &
                            coef3sq3=5.19615242270663_MP_KIND                    ! acc_amps.F90:28
real(DS), parameter      :: coedpi6=0.523598775598299_DS, &
                            coedsq2p=2.506628274631_DS, &
                            coed3sq3=5.19615242270663_DS                         ! acc_amps.F90:32
```

`coedpi6 = π/6`, `coedsq2p = √(2π)`. `gg` = gravitational acceleration; `den_w` = density of water (**1000 kg/m³ in SI in `cal_Coalescence_Efficiency`**, but implicitly **1 g/cm³ in CGS** inside `cal_breakup_dis_LL`, where a drop diameter is `D=(mass/coedpi6)^(1/3)` — i.e. density folded to 1). `T_0` = freezing point; `imw`, `rmat`, `rmat_r` are mass-component indices (ice path only — not reached for token==1).

---

## 2. `cal_breakup_dis_LL` — THE core (Low & List 1982) — `mod_amps_core.F90:12019–12447` VERBATIM

```fortran
  subroutine cal_breakup_dis_LL(g_1,i,j,imin_bk,jmin_bk,bu_tmass,bu_fd,&
                                xD_L,xD_S,xS_T,xS_C,xCKE)
    use mod_amps_utility, only: &
       getznorm2
    ! ************************************************************************************
    ! This formulation is the original version by Low and List (1982), JAS.
    ! ************************************************************************************
!!c    use com_amps
    ! number of bins for the liquid spectrum
    type (Group), intent(in)               :: g_1
    ! catching drop and catched one.
    integer,intent(in)   :: i,j
    integer,intent(in) :: imin_bk,jmin_bk!,imax_bk,jmax_bk
    real(PS),intent(inout) :: bu_fd(2,*),bu_tmass(*)  ! 2014/10 T. Hashino modify for KID
    ! number of breakup drops
    !real(PS), intent(in)                    :: N_break
    ! diameter of large and small drops in m
    real(PS)                    :: xD_L,xD_S
    ! total energy of coalescence and collision kinetic energy (J)
    real(PS),intent(in)                    :: xCKE
    ! total surface energy (J)
    real(PS),intent(in)         :: xS_T,xS_C
    ! decrease in the surface energy (J)
    !real(PS)                    :: xdS_S
    ! diameter of large and small drops in m
    real(DS)                    :: D_L,D_S
    ! decrease in the surface energy (J)
    !real(DS)                    :: dS_S
    ! calculate the surface energy of the spherical equivalent of
    ! the united drop mass
    real(PS) :: S_C
    ! Weber number
    real(DS)                    :: W1,W2
    ! total energy of coalescence and collision kinetic energy (J)
    real(DS)                   :: CKE
    ! total surface energy (J)
    real(DS)         :: S_T
!!c    ! breakup distribution
!!c    real(PS) :: BD(81,2,81,81)

    real(DS),parameter          :: CKE0=8.93d-7
    real(DS),parameter          :: W0=0.86D0
    ! fraction of breakup
    real(DS)                    :: R_f,R_s,R_d


    ! fit parameters
    real(DS) :: mu_lnf,sig_lnf,mu_sf,H_sf,mu_lns,sig_lns,&
         mu_lnd,sig_lnd,sig_ld,mu_ls,sig_ls,mu_ld,&
         H_lf,mu_lf,sig_lf,lin_mu_lnf,lin_mu_lnd,lin_mu_lns

    real(DS) :: H_lnf,H_lns,H_lnd,sig_sf,H_ls,H_ld,P_MODE



    ! average number of fragments
    real(DS) :: F_f,F_s,F_d

    real(DS) :: D_1,D_2,D_log1,D_log2!,dir_del,D_S0
    real(DS),parameter :: app=1.02d+4, bpp=2.83D0
    ! low-diameter cut off related to the resolution of the experiments (cm)
    real(DS),parameter :: D_0=0.01D0
    ! diameter and mass of a coalescenced drop
    real(DS)  :: D_coal, m_coal
    real(DS)  :: mrat!,m_L

    ! total mass of the fragments of a large drop by sheet and disk breakup
    !real(DS) :: M_ls,M_ld

    real(DS),dimension(g_1%N_BIN) :: dmass,dcon !,dmass_s

    real(DS),dimension(3,g_1%N_BIN) :: m_f,m_s,m_d,n_f,n_s,n_d
    real(DS) :: frag_mass!,m_Lfrag
!!c    real(DS) :: getznorm2

    integer :: ibin,ibin_coal,kk,i1d_pair

    real(DS) :: dum,sq_twod
    real(DS) :: c1,c2,a1,x1,x2




    sq_twod=sqrt(2.0_DS)

    D_L=xD_L
    D_S=xD_S
    S_C=xS_C
    CKE=xCKE
    S_T=xS_T
    S_C=xS_C
    ! convert SI units to CGS except W, dS_S, CKE
    D_L=D_L*100.0_DS
    D_S=D_S*100.0_DS

    D_coal=(D_L**3.0+D_S**3.0)**(1.0/3.0)
    m_coal=coedpi6*D_coal**3.0

!!c    if(min(D_L,D_S)<=D_0) return
    if(D_coal<=D_0) return

    ! calculate variables necessary to calculate the break up probability
    ! Weber number
    W1=CKE/S_C
    W2=CKE/S_T


    ! 1. Determine the fraction of collision-breakup types
    ! fraction of the filament breakup
    if(CKE>=CKE0) then
       R_f=1.11e-4_DS*CKE**(-0.654)
    else
       R_f=1.0_DS
    end if
    ! fraction of sheet breakup
    if(W2>=W0) then
       R_s=0.685_DS*(1.0_DS-dexp(-1.63_DS*(W2-W0)))
    else
       R_s=0.0_DS
    end if
    if(R_s+R_f>1.0_DS) then
       R_s=1.0_DS-R_f
       R_d=0.0_DS
    else
       ! fraction of disk breakup
       R_d=dmax1(1.0_DS-R_f-R_s,0.0_DS)
    end if


    ! 2. calculate average number of fragments per a collision
    ! filament breakup
    F_f=(-2.25e+4_DS*(D_L-0.403_DS)*(D_L-0.403_DS)-37.9_DS)*D_S**2.5+&
         9.67_DS*(D_L+0.170_DS)*(D_L+0.170_DS)+4.95_DS
!tmp    F_f=(-2.25e+4_DS*(D_L-0.403_DS)**2.0-37.9_DS)*D_S**2.5+&
!tmp         9.67_DS*(D_L+0.170_DS)**2.0+4.95_DS
!tmp    write(*,*) "ck break",F_f,D_L,D_S,bpp,app
!tempei 2012/10    D_S0=((F_f-2.0)/bpp)**(1.0/app)
    F_f=dmax1(2.0_DS,dmin1(F_f,app*D_S**bpp+2.0_DS))

    ! sheet breakup
    F_s=dmax1(5.0_DS*&
         ! erf((S_T-2.53e-6_DS)/1.85e-6_DS)
         (2.0_DS*getznorm2(sq_twod*(S_T-2.53e-6_DS)/1.85e-6_DS)-1.0_DS)&
         +6.0_DS,2.0_DS)

    ! disk breakup
    F_d=dmax1(297.5_DS+23.76_DS*dlog(CKE),2.0_DS)


!!c    write(*,120) D_L,D_S,F_f,F_S,F_D,R_f,R_s,R_d
120 format("DL,DS,Ff,Fs,Fd,Rf,Rs,Rd",8ES12.3)



    ! 3. Calculate the parameters for parents distribution
    H_lf=50.8_DS*D_L**(-0.718)
    H_sf=4.18_DS*D_S**(-1.17)
    H_ls=100.0_DS*exp(-3.25_DS*D_S)
    H_ld=1.58e-5_DS*CKE**(-1.22)


    mu_lf=D_L
    mu_sf=D_S
    mu_ls=D_L
    mu_ld=D_L*(1.0_DS-exp(-3.70_DS*(3.10_DS-CKE/S_C)))

    call cal_sig_sf(D_0,R_f,1.0_DS,H_lf,mu_lf,sig_lf)
    call cal_sig_sf(D_0,R_f,1.0_DS,H_sf,mu_sf,sig_sf)
    call cal_sig_sf(D_0,R_s,1.0_DS,H_ls,mu_ls,sig_ls)
    call cal_sig_sf(D_0,R_d,1.0_DS,H_ld,mu_ld,sig_ld)

    ! 4. calculate the parameters for fragment distribution
    lin_mu_lnf=0.241_DS*D_S+0.0129_DS
    if(D_S<=D_0) then
       P_mode=1.68e+5_DS*D_S**2.33
    elseif(D_S>=1.2_DS*D_0) then
!tmp       P_mode=(43.4_DS*(D_L+1.81_DS)**2.0-159.0_DS)/D_S&
!tmp            -3870.0_DS*(D_L-0.285_DS)**2.0-58.1_DS
       P_mode=(43.4_DS*(D_L+1.81_DS)*(D_L+1.81_DS)-159.0_DS)/D_S&
            -3870.0_DS*(D_L-0.285_DS)*(D_L-0.285_DS)-58.1_DS
    else
       dum=(D_S-D_0)/(0.2_DS*D_0)
       P_mode=dum*(1.68e+5_DS*D_S**2.33)+(1.0_DS-dum)*&
            ((43.4_DS*(D_L+1.81_DS)*(D_L+1.81_DS)-159.0_DS)/D_S&
            -3870.0_DS*(D_L-0.285_DS)*(D_L-0.285_DS)-58.1_DS)
!tmp            ((43.4_DS*(D_L+1.81_DS)**2.0-159.0_DS)/D_S&
!tmp            -3870.0_DS*(D_L-0.285_DS)**2.0-58.1_DS)
    end if
    call cal_Hmusig(D_0,R_f,lin_mu_lnf,P_mode,F_f-2.0_DS,H_lnf,mu_lnf,sig_lnf)


    lin_mu_lns=0.254_DS*D_S**0.413*exp((3.53_DS*D_S-2.51_DS)*(D_L-D_S))
    P_mode=0.23_DS*D_S**(-3.93)*D_L**(14.2_DS*exp(-17.2_DS*D_S))
    call cal_Hmusig(D_0,R_s,lin_mu_lns,P_mode,F_s-1.0_DS,H_lns,mu_lns,sig_lns)

    lin_mu_lnd=exp((-17.4_DS*D_S-0.671_DS)*(D_L-D_S))*D_S
    if(D_L-D_S<0.5.and.0.007*D_S**(-2.54)>100.0) then
      P_mode=0.0_DS
    else
      P_mode=8.84_DS*D_S**(-2.52)*(D_L-D_S)**(0.007*D_S**(-2.54))
    endif

!!c    write(*,*) "check P_mode",P_mode,D_S,D_L
    call cal_Hmusig(D_0,R_d,lin_mu_lnd,P_mode,F_d-1.0_DS,H_lnd,mu_lnd,sig_lnd)


!!c    ! modification of mean to conserve total mass
!!c    call mod_mean_normal(mu_lf,sig_sf)

    ! calculate the fragment number concentration for each bin of liquid spectrum.

    Do ibin=max(i,j)+1,g_1%N_BIN+1
       if(g_1%binb(ibin)>m_coal) then
          ibin_coal=ibin-1
          exit
       end if
    end do

!!c    if(i==32.and.j==22) then
!!c	  write(*,*) "here"
!!c	end if
!!c    write(*,*) "i,j",i,j


    dcon=0.0_DS
    dmass=0.0_DS
!!c    dmass_s=0.0_DS

    m_f=0.0_DS
    m_s=0.0_DS
    m_d=0.0_DS
    n_f=0.0_DS
    n_s=0.0_DS
    n_d=0.0_DS

    frag_mass=0.0_DS
!!c    Do ibin=1,ibin_coal
    Do ibin=1,g_1%N_BIN
!!c       write(*,*) ibin,ibin_coal
       D_log2=dlog( (g_1%binb(ibin+1)/coedpi6)**(1.0/3.0))
       D_log1=dlog( (g_1%binb(ibin)/coedpi6)**(1.0/3.0))
       D_2=(g_1%binb(ibin+1)/coedpi6)**(1.0/3.0)
       D_1=(g_1%binb(ibin)/coedpi6)**(1.0/3.0)


       ! concentration produced by breakup
       ! ------------ filament breakup -------------
       ! lognormal distribution of the fragment drops
       n_f(1,ibin)=dmax1(0.0_DS,H_lnf*sig_lnf*coedsq2p*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lnf)/sig_lnf)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lnf)/sig_lnf)))

       ! normal distribution of the small parent drop
       n_f(2,ibin)=H_sf*sig_sf*coedsq2p*(&
            getznorm2((D_2-mu_sf)/sig_sf)-&
            getznorm2((D_1-mu_sf)/sig_sf))

       ! normal distribution of the large parent drop
       n_f(3,ibin)=H_lf*sig_lf*coedsq2p*(&
            getznorm2((D_2-mu_lf)/sig_lf)-&
            getznorm2((D_1-mu_lf)/sig_lf))


       ! ------------ sheet breakup -----------------
       ! lognormal distribution of the fragment drops
       n_s(1,ibin)=dmax1(0.0_DS,H_lns*sig_lns*coedsq2p*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lns)/sig_lns)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lns)/sig_lns)))

       ! no distribution for the small parent drop
       n_s(2,ibin)=0.0_DS


       ! normal distribution of the large parent drop
       n_s(3,ibin)=H_ls*sig_ls*coedsq2p*(&
            getznorm2((D_2-mu_ls)/sig_ls)-&
            getznorm2((D_1-mu_ls)/sig_ls))
       ! ------------ disc breakup -----------------
       ! lognormal distribution of the fragment drops
       n_d(1,ibin)=dmax1(0.0_DS,H_lnd*sig_lnd*coedsq2p*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lnd)/sig_lnd)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lnd)/sig_lnd)))

       ! no distribution of the small parent drop
       n_d(2,ibin)=0.0_DS

       ! normal distribution of the large parent drop
       n_d(3,ibin)=H_ld*sig_ld*coedsq2p*(&
            getznorm2((D_2-mu_ld)/sig_ld)-&
            getznorm2((D_1-mu_ld)/sig_ld))




       ! ------------ filament breakup -------------
       ! lognormal distribution of the fragment drops
!       m_f(1,ibin)=dmax1(0.0_DS,coedpi6*H_lnf*sig_lnf*coedsq2p*dexp(4.5_DS*sig_lnf**2.0+3.0_DS*mu_lnf)*(&
       m_f(1,ibin)=dmax1(0.0_DS,coedpi6*H_lnf*sig_lnf*coedsq2p*dexp(4.5_DS*sig_lnf*sig_lnf+3.0_DS*mu_lnf)*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lnf)/sig_lnf-3.0_DS*sig_lnf)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lnf)/sig_lnf-3.0_DS*sig_lnf)))

       ! normal distribution of the small parent drop
       x1=(D_1-mu_sf)/sig_sf
       x2=(D_2-mu_sf)/sig_sf
       c1=sig_sf*(sig_sf*sig_sf*(x1*x1+2.0_DS)+3.0_DS*mu_sf*(sig_sf*x1+mu_sf))/coedsq2p
       c2=sig_sf*(sig_sf*sig_sf*(x2*x2+2.0_DS)+3.0_DS*mu_sf*(sig_sf*x2+mu_sf))/coedsq2p
!       c1=sig_sf*(sig_sf**2.0*(x1**2.0+2.0_DS)+3.0_DS*mu_sf*(sig_sf*x1+mu_sf))/coedsq2p
!       c2=sig_sf*(sig_sf**2.0*(x2**2.0+2.0_DS)+3.0_DS*mu_sf*(sig_sf*x2+mu_sf))/coedsq2p

!       a1=mu_sf*(3.0_DS*sig_sf**2.0+mu_sf**2.0)
       a1=mu_sf*(3.0_DS*sig_sf*sig_sf+mu_sf*mu_sf)
!       m_f(2,ibin)=coedpi6*H_sf*sig_sf*coedsq2p*(-c2*dexp(-x2**2.0/2.0_DS)+c1*dexp(-x1**2.0/2.0_DS)+&
!            a1*(getznorm2(x2)-getznorm2(x1)))
       m_f(2,ibin)=coedpi6*H_sf*sig_sf*coedsq2p*(-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))


       ! normal distribution of the large parent drop
       x1=(D_1-mu_lf)/sig_lf
       x2=(D_2-mu_lf)/sig_lf
!       c1=sig_lf*(sig_lf**2.0*(x1**2.0+2.0_DS)+3.0_DS*mu_lf*(sig_lf*x1+mu_lf))/coedsq2p
!       c2=sig_lf*(sig_lf**2.0*(x2**2.0+2.0_DS)+3.0_DS*mu_lf*(sig_lf*x2+mu_lf))/coedsq2p
!       a1=mu_lf*(3.0_DS*sig_lf**2.0+mu_lf**2.0)
!       m_f(3,ibin)=coedpi6*H_lf*sig_lf*coedsq2p*(-c2*dexp(-x2**2.0/2.0_DS)+c1*dexp(-x1**2.0/2.0_DS)+&
!            a1*(getznorm2(x2)-getznorm2(x1)))

       c1=sig_lf*(sig_lf*sig_lf*(x1*x1+2.0_DS)+3.0_DS*mu_lf*(sig_lf*x1+mu_lf))/coedsq2p
       c2=sig_lf*(sig_lf*sig_lf*(x2*x2+2.0_DS)+3.0_DS*mu_lf*(sig_lf*x2+mu_lf))/coedsq2p
       a1=mu_lf*(3.0_DS*sig_lf*sig_lf+mu_lf*mu_lf)
       m_f(3,ibin)=coedpi6*H_lf*sig_lf*coedsq2p*(-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))


       ! ------------ sheet breakup -----------------
       ! lognormal distribution of the fragment drops
!       m_s(1,ibin)=dmax1(0.0_DS,coedpi6*H_lns*sig_lns*coedsq2p*dexp(4.5_DS*sig_lns**2.0+3.0_DS*mu_lns)*(&
       m_s(1,ibin)=dmax1(0.0_DS,coedpi6*H_lns*sig_lns*coedsq2p*dexp(4.5_DS*sig_lns*sig_lns+3.0_DS*mu_lns)*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lns)/sig_lns-3.0_DS*sig_lns)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lns)/sig_lns-3.0_DS*sig_lns)))
       ! for the small parent drop
       m_s(2,ibin)=0.0_DS

       ! normal distribution of the large parent drop
       x1=(D_1-mu_ls)/sig_ls
       x2=(D_2-mu_ls)/sig_ls
!       c1=sig_ls*(sig_ls**2.0*(x1**2.0+2.0_DS)+3.0_DS*mu_ls*(sig_ls*x1+mu_ls))/coedsq2p
!       c2=sig_ls*(sig_ls**2.0*(x2**2.0+2.0_DS)+3.0_DS*mu_ls*(sig_ls*x2+mu_ls))/coedsq2p
!       a1=mu_ls*(3.0_DS*sig_ls**2.0+mu_ls**2.0)
!       m_s(3,ibin)=coedpi6*H_ls*sig_ls*coedsq2p*(-c2*dexp(-x2**2.0/2.0_DS)+c1*dexp(-x1**2.0/2.0_DS)+&
       c1=sig_ls*(sig_ls*sig_ls*(x1*x1+2.0_DS)+3.0_DS*mu_ls*(sig_ls*x1+mu_ls))/coedsq2p
       c2=sig_ls*(sig_ls*sig_ls*(x2*x2+2.0_DS)+3.0_DS*mu_ls*(sig_ls*x2+mu_ls))/coedsq2p
       a1=mu_ls*(3.0_DS*sig_ls*sig_ls+mu_ls*mu_ls)
       m_s(3,ibin)=coedpi6*H_ls*sig_ls*coedsq2p*&
          (-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))

       ! ------------ disc breakup -----------------
       ! lognormal distribution of the fragment drops
!       m_d(1,ibin)=dmax1(0.0_DS,coedpi6*H_lnd*sig_lnd*coedsq2p*dexp(4.5_DS*sig_lnd**2.0+3.0_DS*mu_lnd)*(&
       m_d(1,ibin)=dmax1(0.0_DS,coedpi6*H_lnd*sig_lnd*&
          coedsq2p*dexp(4.5_DS*sig_lnd*sig_lnd+3.0_DS*mu_lnd)*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lnd)/sig_lnd-3.0_DS*sig_lnd)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lnd)/sig_lnd-3.0_DS*sig_lnd)))

       ! for the small drop
       m_d(2,ibin)=0.0

       ! normal distribution of the large parent drop
       x1=(D_1-mu_ld)/sig_ld
       x2=(D_2-mu_ld)/sig_ld
!       c1=sig_ld*(sig_ld**2.0*(x1**2.0+2.0_DS)+3.0_DS*mu_ld*(sig_ld*x1+mu_ld))/coedsq2p
!       c2=sig_ld*(sig_ld**2.0*(x2**2.0+2.0_DS)+3.0_DS*mu_ld*(sig_ld*x2+mu_ld))/coedsq2p
!       a1=mu_ld*(3.0_DS*sig_ld**2.0+mu_ld**2.0)
!       m_d(3,ibin)=coedpi6*H_ld*sig_ld*coedsq2p*(-c2*dexp(-x2**2.0/2.0_DS)+c1*dexp(-x1**2.0/2.0_DS)+&
!            a1*(getznorm2(x2)-getznorm2(x1)))

       c1=sig_ld*(sig_ld*sig_ld*(x1*x1+2.0_DS)+3.0_DS*mu_ld*(sig_ld*x1+mu_ld))/coedsq2p
       c2=sig_ld*(sig_ld*sig_ld*(x2*x2+2.0_DS)+3.0_DS*mu_ld*(sig_ld*x2+mu_ld))/coedsq2p
       a1=mu_ld*(3.0_DS*sig_ld*sig_ld+mu_ld*mu_ld)
       m_d(3,ibin)=coedpi6*H_ld*sig_ld*coedsq2p*&
            (-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))

    end do

    Do ibin=1,g_1%N_BIN
       dmass(ibin)=R_f*(m_f(1,ibin)+m_f(2,ibin)+m_f(3,ibin))+&
            R_s*(m_s(1,ibin)+m_s(2,ibin)+m_s(3,ibin))+&
            R_d*(m_d(1,ibin)+m_d(2,ibin)+m_d(3,ibin))
       dcon(ibin)=R_f*(n_f(1,ibin)+n_f(2,ibin)+n_f(3,ibin))+&
            R_s*(n_s(1,ibin)+n_s(2,ibin)+n_s(3,ibin))+&
            R_d*(n_d(1,ibin)+n_d(2,ibin)+n_d(3,ibin))

!!c       dmass_s(ibin)=R_f*(m_f(1,ibin)*0.5_PS+m_f(2,ibin))+&
!!c            R_s*(m_s(1,ibin)*0.5_PS+m_s(2,ibin))+&
!!c            R_d*(m_d(1,ibin)+m_d(2,ibin)+m_d(3,ibin))


       if(dcon(ibin)<1.0e-100_DS.or.dmass(ibin)<1.0e-100_DS) then
          dcon(ibin)=0.0_DS
          dmass(ibin)=0.0_DS
       else
          if(dmass(ibin)/dcon(ibin)>g_1%binb(ibin+1).or.&
               dmass(ibin)/dcon(ibin)<g_1%binb(ibin))then
!!c             write(*,*) "here the bin shift needed",ibin,D_L*100.0_RP,D_S*100.0_RP
!!c             return
!!c             dcon(ibin)=2.0_DS*dmass(ibin)/(g_1%binb(ibin+1)+g_1%binb(ibin))
             dmass(ibin)=(g_1%binb(ibin+1)+g_1%binb(ibin))*0.5_DS*dcon(ibin)
          end if
       end if
    end do


    i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2
    bu_tmass(i1d_pair)=m_coal

    mrat=m_coal/sum(dmass)
    Do ibin=1,g_1%N_BIN
       kk=(i1d_pair-1)*g_1%N_BIN+ibin
       bu_fd(2,kk)=mrat*dcon(ibin)
       bu_fd(1,kk)=mrat*dmass(ibin)
!!c       bu_fd(3,kk)=
!!c       BD(ibin,2,i,j)=mrat*dcon(ibin)
!!c       BD(ibin,1,i,j)=mrat*dmass(ibin)
    end do

!tmp    write(*,*) "i1d_pair,kk",i1d_pair,kk

  end subroutine cal_breakup_dis_LL
```

**What it reads from the `Group`:** only `g_1%binb(1:N_BIN+1)` (bin mass boundaries, CGS grams) and `g_1%N_BIN`. Output: `bu_tmass(i1d_pair)` (= coalesced mass) and `bu_fd(1:2, kk)` (fragment mass/con per destination bin). The `i1d_pair` / `kk` index arithmetic must be reproduced exactly (matches the fact-file `i1d_pair_max`/`kk_max` sizing in the driver).

---

## 3. `cal_Coalescence_Efficiency` — token==1 (liquid-liquid) branch — `mod_amps_core.F90:11796–12017` VERBATIM

The full routine is quoted; **only the `g_1%token==1 .and. g_2%token==1` branch (lines 11842–11881) is on the breakup path** — the `token==2`/ice branches (11882–12016) are dead for cal_breakfragment (liquid group, token=1) and can be dropped from the port.

```fortran
  subroutine cal_Coalescence_Efficiency(g_1,i,g_2,j,n,th_var,E_coal,&
       D_L,D_S,S_T,S_C,DS_S,CKE)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the coalescence efficiency.
    !
    ! Assumptions
    ! For liquid hydros: formula by Low and List (1982)
    ! For liquid and solid: 1.0
    ! For solids: 1.0 ( need to incorporate "sticking mechanism")
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Group), intent(in)               :: g_1, g_2
    integer, intent(in)                    :: i, j,n
    type (Thermo_Var) :: th_var
    ! total coalescence efficiency
    real(PS)                    :: E_coal
    ! interlocking coalscence efficiency and sticking coalescence efficiency
    real(PS)                    :: E_int,E_stick
    ! diameter and velocity of large and small drops in SI units
    real(PS)                    :: D_L, D_S, V_L, V_S
    ! total energy of coalescence and collision kinetic energy in SI units
    real(PS)                    :: E_T, CKE
    real(PS)                    :: S_T, S_C, dS_S
    ! density of water (kg/m^3)
    real(PS), parameter         :: den_w = 1000.0D0
    real(PS), parameter         :: a = 0.778D0
    real(PS), parameter         :: b = 2.61d+6
    ! bulk density of hexagonal ice crystal calculated with circumscribing sphere
    !real(PS)        :: den_hex1,den_hex2
    ! low-diameter cut off related to the resolution of the experiments (cm)
    real(PS),parameter :: D_0=0.01D0
    ! mass limits for interlocking and stickiness efficiencies.
    ! This practically shut down aggregation process for massive particles.
!!c    real(PS),parameter :: mlimit_int=1.0,mlimit_stick=1.0
!!c    real(PS),parameter :: mlimit_int=1.0e-2,mlimit_stick=1.0
!!c    real(PS),parameter :: mlimit_int=1.0e-2,mlimit_stick=1.0e-2
    real(PS),parameter :: mlimit_int=1.0d-3,mlimit_stick=1.0d-2

    real (PS), parameter          :: stick_lmt = 0.5_PS

    ! equivalent lengths
    !real (PS) :: d1,d2

    ! reference density for polycrystals assuming phi=0.25
    real(PS),parameter :: den_hex_pol=0.259572D0


    if( g_1%token == 1 .and. g_2%token == 1 ) then
       ! define large diameter and small diameter (m)
       D_L = max( g_1%MS(i,n)%len, g_2%MS(j,n)%len ) * 1.0e-2
       D_S = min( g_1%MS(i,n)%len, g_2%MS(j,n)%len ) * 1.0e-2

       if(min(D_L,D_S)*100.0<D_0) then
         E_coal=1.0_PS

       else

         ! define large and small terminal velocity (m/s)
         V_L = max( g_1%MS(i,n)%vtm, g_2%MS(j,n)%vtm ) * 1.0e-2
         V_S = min( g_1%MS(i,n)%vtm, g_2%MS(j,n)%vtm ) * 1.0e-2

         ! calculate the total surface energy of incident drops
         S_T = PI*th_var%sig_wa*1.0e-3_DS*( D_L**2.0 + D_S**2.0)

         ! calculate the surface energy of the spherical equivalent of
         ! the united drop mass
         S_C = PI*th_var%sig_wa*1.0e-3_DS*( D_L**3.0 + D_S**3.0)**(2.0/3.0)

         ! decrease in surface energy
         dS_S = S_T - S_C

         ! calculate collision kinetic energy
         CKE = (den_w*PI/12.0_PS)*((V_L-V_S)**2.0)*&
              ((D_L*D_S)**3.0)/(D_L**3.0+D_S**3.0)

         ! calculate the total energy of coalescence
         E_T = CKE + dS_S

         if( E_T < 5.0e-6_PS ) then
            E_coal = a * ((1.0_PS + D_S/D_L)**(-2.0)) * &
               exp( - b*th_var%sig_wa*1.0e-3_DS*(E_T**2.0)/S_c)
         else
            E_coal = 0.0_PS
         end if

       end if

    else if(( g_1%token == 1 .and.  g_2%token == 2 ) .or. &
    ...  [ice/mixed branches 11882–12016 — NOT reached for cal_breakfragment] ...
    end if
  end subroutine cal_Coalescence_Efficiency
```

**token==1 inputs:** `%len` (bin diameter, CGS cm), `%vtm` (drop terminal velocity, CGS cm/s), and `th_var%sig_wa` (surface tension of water, dyn/cm — one scalar, function of T only). **Outputs used by the driver:** `CKE`, `S_T`, `S_C`, `D_L`, `D_S` (in **SI** — note `*1.0e-2` cm→m and `sig_wa*1.0e-3` conversions; `cal_breakup_dis_LL` re-converts `D_L,D_S` back to CGS internally via `*100`). `E_coal`/`DS_S` are computed but the driver only gates on `CKE<=1.0e-20`.

---

## 4. Sigma solvers — `mod_amps_core.F90:27254–27354` VERBATIM

```fortran
  subroutine cal_sig_sf(D_0,R,Nx,H_s,mu_s,sig_s)
    use mod_amps_utility, only: &
       zbrent
    real(DS), intent(in) :: D_0,R,H_s,mu_s,Nx!,D_coal
    real(DS), intent(inout) :: sig_s
!!c    real(DS) :: zbrent,bsig
    !real(DS) :: bsig
    integer,parameter :: ITER=10
    integer :: i

    sig_s=1.0_DS/(H_s*coedsq2p)

    if(R<1.0e-20_DS) return


!!c    do i=1,ITER
!!c       bsig=sig_s
!!c       sig_s=1.0_DS/(H_s*coedsq2p)/&
!!c           (getznorm2((D_coal-mu_s)/sig_s )-getznorm2((D_0-mu_s)/sig_s ))
!!c       if(abs((bsig-sig_s)/sig_s)<0.01_DS) then
!!c          exit
!!c       end if
!!c    end do
    do i=1,ITER
!!c       sig_s=zbrent(FSIG,sig_s*1.0e-8,sig_s*10.0_RP,1.0e-4_DS)
       sig_s=zbrent(1,0.0_DS,Nx,H_s,D_0,mu_s,sig_s,0.0_DS,0.0_DS,&
            sig_s*1.0e-8_DS,sig_s*10.0_DS,1.0e-4_DS)

       if(sig_s==-999.9_DS) then
          sig_s=1.0_DS/(H_s*coedsq2p)/10.0**real(i,PS_KIND)
       else
          exit
       end if
    end do

    [commented-out FSIG contains-function omitted — see source :27290-27301]
  end subroutine cal_sig_sf

  subroutine cal_Hmusig(D_0,R,lin_mu,P_mode,Nx,H,mu,sig)
    use mod_amps_utility, only: &
       zbrent
    real(DS),intent(in) :: D_0,R,lin_mu,P_mode,Nx
    real(DS),intent(inout) :: H,mu,sig
    integer,parameter :: ITER=10
    integer :: i
!!c    real(DS) :: zbrent,osig
    real(DS) :: osig

    sig=10.0_DS*lin_mu
    if(Nx<=1.0e-20.or.P_mode<=1.0e-20.or.R<=1.0e-20) then
       H=0.0_DS
       mu=0.0_DS
       return
    end if

    do i=1,ITER
       osig=sig
!!c       sig=zbrent(FSIG2,sig*1.0e-5,sig*10.0_RP,1.0e-4_DS)
       sig=zbrent(2,P_mode,Nx,H,D_0,lin_mu,sig,0.0_DS,0.0_DS,&
            sig*1.0e-5,sig*10.0_RP,1.0e-4_DS)
       if(sig==-999.9_DS) then
          sig=osig*1.5_RP
       else
          exit
       end if
    end do

    H=P_mode*lin_mu*exp(0.5_DS*sig**2)
    mu=log(lin_mu)+sig**2

    [commented-out FSIG2 contains-function omitted — see source :27339-27351]
  end subroutine cal_Hmusig
```

Note: the `-999.9` sentinel is the `zbrent` non-bracketing fallback; both routines retry up to `ITER=10` times, shrinking (`cal_sig_sf`) or growing (`cal_Hmusig` ×1.5) the sigma seed. `cal_sig_sf` passes `iwhich=1` to `zbrent`; `cal_Hmusig` passes `iwhich=2`. (`iwhich=3` is used only by the ice-mass path, not by cal_breakfragment.)

---

## 5. `zbrent` (Brent root-finder + internal `func`) — `mod_amps_utility.F90:12778–12907` VERBATIM

```fortran
  function zbrent(iwhich,p_mode,nx,h,d_0,mu,sig,m_l,d_coal,x1,x2,tol) result(out)
    implicit none
    integer :: iwhich
    integer :: itmax
    real(ds) :: out,tol,x1,x2,eps,p_mode,nx,h,d_0,mu,sig,m_l,d_coal

    parameter (itmax=100,eps=3.e-8)
    !  using brent's method, find the root of a function func known to lie between x1 and x2.
    !  the root, returned as zbrent, will be refined until its accuracy is tol.

    !  parameters: maximum allowed number of iterations, and machine floating-point precision.
    integer :: iter
    real(ds) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    if((fa.gt.0.0_ds.and.fb.gt.0.0_ds).or.(fa.lt.0.0_ds.and.fb.lt.0.0_ds)) then
!!c     write(*,*) "root must be bracketed for zbrent"
       out=-999.9_ds
       return
    end if

    c=b
    fc=fb
    do iter=1,itmax
       if((fb.gt.0.0_ds.and.fc.gt.0.0_ds).or.(fb.lt.0.0_ds.and.fc.lt.0.0_ds))then
          c=a ! rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2.0_ds*eps*abs(b)+0.5_ds*tol ! convergence check.
       xm=0.5_ds*(c-b)
       if(abs(xm).le.tol1 .or. fb.eq.0.)then
          out=b
          return
       endif
       if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa ! attempt inverse quadratic interpolation.
          if(a.eq.c) then
             p=2.0_ds*xm*s
             q=1.0_ds-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.0_ds))
             q=(q-1.0_ds)*(r-1.0_ds)*(s-1.0_ds)
          endif
          if(p.gt.0.) q=-q ! check whether in bounds.
          p=abs(p)
          if(2.0_ds*p .lt. min(3.0_ds*xm*q-abs(tol1*q),abs(e*q))) then
             e=d ! accept interpolation.
             d=p/q
          else
             d=xm ! interpolation failed, use bisection.
             e=d
          endif
       else ! bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b ! move last best guess to a.
       fa=fb
       if(abs(d) .gt. tol1) then ! evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       fb=func(b)
    enddo
    if(debug) write(*,*) "brent exceeding maximum iterations"

    out=b
    return
  contains
    function func(x) result(out)
      real(ds) :: x,out
      real(ds) :: c1,c2,a1,mass,z_coal,z_0
      if(iwhich==1) then
!!c      out=1.0_ds-(getznorm2((d_coal-mu_s)/x )-getznorm2((d_0-mu_s)/x))
!!c         out=1.0_ds-getznorm2((d_coal-mu_s)/x )

         out=x-nx/h/coedsq2p/(&
              1.0_ds-&
              min(0.99999_DS,getznorm2((d_0-mu)/max(x,1.0e-20_DS))))
!org              min(0.99999999_DS,getznorm2((d_0-mu)/max(x,1.0e-20_DS))))
!tempei bug              getznorm2((d_0-mu)/x))
      elseif(iwhich==2) then
         h=p_mode*mu*exp(0.5_ds*x**2)
         c1=log(mu)+x**2
         out=x-nx/h/coedsq2p/(&
              1.0_ds-&
!tempei bug              getznorm2((dlog(d_0)-c1)/x))
              min(0.99999999_DS,getznorm2((dlog(d_0)-c1)/max(x,1.0e-20_DS))))
      elseif(iwhich==3) then
         z_coal=(d_coal-x)/sig
         z_0=-x/sig

         c1=sig*(sig**2*(z_0**2+2.0_ds)+3.0_ds*x*(sig*z_0+x))/coedsq2p
         c2=sig*(sig**2*(z_coal**2+2.0_ds)+3.0_ds*x*(sig*z_coal+x))/coedsq2p
         a1=x*(3.0_ds*sig**2.0+x**2)

         if(abs(z_coal)>1.0e+30_DS) then
            c2=0.0_ds
         else
            c2=c2*dexp(-z_coal**2/2.0_ds)
         end if

         if(abs(z_0)>1.0e+30_DS) then
            c1=0.0_ds
         else
            c1=c1*dexp(-z_0**2/2.0_ds)
         end if

         mass=-c2+c1+a1*(getznorm2(z_coal)-getznorm2(z_0))
         out=(pi/6.0_ds)*mass-m_l*(getznorm2(z_coal)-getznorm2(z_0))

      end if
    end function func

  end function zbrent
```

**Port note:** `func` mutates `h` in-place for `iwhich==2` (`h=p_mode*mu*exp(0.5*x²)`). In `cal_Hmusig`, after `zbrent` returns, `H` is recomputed once more from the converged `sig` (`H=P_mode*lin_mu*exp(0.5*sig²)`), so the in-loop mutation is transient — the port can pass `h` by value into the objective and recompute `H` at the end. `iwhich==3` is dead for cal_breakfragment.

---

## 6. Normal-CDF special functions — `mod_amps_utility.F90` VERBATIM

### 6.1 `getznorm2` (:369–375) and hot path

```fortran
      function getznorm2(x) result(p)
        real(DS) :: x,mean,sd,q,bound,p!,out
        integer :: status
        mean=0.0_DS
        sd=1.0_DS
        call cdfnor (1, p, q, x, mean, sd, status, bound )
      end function getznorm2
```

`getznorm2(x)` = standard-normal CDF Φ(x). Always calls `cdfnor(which=1)` → `cumnor`. **This is the entire breakup hot path for the normal CDF.** A direct port can replace the whole chain with `0.5*erfc(-x/√2)` if bit-agreement with `cumnor`'s Cody rational approximation is not required; if bit-agreement IS required (recommended for LUT reproducibility), port `cumnor` verbatim (§6.3).

### 6.2 `cdfnor` (:10033–10230)

Full routine at `mod_amps_utility.F90:10033`. For `which==1` (the only value `getznorm2` ever passes) the body reduces to:

```fortran
  if ( which == 1 ) then
    z = ( x - mean ) / sd
    call cumnor ( z, p, q )
  ...
```

All the argument-range checks (lines 10126–10197) and the `which==2/3/4` inverse branches (which call `dinvnr`) are **not exercised** by `getznorm2`. Port `cumnor` and a thin `which==1` wrapper; the range checks and inverse branches can be omitted for the breakup port.

### 6.3 `cumnor` (:10231–10435) — Cody rational-approximation CDF — VERBATIM (the routine to port for bit-agreement)

```fortran
subroutine cumnor ( arg, cum, ccum )
  implicit none
  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287D+00, 1.6102823106855587881D+02, &
    1.0676894854603709582D+03, 1.8154981253343561249D+04, &
    6.5682337918207449113D-02 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870D+01, 9.7609855173777669322D+02, &
    1.0260932208618978205D+04, 4.5507789335026729956D+04 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764D-01, 8.8831497943883759412D+00, &
    9.3506656132177855979D+01, 5.9727027639480026226D+02, &
    2.4945375852903726711D+03, 6.8481904505362823326D+03, &
    1.1602651437647350124D+04, 9.8427148383839780218D+03, &
    1.0765576773720192317D-08 /)
  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691D+01, 2.3538790178262499861D+02, &
    1.5193775994075548050D+03, 6.4855582982667607550D+03, &
    1.8615571640885098091D+04, 3.4900952721145977266D+04, &
    3.8912003286093271411D+04, 1.9685429676859990727D+04 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  integer i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699D-01, 1.274011611602473639D-01, &
    2.2235277870649807D-02, 1.421619193227893466D-03, &
    2.9112874951168792D-05, 2.307344176494017303D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121D+00, 4.68238212480865118D-01, &
    6.59881378689285515D-02, 3.78239633202758244D-03, &
    7.29751555083966205D-05 /)
  real ( kind = 8 ), parameter :: root32 = 5.656854248D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: sqrpi = 3.9894228040143267794D-01
  real ( kind = 8 ), parameter :: thrsh = 0.66291D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) xsq

  eps = epsilon ( 1.0D+00 ) * 0.5D+00

  x = arg
  y = abs ( x )

  if ( y <= thrsh ) then
    if ( eps < y ) then
      xsq = x * x
    else
      xsq = 0.0D+00
    end if
    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do
    cum = x * ( xnum + a(4) ) / ( xden + b(4) )
    temp = cum
    cum = 0.5D+00 + temp
    ccum = 0.5D+00 - temp
  else if ( y <= root32 ) then
    xnum = c(9) * y
    xden = y
    do i = 1, 7
      xnum = ( xnum + c(i) ) * y
      xden = ( xden + d(i) ) * y
    end do
    cum = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum
    if ( 0.0D+00 < x ) then
      call d_swap ( cum, ccum )
    end if
  else
    cum = 0.0D+00
    xsq = 1.0D+00 / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    do i = 1, 4
      xnum = ( xnum + p(i) ) * xsq
      xden = ( xden + q(i) ) * xsq
    end do
    cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    cum = ( sqrpi - cum ) / y
    xsq = aint ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( - del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum
    if ( 0.0D+00 < x ) then
      call d_swap ( cum, ccum )
    end if
  end if

  if ( cum < tiny ( cum ) ) then
    cum = 0.0D+00
  end if
  if ( ccum < tiny ( ccum ) ) then
    ccum = 0.0D+00
  end if

  return
end subroutine cumnor
```

`d_swap` (:10436–10466) is a trivial two-value swap — in Python just `cum, ccum = ccum, cum`.

### 6.4 Inverse-CDF chain — NOT on the breakup hot path (quoted for completeness)

`dinvnr` (:10467–10548), `stvaln` (:11698–11764), `eval_pol` (:11765–11808). These are reached only via `cdfnor(which=2/3/4)`, which `getznorm2` never calls. The breakup port does **not** need them. Verbatim source:

`stvaln` uses two coefficient arrays and calls `eval_pol` (Horner scheme):
```fortran
  xden = (/ 0.993484626060D-01, 0.588581570495D+00, 0.531103462366D+00, &
            0.103537752850D+00, 0.38560700634D-02 /)   ! dimension(0:4)
  xnum = (/ -0.322232431088D+00, -1.000000000000D+00, -0.342242088547D+00, &
            -0.204231210245D-01, -0.453642210148D-04 /) ! dimension(0:4)
  ! p<=0.5: sgn=-1, z=p ; else sgn=+1, z=1-p
  y = sqrt(-2*log(z));  stvaln = sgn*( y + eval_pol(xnum,4,y)/eval_pol(xden,4,y) )
```
`eval_pol(a,n,x)` = Horner: `term=a(n); do i=n-1,0,-1: term=term*x+a(i)`.
`dinvnr(p,q)`: `pp=min(p,q); strtx=stvaln(pp)`; Newton iteration (maxit=100, eps=1e-13) `dx=(cum-pp)/(r2pi*exp(-0.5*xcur²))`, `r2pi=0.3989422804014326`, returning `±xcur` by sign of `p<=q`. Full source at `mod_amps_utility.F90:10467–10548` / `:11698–11808`.

---

## 7. Data flow into the LL table (from cal_breakfragment driver)

The driver (fact-file §A.1) does, for each bin pair `(i,j)`, `imin_bk<=i<=imax_bk=NRBIN`, `jmin_bk<=j<=i-1`:

1. `cal_Coalescence_Efficiency(liquid,i,liquid,j,1, steady%TV(1), E_coal, D_L,D_S,S_T,S_C,DS_S,CKE)` → needs `liquid%MS(i,1)%len,%vtm`, `liquid%MS(j,1)%len,%vtm`, and `sig_wa`.
2. `if (CKE<=1.0e-20) cycle`.
3. `cal_breakup_dis_LL(liquid,i,j,imin_bk,jmin_bk, bu_tmass,bu_fd, D_L,D_S,S_T,S_C,CKE)` → needs `liquid%binb(:)`, `liquid%N_BIN`, and the `(D_L,D_S,S_T,S_C,CKE)` from step 1.

Index bookkeeping (from driver): `jmin_bk` = first bin with `%len>=D_0` (`D_0=0.01 cm`); `imin_bk=jmin_bk+1`; `imax_bk=NRBIN`; `jmax_bk=NRBIN-1`. `i1d_pair_max=(imax_bk-1)-jmin_bk+1 + (imax_bk-imin_bk)*(1+imax_bk-imin_bk)/2`; `kk_max=i1d_pair_max*N_BIN`. Fixed air state: `T=278.6795 K`, `PT=850 hPa`, `W=0`, `phase=1`, `RH=100%`.

---

## 8. MINIMAL `%len` / `%vtm` — replace `diag_pq` (Cluster 2)

`cal_breakfragment` currently calls `diag_pq` (mod_amps_core.F90:12552, ~648 lines + the whole property-vec stack) **only** to fill `liquid%MS(:,1)%len` and `%vtm`. For a liquid (`token==1`) group these reduce to two closed-form scalars per bin. **M2b should compute these two directly and skip `diag_pq` entirely.**

### 8.1 `%len` — geometric equivalent diameter (CGS, cm)

For a pure liquid drop the general `diag_pq` length (`cal_den_aclen_vec`, `class_Group.F90:11253`)

```fortran
g%MS(i,n)%len = (6.0_PS*(g%MS(i,n)%mean_mass/(PI*g%MS(i,n)%den)))**(1.0/3.0)
```

collapses, with drop density `den = den_w = 1 g/cm³` (CGS) and `mean_mass` = the bin's representative mass, to

```
len_i = ( mean_mass_i / (π/6) )^(1/3)        [cm]   ≡ (mean_mass_i / coedpi6)^(1/3)
```

This is **exactly the same diameter formula `cal_breakup_dis_LL` uses internally** on the bin boundaries (`D=(binb/coedpi6)^(1/3)`, source lines 12258–12261), so it is self-consistent by construction. The `max(r_n*1.05, len)` aerosol-core floor (`:11259`) is negligible for the breakable rain bins (`len>=D_0=0.01 cm = 100 µm`) and can be dropped for the table. `mean_mass_i` is the bin representative mass (bin geometric-center mass); for the fixed cloudlab grid this is deterministic from `binb`/`srat_r`/`minmass_r`.

### 8.2 `%vtm` — liquid-drop terminal velocity — USE THE ALREADY-PORTED FUNCTION

The `%vtm` that `cal_Coalescence_Efficiency` reads is the **per-bin (mean-mass-point) Beard-regime fall speed**, `cal_terminal_vel_vec` **phase==1**, `class_Group.F90:8069–8119`:

```fortran
    if( phase == 1 ) then
       do n = 1, g%L
       do i = 1, g%N_BIN
        if(icond1(i,n)==0) then
          rad=g%MS(i,n)%len*0.5_PS
          if(rad<0.5e-4_PS) then                          ! bug found 2016/05/30
            g%MS(i,n)%vtm=0.0_PS
          elseif(0.5e-4_PS<=rad.and.rad<10.0e-4_PS) then  ! Stokes regime
            U_S=rad**2.0*gg*(den_w-ag%TV(n)%den_a)/4.5_PS/ag%TV(n)%d_vis
            lambda_a=6.6e-6_PS*(ag%TV(n)%d_vis/1.818e-4_PS)*(1013250.0_PS/ag%TV(n)%P)*(ag%TV(n)%T/293.15_PS)
            g%MS(i,n)%vtm=(1.0_PS+1.26_PS*lambda_a/rad)*U_S
          elseif(rad<535.0e-4_PS) then                    ! Beard intermediate (Davies number)
            X=log(32.0_PS*rad**3.0*(den_w-ag%TV(n)%den_a)*ag%TV(n)%den_a*gg/3.0_PS/ag%TV(n)%d_vis**2.0)
            Y=-0.318657e+1_PS+0.992696_PS*X-0.153193e-2_PS*X*X-0.987059e-3_PS*X*X*X &
               -0.578878e-3_PS*X*X*X*X+0.855176e-4*X*X*X*X*X-0.327815e-5*X*X*X*X*X*X
            g%MS(i,n)%Nre=exp(Y)
            g%MS(i,n)%vtm=g%MS(i,n)%Nre*ag%TV(n)%d_vis/(2.0_PS*rad*ag%TV(n)%den_a)
          else                                            ! Beard large-drop (Bond/physical-property number)
            rad=min(rad,3500.0e-4_RP)
            NBO=gg*(den_w-ag%TV(n)%den_a)*rad**2.0/ag%TV(n)%sig_wa
            NP=ag%TV(n)%sig_wa**3.0*ag%TV(n)%den_a**2.0/ag%TV(n)%d_vis**4.0/gg
            X=log(NBO*NP**(1.0/6.0)*16.0_PS/3.0_PS)
            Y=-0.500015e+1_PS+0.523778e+1_PS*X-0.204914e+1_PS*X*X+0.475294_PS*X*X*X &
               -0.542819e-1_PS*X*X*X*X+0.238449e-2_PS*X*X*X*X*X
            g%MS(i,n)%Nre=NP**(1.0/6.0)*exp(Y)
            g%MS(i,n)%vtm=g%MS(i,n)%Nre*ag%TV(n)%d_vis/(2.0_PS*rad*ag%TV(n)%den_a)
          end if
        endif
      enddo
      enddo
```

**CONFIRMED: `core/liquid_diag.py::_terminal_velocity` IS the correct `%vtm`.** Per `docs/superpowers/facts/m1/fortran-mapping.md:63`, `core/liquid_diag.py` (`diag_pq_liquid`/`LiquidDiag`) ports exactly `cal_terminal_vel_vec` **phase==1** at `class_Group.F90:8069–8119` — the same lines quoted above — and the fact explicitly distinguishes it from the bin-integrated Böhm `cal_wterm_vel_v3_vec` sedimentation fall speed (`w_terminal_vel`, M2c scope), which is a *different* function. The breakup `%vtm` is the per-bin diag one → the ported `_terminal_velocity` is the right call.

**Caveat for M2b — the air-state arguments.** `_terminal_velocity` needs `den_a` (dry-air density), `d_vis` (dynamic viscosity), `P`, `T`, `sig_wa`, all from `ag%TV(n)` = the single `steady` thermo state built in the driver at the fixed `T=278.6795 K, PT=850 hPa, RH=100%` (via `make_AirGroup_2`/`make_thermo_var3_2`). M2b must supply these four/five scalars evaluated at that fixed state (they are constants for the table — one `ThermoState`), NOT advected fields. `sig_wa` at 278.68 K is also the surface tension consumed by `cal_Coalescence_Efficiency`. So the entire Cluster-2 dependency of the table collapses to: **one fixed `ThermoState` (T, P, den_a, d_vis, sig_wa) + the closed-form `len_i` + the already-ported `_terminal_velocity(len_i, thermo)`**.

### 8.3 Net port shape for M2b

```
for each rain bin i:  len_i  = (mean_mass_i / coedpi6)**(1/3)                 # §8.1, CGS cm
                      vtm_i  = liquid_diag._terminal_velocity(len_i, thermo0) # §8.2, ALREADY PORTED
for each pair (i>j>=jmin_bk):
    D_L,D_S,S_T,S_C,CKE = coalescence_eff_token1(len_i,vtm_i,len_j,vtm_j, sig_wa0)   # §3, token==1 only
    if CKE <= 1e-20: continue
    breakup_dis_LL(binb, N_BIN, i,j, imin_bk,jmin_bk, D_L,D_S,S_T,S_C,CKE) -> bu_tmass, bu_fd  # §2
        # internally: cal_sig_sf/cal_Hmusig (§4) -> zbrent (§5) -> getznorm2->cumnor (§6)
```

where `thermo0` / `sig_wa0` are the single fixed-state thermo scalars (T=278.6795 K, P=850 hPa, RH=100%). **No `diag_pq`, no ventilation/capacitance/vapor-dep/surface-temp/ice-shape.**

---

## 9. Numerical-fidelity risks M2b must budget for

- `cumnor` is a **piecewise Cody rational approximation** with three branches (`|x|<=0.66291`, `<=√32`, `>√32`) and an `aint(x*16)/16` range-reduction trick. Bit-agreement with `scipy.special.ndtr`/`erfc` is **not** guaranteed — port `cumnor` verbatim if the dumped-vs-computed LUT must match the reference AMPS run to tight tolerance.
- `zbrent` returns the `-999.9` sentinel on non-bracketing; `cal_sig_sf`/`cal_Hmusig` retry with rescaled seeds. The **exact retry schedule** (`sig/10^i` shrink vs `×1.5` grow, `ITER=10`) must be reproduced or the converged sigmas — hence the whole fragment distribution — diverge.
- `cal_breakup_dis_LL` sums lognormal/normal moment integrals with `dexp(-x²/2)` terms and `dmax1(0, …)` clamps; the `dmass(ibin)/dcon(ibin)` bin-boundary reprojection (`:12421-12427`) and the final `mrat=m_coal/sum(dmass)` mass-normalization (`:12435`) are order-sensitive.
- This is exactly the catastrophic-cancellation / summation-order divergence class the fact file (§A.4) cites as the reason it recommends **dumping** `bu_fd`/`bu_tmass` rather than porting — the user has overridden that toward a port, so these three items are the fidelity budget. If M2b ports, the recommended validation is a bit-diff of `bu_fd`/`bu_tmass` against a dumped reference AMPS run at the same `NRBIN`/bin grid.