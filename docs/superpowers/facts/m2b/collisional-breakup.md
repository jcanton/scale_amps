All sections extracted. Here is the consolidated document.

---

# AMPS Warm-Rain Collisional Breakup (Low–List 1982) — Verbatim Extraction

File (all snippets): `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_core.F90`
Table build / config: `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_lib.F90`, `com_amps.F90`

## Architecture summary (how the flag reaches the code)

`micexfg(18)` is the KID collisional-breakup switch (`com_amps.F90:100-101`, `class_Cloud_Micro.F90:1017` passes `CM%micexfg(18)`). When it is `1`, the fragment lookup tables `bu_fd` / `bu_tmass` are **precomputed once** by `cal_breakup_dis_LL` (called only from `mod_amps_lib.F90:1992`, never at runtime), and at runtime the collision routine sets the internal `ibreak` flag which gates every branch below. `cal_breakup_dis_LL` itself is a table builder; `add_fragments_col_vec`, and the `i1d_pair`/`kk` index math, are the runtime consumers.

---

## 1. How breakup couples into coalescence (the `ibreak` path)

### 1a. Subroutine signature and table declarations (`mod_amps_core.F90:1061-1100`)

```fortran
  subroutine coalescence(g_1, g_2, ag, level,col_level,mes_rc,ibreak &
              ,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd,ID,JD,KD &
```
```fortran
    integer,intent(in)  :: ibreak
```
```fortran
    real(PS),intent(in) :: bu_fd(2,*),bu_tmass(*)  ! 2014/10 T. Hashino modify for KID
```

### 1b. Kernel call — `ibreak` forwarded into the kernel (`mod_amps_core.F90:1600-1604`)

```fortran
       do n=1,g_1%L
          if(icycle_n(n)==0) then
            call cal_collision_kernel_func(g_1,g_2,ag%TV(n),col_level,ibreak,&
               imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,&
               n,KC(:,:,n),E_coal(:,:,n), &
```

### 1c. Consumed number/mass by breakup — the fraction `(1 - E_coal)` of collisions becomes fragments (`mod_amps_core.F90:1668-1707`)

```fortran
       ! calculate consumed con and mass by breakup process
       if(ibreak==1) then
         do i=1,g_1%N_BIN
            do n = 1, g_2%L
            do j = 1, g_2%N_BIN
             if(icycle_ijn(i,j,n)==0) then
               N_bk=N_col(i,j,n)*max(0.0_PS,1.0_PS-E_coal(i,j,n))
               used_N_b(j,n)=used_N_b(j,n)+N_bk
             endif
           enddo
           enddo
         enddo
         do j=1,g_2%N_BIN
            do n = 1, g_1%L
            do i = 1, g_1%N_BIN
             if(icycle_ijn(i,j,n)==0) then
               N_bk=N_col(i,j,n)*max(0.0_PS,1.0_PS-E_coal(i,j,n))
               used_N_b(i,n)=used_N_b(i,n)+N_bk
             endif
           enddo
           enddo
         enddo

         do n = 1, g_1%L
         do i = 1, g_1%N_BIN
           if(icycle_n(n)==0) then
             used_N_2(i,n)=used_N_2(i,n)+used_N_b(i,n)
             used_M_2(i,n)=used_M_2(i,n)+used_N_b(i,n)*g_2%MS(i,n)%mean_mass
           endif
         enddo
         enddo
       end if
```

The same `N_bk=N_col*max(0.0_PS,1.0_PS-E_coal)` accounting is repeated inside the over-depletion fixer under `if(ibreak==1)` at lines 1754-1798 (concentration fix) and 1839-1883 (mass fix).

### 1d. The dispatch to the fragment adder (`mod_amps_core.F90:2884-2893`)

```fortran
         ! ********************************************************************
         ! calculation of concentration and mass transfer
         ! by collisional breakup
         ! ********************************************************************
         if(pro_type==2.and.ibreak==1) then

           call add_fragments_col_vec(new_N_1,new_M_1, &
                    i,g_1,g_2,icond1,used_marker,E_coal,N_col, &
                    imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd)
         endif
```

So: a rain-rain collision produces a coalesced drop with probability `E_coal`; the complementary fraction `(1-E_coal)` is removed as `used_N_b`/`used_M_2` (breakup consumption) and re-injected as a fragment spectrum via `add_fragments_col_vec` using the precomputed `bu_fd`/`bu_tmass` tables.

---

## 2. `cal_breakup_dis_LL` — FULL (Low–List 1982 filament/sheet/disk) (`mod_amps_core.F90:12019-12447`)

This routine is > 300 lines; quoted below are the full argument list + declarations, the fraction/number-of-fragments logic, the three fragment-mode branches with their exact distribution formulas and coefficients, and the bin mapping.

### 2a. Signature + declarations (`12019-12097`)

```fortran
  subroutine cal_breakup_dis_LL(g_1,i,j,imin_bk,jmin_bk,bu_tmass,bu_fd,&
                                xD_L,xD_S,xS_T,xS_C,xCKE)
    use mod_amps_utility, only: &
       getznorm2
    ! ************************************************************************************
    ! This formulation is the original version by Low and List (1982), JAS.
    ! ************************************************************************************
    ! number of bins for the liquid spectrum
    type (Group), intent(in)               :: g_1
    ! catching drop and catched one.
    integer,intent(in)   :: i,j
    integer,intent(in) :: imin_bk,jmin_bk!,imax_bk,jmax_bk
    real(PS),intent(inout) :: bu_fd(2,*),bu_tmass(*)  ! 2014/10 T. Hashino modify for KID
    ! diameter of large and small drops in m
    real(PS)                    :: xD_L,xD_S
    ! total energy of coalescence and collision kinetic energy (J)
    real(PS),intent(in)                    :: xCKE
    ! total surface energy (J)
    real(PS),intent(in)         :: xS_T,xS_C
    ! diameter of large and small drops in m
    real(DS)                    :: D_L,D_S
    real(PS) :: S_C
    ! Weber number
    real(DS)                    :: W1,W2
    real(DS)                   :: CKE
    real(DS)         :: S_T

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

    real(DS) :: D_1,D_2,D_log1,D_log2
    real(DS),parameter :: app=1.02d+4, bpp=2.83D0
    ! low-diameter cut off related to the resolution of the experiments (cm)
    real(DS),parameter :: D_0=0.01D0
    ! diameter and mass of a coalescenced drop
    real(DS)  :: D_coal, m_coal
    real(DS)  :: mrat

    real(DS),dimension(g_1%N_BIN) :: dmass,dcon

    real(DS),dimension(3,g_1%N_BIN) :: m_f,m_s,m_d,n_f,n_s,n_d
    real(DS) :: frag_mass

    integer :: ibin,ibin_coal,kk,i1d_pair

    real(DS) :: dum,sq_twod
    real(DS) :: c1,c2,a1,x1,x2
```

### 2b. Unit conversion, coalesced drop, Weber numbers (`12102-12123`)

```fortran
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

    if(D_coal<=D_0) return

    ! Weber number
    W1=CKE/S_C
    W2=CKE/S_T
```

### 2c. Mode fractions R_f / R_s / R_d (`12126-12145`)

```fortran
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
```

### 2d. Average number of fragments F_f / F_s / F_d (`12148-12165`)

```fortran
    ! filament breakup
    F_f=(-2.25e+4_DS*(D_L-0.403_DS)*(D_L-0.403_DS)-37.9_DS)*D_S**2.5+&
         9.67_DS*(D_L+0.170_DS)*(D_L+0.170_DS)+4.95_DS
    F_f=dmax1(2.0_DS,dmin1(F_f,app*D_S**bpp+2.0_DS))

    ! sheet breakup
    F_s=dmax1(5.0_DS*&
         ! erf((S_T-2.53e-6_DS)/1.85e-6_DS)
         (2.0_DS*getznorm2(sq_twod*(S_T-2.53e-6_DS)/1.85e-6_DS)-1.0_DS)&
         +6.0_DS,2.0_DS)

    ! disk breakup
    F_d=dmax1(297.5_DS+23.76_DS*dlog(CKE),2.0_DS)
```

### 2e. Parent-distribution parameters (`12173-12188`)

```fortran
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
```

### 2f. Fragment-distribution parameters — filament / sheet / disk lognormal modes (`12190-12222`)

```fortran
    ! 4. calculate the parameters for fragment distribution
    lin_mu_lnf=0.241_DS*D_S+0.0129_DS
    if(D_S<=D_0) then
       P_mode=1.68e+5_DS*D_S**2.33
    elseif(D_S>=1.2_DS*D_0) then
       P_mode=(43.4_DS*(D_L+1.81_DS)*(D_L+1.81_DS)-159.0_DS)/D_S&
            -3870.0_DS*(D_L-0.285_DS)*(D_L-0.285_DS)-58.1_DS
    else
       dum=(D_S-D_0)/(0.2_DS*D_0)
       P_mode=dum*(1.68e+5_DS*D_S**2.33)+(1.0_DS-dum)*&
            ((43.4_DS*(D_L+1.81_DS)*(D_L+1.81_DS)-159.0_DS)/D_S&
            -3870.0_DS*(D_L-0.285_DS)*(D_L-0.285_DS)-58.1_DS)
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

    call cal_Hmusig(D_0,R_d,lin_mu_lnd,P_mode,F_d-1.0_DS,H_lnd,mu_lnd,sig_lnd)
```

### 2g. Locate the coalesced-drop bin (`12230-12235`)

```fortran
    Do ibin=max(i,j)+1,g_1%N_BIN+1
       if(g_1%binb(ibin)>m_coal) then
          ibin_coal=ibin-1
          exit
       end if
    end do
```

### 2h. Per-bin NUMBER distributions — the three modes, each with fragment (lognormal), small parent, large parent (`12256-12308`)

```fortran
    Do ibin=1,g_1%N_BIN
       D_log2=dlog( (g_1%binb(ibin+1)/coedpi6)**(1.0/3.0))
       D_log1=dlog( (g_1%binb(ibin)/coedpi6)**(1.0/3.0))
       D_2=(g_1%binb(ibin+1)/coedpi6)**(1.0/3.0)
       D_1=(g_1%binb(ibin)/coedpi6)**(1.0/3.0)

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
```

### 2i. Per-bin MASS distributions — filament mode (fragment lognormal 3rd-moment + small/large parent normals) (`12313-12349`)

```fortran
       ! ------------ filament breakup -------------
       ! lognormal distribution of the fragment drops
       m_f(1,ibin)=dmax1(0.0_DS,coedpi6*H_lnf*sig_lnf*coedsq2p*dexp(4.5_DS*sig_lnf*sig_lnf+3.0_DS*mu_lnf)*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lnf)/sig_lnf-3.0_DS*sig_lnf)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lnf)/sig_lnf-3.0_DS*sig_lnf)))

       ! normal distribution of the small parent drop
       x1=(D_1-mu_sf)/sig_sf
       x2=(D_2-mu_sf)/sig_sf
       c1=sig_sf*(sig_sf*sig_sf*(x1*x1+2.0_DS)+3.0_DS*mu_sf*(sig_sf*x1+mu_sf))/coedsq2p
       c2=sig_sf*(sig_sf*sig_sf*(x2*x2+2.0_DS)+3.0_DS*mu_sf*(sig_sf*x2+mu_sf))/coedsq2p

       a1=mu_sf*(3.0_DS*sig_sf*sig_sf+mu_sf*mu_sf)
       m_f(2,ibin)=coedpi6*H_sf*sig_sf*coedsq2p*(-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))

       ! normal distribution of the large parent drop
       x1=(D_1-mu_lf)/sig_lf
       x2=(D_2-mu_lf)/sig_lf
       c1=sig_lf*(sig_lf*sig_lf*(x1*x1+2.0_DS)+3.0_DS*mu_lf*(sig_lf*x1+mu_lf))/coedsq2p
       c2=sig_lf*(sig_lf*sig_lf*(x2*x2+2.0_DS)+3.0_DS*mu_lf*(sig_lf*x2+mu_lf))/coedsq2p
       a1=mu_lf*(3.0_DS*sig_lf*sig_lf+mu_lf*mu_lf)
       m_f(3,ibin)=coedpi6*H_lf*sig_lf*coedsq2p*(-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))
```

### 2j. Per-bin MASS distributions — sheet mode (`12352-12373`)

```fortran
       ! ------------ sheet breakup -----------------
       ! lognormal distribution of the fragment drops
       m_s(1,ibin)=dmax1(0.0_DS,coedpi6*H_lns*sig_lns*coedsq2p*dexp(4.5_DS*sig_lns*sig_lns+3.0_DS*mu_lns)*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lns)/sig_lns-3.0_DS*sig_lns)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lns)/sig_lns-3.0_DS*sig_lns)))
       ! for the small parent drop
       m_s(2,ibin)=0.0_DS

       ! normal distribution of the large parent drop
       x1=(D_1-mu_ls)/sig_ls
       x2=(D_2-mu_ls)/sig_ls
       c1=sig_ls*(sig_ls*sig_ls*(x1*x1+2.0_DS)+3.0_DS*mu_ls*(sig_ls*x1+mu_ls))/coedsq2p
       c2=sig_ls*(sig_ls*sig_ls*(x2*x2+2.0_DS)+3.0_DS*mu_ls*(sig_ls*x2+mu_ls))/coedsq2p
       a1=mu_ls*(3.0_DS*sig_ls*sig_ls+mu_ls*mu_ls)
       m_s(3,ibin)=coedpi6*H_ls*sig_ls*coedsq2p*&
          (-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))
```

### 2k. Per-bin MASS distributions — disc mode (`12375-12402`)

```fortran
       ! ------------ disc breakup -----------------
       ! lognormal distribution of the fragment drops
       m_d(1,ibin)=dmax1(0.0_DS,coedpi6*H_lnd*sig_lnd*&
          coedsq2p*dexp(4.5_DS*sig_lnd*sig_lnd+3.0_DS*mu_lnd)*(&
            getznorm2((dmin1(D_log2,dlog(D_coal))-mu_lnd)/sig_lnd-3.0_DS*sig_lnd)-&
            getznorm2((dmax1(D_log1,dlog(D_0))-mu_lnd)/sig_lnd-3.0_DS*sig_lnd)))

       ! for the small drop
       m_d(2,ibin)=0.0

       ! normal distribution of the large parent drop
       x1=(D_1-mu_ld)/sig_ld
       x2=(D_2-mu_ld)/sig_ld
       c1=sig_ld*(sig_ld*sig_ld*(x1*x1+2.0_DS)+3.0_DS*mu_ld*(sig_ld*x1+mu_ld))/coedsq2p
       c2=sig_ld*(sig_ld*sig_ld*(x2*x2+2.0_DS)+3.0_DS*mu_ld*(sig_ld*x2+mu_ld))/coedsq2p
       a1=mu_ld*(3.0_DS*sig_ld*sig_ld+mu_ld*mu_ld)
       m_d(3,ibin)=coedpi6*H_ld*sig_ld*coedsq2p*&
            (-c2*dexp(-x2*x2/2.0_DS)+c1*dexp(-x1*x1/2.0_DS)+&
            a1*(getznorm2(x2)-getznorm2(x1)))

    end do
```

### 2l. Combine modes weighted by R_f/R_s/R_d, mean-mass consistency clamp (`12404-12429`)

```fortran
    Do ibin=1,g_1%N_BIN
       dmass(ibin)=R_f*(m_f(1,ibin)+m_f(2,ibin)+m_f(3,ibin))+&
            R_s*(m_s(1,ibin)+m_s(2,ibin)+m_s(3,ibin))+&
            R_d*(m_d(1,ibin)+m_d(2,ibin)+m_d(3,ibin))
       dcon(ibin)=R_f*(n_f(1,ibin)+n_f(2,ibin)+n_f(3,ibin))+&
            R_s*(n_s(1,ibin)+n_s(2,ibin)+n_s(3,ibin))+&
            R_d*(n_d(1,ibin)+n_d(2,ibin)+n_d(3,ibin))

       if(dcon(ibin)<1.0e-100_DS.or.dmass(ibin)<1.0e-100_DS) then
          dcon(ibin)=0.0_DS
          dmass(ibin)=0.0_DS
       else
          if(dmass(ibin)/dcon(ibin)>g_1%binb(ibin+1).or.&
               dmass(ibin)/dcon(ibin)<g_1%binb(ibin))then
             dmass(ibin)=(g_1%binb(ibin+1)+g_1%binb(ibin))*0.5_DS*dcon(ibin)
          end if
       end if
    end do
```

### 2m. Fragments → table (the write-side index math + mass renormalization `mrat`) (`12432-12447`)

```fortran
    i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2
    bu_tmass(i1d_pair)=m_coal

    mrat=m_coal/sum(dmass)
    Do ibin=1,g_1%N_BIN
       kk=(i1d_pair-1)*g_1%N_BIN+ibin
       bu_fd(2,kk)=mrat*dcon(ibin)
       bu_fd(1,kk)=mrat*dmass(ibin)
    end do

  end subroutine cal_breakup_dis_LL
```

So `bu_fd(1,kk)` = mass in bin, `bu_fd(2,kk)` = number in bin, both scaled so total fragment mass = coalesced mass `m_coal` stored in `bu_tmass(i1d_pair)`.

---

## 3. Runtime consumers — FULL

### 3a. `add_fragments_col_vec` (`mod_amps_core.F90:15659-15774`)

```fortran
  subroutine add_fragments_col_vec(new_N_1,new_M_1, &
                    i,g_1,g_2,icond1,used_marker,E_coal,N_col, &
                    imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd)

    type (Group), intent(in)   :: g_1, g_2
    real(8), dimension(mxnbin,*),intent(inout)                :: new_N_1
    real(8), dimension(mxnbin,g_1%L,1+mxnmasscomp),intent(inout)  :: new_M_1
    ! catching bin
    integer,intent(in) :: i
    integer,dimension(mxnbin,*),intent(in) :: icond1
    integer,dimension(mxnbin,*),intent(in)      :: used_marker
    ! coalescence efficiency
    real(PS), dimension(mxnbin,mxnbin,*),intent(in)          :: E_coal
    ! total number of collected hydrometeors in j bin by i bin, N_col(i,j)
    real(PS),dimension(mxnbin,mxnbin,*),intent(in)   ::  N_col
    integer,intent(in) :: imin_bk,imax_bk,jmin_bk,jmax_bk
    real(PS),intent(in) :: bu_fd(2,*),bu_tmass(*)  ! 2014/10 T. Hashino modify for KID

    ! number of collitional breakup
    real(PS)         :: N_bk
    real(PS)      :: ratio_M

    integer            :: j,k,kk,n,l,i1d_pair
    real(PS) :: mod_rat
    real(PS) :: mass_add

    if(i<imin_bk.or.i>imax_bk) return

    do j=1,i-1
      if(j<jmin_bk.or.j>jmax_bk) cycle

      i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2
      if(bu_tmass(i1d_pair)<=1.0e-30_PS) cycle

      do n = 1, g_1%L
      do k = 1, g_1%n_bin

         if(used_marker(i,n)==1.or.icond1(i,n)==1) then
            kk=(i1d_pair-1)*g_1%N_BIN+k

            N_bk=max(0.0_PS,1.0_PS-E_coal(i,j,n))*N_col(i,j,n)
            mod_rat=(g_1%MS(i,n)%mean_mass+g_2%MS(j,n)%mean_mass)/bu_tmass(i1d_pair)

            new_N_1(k,n)=new_N_1(k,n)+bu_fd(2,kk)*N_bk*mod_rat
            mass_add=bu_fd(1,kk)*N_bk*mod_rat
            new_M_1(k,n,1)=new_M_1(k,n,1)+mass_add

            ! assume that the breakup fragments have the average property
            ! this is not really physical.
            ! In reality, the fragments distribution is made up of mass from small and large drop.
            ! The ratio of mass from two drops depends on CKE and other physical parameters.
            ! Doing this take more memory in lookup tables.

            do l=1,g_1%N_masscom
               ratio_M=0.0_PS
               if(l==rmat_m) then
                  ratio_M=(g_1%MS(i,n)%mass(rmat)+g_1%MS(j,n)%mass(rmat))/&
                          (g_1%MS(i,n)%mass(rmt)+g_1%MS(j,n)%mass(rmt))
               elseif(l==rmas_m) then
                  ratio_M=(g_1%MS(i,n)%mass(rmas)+g_1%MS(j,n)%mass(rmas))/&
                          (g_1%MS(i,n)%mass(rmt)+g_1%MS(j,n)%mass(rmt))
               elseif(l==rmai_m) then
                  ratio_M=max(0.0_PS,&
                          (g_1%MS(i,n)%mass(rmat)+g_1%MS(j,n)%mass(rmat))/&
                          (g_1%MS(i,n)%mass(rmt)+g_1%MS(j,n)%mass(rmt)) -&
                          (g_1%MS(i,n)%mass(rmas)+g_1%MS(j,n)%mass(rmas))/&
                          (g_1%MS(i,n)%mass(rmt)+g_1%MS(j,n)%mass(rmt)) &
                          )
               endif
               new_M_1(k,n,1+l)=new_M_1(k,n,1+l)+ratio_M*mass_add
            end do
         endif
      end do
      end do

    end do

  end subroutine add_fragments_col_vec
```

### 3b. `P_breakup` (`mod_amps_core.F90:19873-19904`)

```fortran
  function P_breakup(phase, a_star, max_dim) result(out)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate probability for the parent drop to break up in time dt
    ! based on Kombayashi et al. (1964)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of water
    ! 1: water, 2: ice
    integer, intent(in)    :: phase
    ! radius of parent hydrometeor in mm
    real(PS), intent(in)  :: a_star
    ! maximum dimension, (radius), of rain or aggregate (cm)
    real(PS), intent(in)             :: max_dim
    real(PS)             :: k, out
    if( phase == 1 ) then
       ! case of liquid
       ! the probability becomes about 1 at 4.5 mm
       if( a_star >= max_dim ) then
          out = 1.0_PS
       else
          out = min( 2.94e-7_PS*exp(3.4_PS*a_star*10.0_PS), 1.0_PS )
       end if
    else if( phase == 2 ) then
       ! case of ice
       ! NOTE: assume that one coefficient is the same as liquid case
       if( a_star >= max_dim ) then
          out = 1.0_PS
       else
          k = 15.03968607_PS/(max_dim*10.0_PS)
          out = min(2.94e-7_PS*exp(k*a_star*10.0_PS),1.0_PS)
       end if
    end if
  end function P_BREAKUP
```

### 3c. `Q_breakup2` (`mod_amps_core.F90:19953-19997`)

```fortran
  function Q_breakup2(phase, a_star, a, m, AA, BB,switch) result(out)
    use mod_amps_utility, only: &
       fGM
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate (switch = 1): number of drops of mass between m and m+dm
    !           (switch = 2): total number of drops between a and m (#/cm^3)
    !                         a < m cm
    !           (switch = 3): total mass of drops between a and m (a < m cm)
    ! based on Srivastava (1971)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! phase of water
    ! 1: water, 2: ice
    integer, intent(in)    :: phase
    real(PS),intent(in)  :: a_star, a, m,AA,BB
    integer, intent(in)  :: switch
    real(PS)             :: out, x1, x2

    if( phase == 1 ) then
       if( switch == 1 ) then
          out = (AA*BB/3.0_PS/m)*(a/a_star)*exp(-BB*a/a_star)
       else if( switch == 2 ) then
          out = -AA*( exp(-BB*m/a_star) - exp(-BB*a/a_star))
       else if( switch == 3 ) then
          !     a_2
          ! M = \  (4*PI/3)*den_w*a^3 Q(a*,a) da
          !     a_1
          !
          ! +++ change radius and to parameter +++
          !     x = B a/a_star
          x1 = BB*a/a_star
          x2 = BB*m/a_star

          out = (4.0_PS*PI/3.0_PS)*AA*((a_star/BB)**3.0)*&
               fGM(4.0_PS, x1, x2)

       end if
    else if( phase == 2 ) then
       if( switch == 1 ) then
          out = (AA*BB/3.0_PS/m)*(a/a_star)*exp(-BB*a/a_star)
       else if( switch == 2 ) then
          out = -AA*( exp(-BB*m/a_star) - exp(-BB*a/a_star))
       end if
    end if

  end function Q_BREAKUP2
```

Note: `Q_breakup2` (Srivastava 1971) and `P_breakup` (Komabayasi 1964) are used by the **spontaneous / large-drop breakup** paths (`mod_amps_core.F90:10057-10182`, `10424-10495`, `11526-11528`), a separate mechanism from the Low–List **collisional** breakup tables; both are included per the request.

---

## 4. The fragment-table indexing at runtime (the `i1d_pair` / `kk` math)

The single formula tying a `(bin_i, bin_j)` collision pair to the precomputed table (a packed lower-triangular index over `i>j`, then flattened over the fragment bin `k`) appears identically in all four places — build side and both runtime consumers:

Build side (`cal_breakup_dis_LL`, `mod_amps_core.F90:12432,12437`):
```fortran
    i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2
    ...
       kk=(i1d_pair-1)*g_1%N_BIN+ibin
```

Runtime — fragment injection (`add_fragments_col_vec`, `mod_amps_core.F90:15699,15709`):
```fortran
      i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2
      if(bu_tmass(i1d_pair)<=1.0e-30_PS) cycle
      ...
            kk=(i1d_pair-1)*g_1%N_BIN+k
```

Runtime — kernel efficiency gate (`cal_collision_kernel_func`, `mod_amps_core.F90:16341-16359`):
```fortran
    if(ibreak==1) then
       do j = 1, g_2%N_BIN
       do i = 1, g_1%N_BIN
        if(icond1(i,j)==0) then
          if(i<imin_bk.or.i>imax_bk.or.j<jmin_bk.or.j>jmax_bk) then
          else
            i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2

            if(bu_tmass(i1d_pair)<=1.0e-30_PS) then
              E_coal(i,j)=1.0_PS
            end if
          end if
        end if
      end do
      end do
    endif
```
(If no fragment mass is tabulated for a pair, `E_coal` is forced to 1 — i.e. pure coalescence, no breakup.)

### Index bounds / table sizing (`mod_amps_lib.F90:1949-1962` and `438-449`)

```fortran
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
```
```fortran
      if (nbr == 40) then
         allocate(bu_fd(2,17835))
         allocate(bu_tmass(435))
      else if (nbr == 80) then
         allocate(bu_fd(2,62400))
         allocate(bu_tmass(780))
      endif
      bu_fd = 0.0_PS
```

`bu_tmass` is indexed by `i1d_pair` (one entry per collidable `(i>j)` pair, triangular count = 435 for 40 rain bins / 780 for 80); `bu_fd(1:2, kk)` is indexed by `kk = (i1d_pair-1)*N_BIN + bin`, storing (mass, number) of fragments deposited into each fragment bin.