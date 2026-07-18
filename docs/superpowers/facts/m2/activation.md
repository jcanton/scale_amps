# AMPS CCN Activation + Vapor Advancement — Production Path (`act_type=1` → `cal_aptact_var8_kc04dep`)

All code verbatim from `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_core.F90` (subroutine spans lines **5776–9388**). Helper functions from `class_Mass_Bin.F90`.

---

## 0. Signature, key declarations, module constants

**Argument list** (`mod_amps_core.F90:5776–5835`):
```fortran
  subroutine cal_aptact_var8_kc04dep(level,ag,ncat_a,ga,gr,gs &
                          ,flagp_a,flagp_r,flagp_s,iflg_inuc,iflg_dep,iflg_dhf & ! iflg: 10, 13 ,19
                          ,mes_rc &
                          ,qtp,thil,nu_aps,phi_aps,M_aps &
                          ,ap_sig_cp,ap_mean_cp,CCNMAX,CRIC_RN_IMM,frac_dust& !,cdf_cp_180m0
                          ,snrml,ID,JD,KD)
```
Groups: `ag` (air), `ga(ncat_a)` (aerosol categories — only cat 1 for CCN, cat 2 for deposition IN), `gr` (rain/liquid), `gs` (ice). Flags: `iflg_inuc` (=10), `iflg_dep` (=13), `iflg_dhf` (micexfg =19). Category 2 is the deposition-nucleation dust category.

**Solver control constants** (`5978–5979`, and the inner Brent's own set at `9020–9024`):
```fortran
    INTEGER,parameter :: ITMAX=200
    real(PS),PARAMETER :: alim1=1.0e-5,alim2=1.0e-5
```
```fortran
    real(PS),parameter :: max_sw=0.1,min_ac=1.0e-7,min_sw=0.0001   ! 6036
    real(PS),parameter :: max_r_n=10.0                              ! 6039
    real(PS),parameter :: max_r_n_DHF=10.0e-4                       ! 6043
    real(PS),parameter :: INMAX=1.7e-3_PS                           ! 6054
    real(PS),parameter :: max_fact=1.0_PS                           ! 6057
```
Key limiter/flag arrays (`5990`, `6019–6022`): `noindep,noccnt,nodhft`; `akk_lmt`, `akk_lmt_DHF`, `akk`, `akk_DHF`. Bin partition `zn` / `fact_c` (`6307–6311`). Coefficients (`6314–6321`):
```fortran
    coef_a=2.0_PS/(den_w*R_v)
    coef_al1=g*L_e*Rdvchiarui*M_a/(c_pa*R_u)
    coef_al2=-g*M_a/R_u
    rm_frg1=coef4pi3*r_frg1**3
```

---

## 1. `cal_aptact_var8_kc04dep` — OUTER STRUCTURE

### 1a. Grid-box skip mask (`6323–6408`)
Loop `do n=1,ag%L` sets `icycle_n(n)`; a box is processed only if water is present and either liquid (`flagp_r>0 & s_v(1)>0`) or ice-nucleation conditions hold:
```fortran
      if( mes_rc(n) == 0 ) then
        icycle_n(n)=1
      elseif(mes_rc(n)==1) then
        if( qtp(n)>1.0e-20_PS.and.&
              ((flagp_r>0.and.ag%TV(n)%s_v(1)>0.0_PS).or.&
               (flagp_s>0.and.iflg_inuc>0.and.(&
                   (iflg_dep>0.and.ag%TV(n)%s_v(2)>0.0_PS.and.ag%TV(n)%T<T_0).or.&
                   (iflg_dhf>0.and.ag%TV(n)%s_v(2)>0.0_PS.and.ag%TV(n)%T<T_0)  ) &
                  )&
               ) ) then
        else
          icycle_n(n)=1
        endif
      endif
```
`if(all(icycle_n(1:ag%L)==1)) ... return` (early exit, `6344–6408`).

### 1b. Per-box init (`6410–6547`)
Zeroes all `used_*`, `gain_*`, sets `akk_lmt(n)=1`, `akk_lmt_DHF(n)=1`, `sw(n)=ag%TV(n)%s_v(1)`, `sw_n=sw_m=sw_b=sw`, `T_a_n=T_a_b=ag%TV(n)%T`, `noccnt=nodhft=noindep=1`, `imethod=0`, `iphase=0`. Accumulates `qr_0,qi_0,nr_0,ni_0` and freezing/riming gains over rain bins (`6472–6489`) and ice bins (`6491–6516`, incl. `cal_coef_Ts3` surface-temp coefficients `TS_A1,TS_B11..13,TS_D1,phase2,Tmax`). For `T<T_0` calls `cal_coef_svsteady_init` (`6544`).

### 1c. CCN critical-supersaturation per aerosol mass bin (`6564–6695`)
Triple loop `n / ica / i=1,N_bin_a`. For `ica/=2`: builds mass-bin boundaries from lognormal params, integrates number/mass fractions via `interp_data1d_lut_big(snrml,...)`, then Köhler activation using **`get_critrad_anal`** and **`get_hazerad_anal`** (`6642–6690`):
```fortran
          r_n=1.0e+4_PS*(M_act(i,ica,n)/N_act(i,ica,n)/ga(ica)%MS(1,n)%den/coef4pi3)**(1.0/3.0)
          sb=nu_aps(ica)*ga(ica)%MS(1,n)%eps_map*M_W*ga(ica)%MS(1,n)%den/(M_aps(ica)*den_w)*phi_aps(ica)
          beta=0.5
          rd_c=1.0e+4_PS*get_critrad_anal(AA/ag%TV(n)%T,sb,beta,r_n*1.0e-4_PS)  ! [micron]
          a_c(i,ica,n)=rd_c*1.0e-4_PS  ! [cm]
          rd_h=get_hazerad_anal(AA/ag%TV(n)%T,sb,beta,1.0+sw(n),r_n*1.0e-4_PS)
          s_c(i,ica,n)=dexp( real(AA/ag%TV(n)%T*1.0e+4_PS/rd_c-&
                          sb*r_n**(2.0*(1.0+beta))/&
                          (rd_c**3-r_n**3),8) ) -1.0d+0
          if(s_c(i,ica,n)<=sw(n).and.r_n<=max_r_n.and.s_c(i,ica,n)<=sw_allow) then
            mean_mass_ap=M_act(i,ica,n)/N_act(i,ica,n)
            mean_mass(i,ica,n)=max(gr%binb(1)*1.05_PS,mean_mass_ap*1.05_PS)
              noccn(i,ica,n)=0
          endif
```
`noccn=0` marks activatable-to-liquid; the "grow 5% by mass" comment explains `mean_mass`.

### 1d. `DHF_IF1` block — Deliquescence-Heterogeneous Freezing precompute (`6707–6908`)
`if(iflg_inuc>0.and.iflg_dhf>0)`. For each aerosol mass bin with `si>0, Tc<0, ica/=2, noccn==1`: computes haze radius, molality, surface tensions (`sigma_iv,sigma_sv,sigma_is`), effective latent heat `Lmef`, `Gn`, `Hvfr`, freezing critical radius `r_cr`, activation energy `dFact`, insoluble radius `r_d`. Then a **40-point Gauss-quadrature θ-PDF integral** (Savre & Ekman 2014) producing `f_gq`, then `N_DHF/M_DHF` scaled by `frac_dust` (`6849–6852`), threshold `Swcr`, and sets `noccn=2` when `Swcr<=1+sw`. `s_c_dhfmin(n)` = min critical S over bins (`6898–6906`).

### 1e. `Dep_IF1` block — Deposition nucleation precompute (`6912–7035`)
`iflg_dep==1` (classical, `dis_type==4` monodisperse): 40-pt Gauss-quadrature over θ giving `f_gq`, `N_dep(i,n)=con*f_gq`, `M_dep=max(gs%binb(1)*1.05,mean_mass_ap*1.05)*N_dep`. `iflg_dep==2`: Meyers scheme `N_dep=min(get_inact(si),con)`. Accumulates `used_Ma_dep(n)`.

### 1f. Activated aerosol totals + akk_lmt (`7039–7093`)
```fortran
          if(noccn(i,ica,n)==0) then
            used_Ma_act(n)=used_Ma_act(n)+max(0.0_PS,N_act(i,ica,n)*mean_mass(i,ica,n)-M_act(i,ica,n))
            used_Na_act(n)=used_Na_act(n)+N_act(i,ica,n)
          elseif(noccn(i,ica,n)==2) then
            ...used_Ma_DHF / used_Na_DHF...
```
`sw_allact`, `ds_allDHF`, `si_alldep` are the "all-activated" saturation offsets (`7055–7079`). CCN number cap (`7085–7091`):
```fortran
        akk_lmt(n)=max(0.0_PS,(min(CCNMAX,used_Na_act(n))-nr_0(n))&
                  /max(1.0e-30_PS,used_Na_act(n)))
```
(The `akk_lmt_DHF` update is present but commented out.)

### 1g. T_il consistency + build effective grid list (`7095–7122`)
```fortran
      til(n)=thil(n)*(ag%TV(n)%P/p00)**(Racp)
      T_a_r=cal_air_temp(Til(n),qr_0(n),qi_0(n))
```
Builds `mbx(1:Lbx)` = boxes with `icycle_n==0 & imethod==0`; copies to `mbx2/Lbx2`.

### 1h. **Vapor-advance stepping — Backward-Euler Til–Qt loop** (`7143–7259`)
This is the primary iterative supersaturation solve:
```fortran
    iter1_loop: do iter=1,ITMAX
      call func_liqvap_vec(used_Mr_vap,used_Mr_act,liq_left,noccnt &
                          ,mbx2,Lbx2,sw_b)!,T_a_b
      call func_icevap_vec(used_Mi_vap,used_Mi_vapliq,used_Mi_act,loss_Mi_mlt,ice_left &
                          ,noindep,nodhft &
                          ,mbx2,Lbx2,sw_b,T_a_b)
      do m=1,Lbx2
        n=mbx2(m)
        trans_Mi(n)=loss_Mi_mlt(n)-min(liq_left(n),gain_Mi_rim(n)+gain_Mi_frn(n))
        qr(n)=max(0.0_PS,qr_0(n)+(used_Mr_act(n)+used_Mr_vap(n)+trans_Mi(n))/ag%TV(n)%den)
        qi(n)=max(0.0_PS,qi_0(n)+(used_Mi_vap(n)+used_Mi_vapliq(n)+used_Mi_act(n)-trans_Mi(n))/ag%TV(n)%den)
        em(n)=0
        if(qr(n)+qi(n)>=qtp(n)) then
          em(n)=1
        endif
        T_a_n(n)=cal_air_temp(Til(n),qr(n),qi(n))
        e_satw=get_sat_vapor_pres_lk(1,T_a_n(n),ag%estbar,ag%esitbar)
        qv_n(n)=max(0.0_PS,qtp(n)-qr(n)-qi(n))
        sw_n(n)=ag%TV(n)%P*qv_n(n)/(Rdvchiarui+qv_n(n))/e_satw-1.0_PS
        e_sati=get_sat_vapor_pres_lk(2,min(T_0,T_a_n(n)),ag%estbar,ag%esitbar)
        r_e=e_satw/e_sati
        si_n(n)=r_e*(sw_n(n)+1.0_PS)-1.0_PS
        dif=abs((qv_n(n)-qv(n))/qv(n))
        difT=abs((T_a_n(n)-T_a_b(n))/T_a_b(n))
        if(&
           (&
           ( (qi(n)>1.0e-20_PS.and.&
             (qi_b(n)<qi(n).and.qi_b2(n)>qi_b(n).and.T_a_b(n)<T_0)).or.&
             (qr(n)>1.0e-20_PS.and.&
             (qr_b(n)<qr(n).and.qr_b2(n)>qr_b(n)))&
             ).and.iter>=10)   .or.&
             (em(n)>=1) &
               )   then
          imethod(n)=1        ! oscillating / over-consumed -> hand to zbrent
          icond1(m)=0
        elseif(dif<=alim1.and.difT<=alim2) then
          sw_n(n)=sw_b(n); T_a_m(n)=T_a_n(n); T_a_n(n)=T_a_b(n)
          qv_n(n)=qtp(n)-qr(n)-qi(n)
          imethod(n)=-1       ! converged
          icond1(m)=0
        else
          sw_b2(n)=sw_b(n); T_a_b2=T_a_b(n)
          sw_b(n)=sw_n(n); si_b(n)=si_n(n); T_a_b(n)=T_a_n(n); qv(n)=qv_n(n)
          qr_b2(n)=qr_b(n); qi_b2(n)=qi_b(n); qr_b(n)=qr(n); qi_b(n)=qi(n)
        endif
      enddo
      ... compact mbx2 to unconverged boxes; exit when none left ...
    enddo iter1_loop
```
Remaining unconverged boxes (`7261–7267`) get `imethod(n)=1`.

### 1i. **Zbrent stage 1** (`7272–7404`)
Builds `mbx` of `imethod==1` boxes, sets brackets, `iphase(n)=2` if pure ice saturation else `1`:
```fortran
    call zbrent_act_vec(sw_n,iphase,1,mbx,Lbx,sw_o,T_a_o,sw_b,T_a_b)
    ... qr_b(n)=qr2(n); qi_b(n)=qi2(n) ...
    call func_vec(sw_m,T_a_m,si_m,fret,0,mbx,Lbx,sw_n,T_a_b)   ! oscillation residual
```
Oscillation check: if `abs(fret(n))>0.1`, reset to initial state, flag `gr%mark_er(n)=1`, re-run `func_vec` on the collected subset.

### 1j. **Zbrent stage 2 — water-saturation adjustment** (`7406–7474`)
For boxes with `qr_0>1e-10` and supersaturation exceeding `sw_allow`, set `iphase=1`, `iqvlmt=2`, and:
```fortran
    call zbrent_act_vec(sw_n,iphase,2,mbx,Lbx,sw,T_a_o,sw_b,T_a_b)
    call func_vec(sw_m,T_a_m,si_m,fret,0,mbx,Lbx,sw_n,T_a_b)
```

### 1k. Placement of activated droplets into liquid bins (`7476–7857`)
`if(sum(noccnt(1:ag%L))==ag%L) goto 5555` (skip if no liquid activation). `category_loop1` over `ica`; for `noccn(i,ica,n)==0` builds shifted bin boundaries `binb3d`, computes the activation fraction `akk` limited by `akk_lmt` and `sw_n/sw_allact` (`7609–7621`):
```fortran
          if(sw_allact(n)>1.0e-25_PS) then
            akk(n)=min(1.0_PS,akk_lmt(n),max(0.0_PS,sw_n(n)/sw_allact(n)))
          else
            akk(n)=min(1.0_PS,max(0.0_PS,akk_lmt(n)))
          endif
          Mp(i,n)=N_act(i,ica,n)*akk(n)*mean_mass(i,ica,n)
          Np(i,n)=N_act(i,ica,n)*akk(n)
```
mass ratios `ratio_Mp(...,rmat_m/rmas_m/rmai_m)`, then:
```fortran
      call add_simple_vec(N_bin_a,ag%L,icond4,gr,Np,Mp,ratio_Mp,new_N,new_M,new_Q,Qp)     ! 7676
      call cal_lincubprms_vec(mxnbin+1,N_bin_a,gr%L,Npd,Mpd,binb3d,a2d,error_number,...)   ! 7682
      ... fallback cal_linprms_vec_s for error bins ...
      call cal_transbin_vec(gr%token,...,new_N,new_M,new_Q,new_mtend,...)                  ! 7745
```
Tendencies `gr%MS(i,n)%dcondt(pro_type)=new_N(i,n)/gr%dt`, `dmassdt(rmt,pro_type)=new_M(i,n,rmt)/gr%dt` (`7796–7842`), with mass-component clamp (`7836–7838`).

### 1l. `5555 continue` → `DHF_IF2` transfer to ice bins (`7859–8069`)
`if(sum(nodhft(1:ag%L))==ag%L) goto 6666`. `DHF_IF2: if(iflg_inuc==1.and.iflg_dhf==1)`. `category_loop2` (skips `ica==2`). For `noccn==2` bins builds frozen-particle geometry `Qp(iacr/iccr/ivcs)`, applies DHF activation fraction:
```fortran
          if(ds_allDHF(n)>1.0e-25_PS) then
            akk_DHF(n)=max(0.0_PS,min(1.0_PS,akk_lmt_DHF(n), &
                          (sw_n(n)-s_c_dhfmin(n))/ds_allDHF(n)))
          else
            akk_DHF(n)=0.0_PS
          endif
          if(akk_DHF(n)>1.0e-6) then
            Mpd(i,n)=M_DHF(i,ica,n)*akk_DHF(n); Npd(i,n)=N_DHF(i,ica,n)*akk_DHF(n)
          endif
```
Finds target ice bin `icond3` and adds `dmassdt/dcondt/dvoldt` under `pro_type_dhf_ice=3` / aerosol sink `pro_type_dhf_aer=6`.

### 1m. `6666 continue` → `DEP_IF2` deposition transfer to ice bins (`8072–8298`)
`if(sum(noindep(1:ag%L))==ag%L) goto 7777`. `DEP_IF2: if(iflg_inuc==1.and.iflg_dep>0)`, `ica=2`. Deposition activation fraction (`8168–8182`):
```fortran
          si=r_e*(sw_n(n)+1.0_PS)-1.0_PS
          if(si_alldep(n)>1.0e-25_PS) then
            akk_dep=min(1.0_PS,si/si_alldep(n))
          else
            akk_dep=0.0_PS
          endif
          if(akk_dep>1.0e-6) then
            Mpd(i,n)=M_dep(i,n)*akk_dep; Npd(i,n)=N_dep(i,n)*akk_dep
          endif
```
Adds ice tendencies under `pro_type_dep_ice=7` / aerosol sink `pro_type_dep_aer=7`.

### 1n. `7777 continue` → finalize thermo (`8301–8402`)
```fortran
        ag%TV(n)%e_sat_n(1)=get_sat_vapor_pres_lk(1,T_a_n(n),...)
        ag%TV(n)%e_sat_n(2)=get_sat_vapor_pres_lk(2,min(T_0,T_a_n(n)),...)
        r_e= ag%TV(n)%e_sat_n(1)/ag%TV(n)%e_sat_n(2)
        ag%TV(n)%T_n=T_a_n(n); ag%TV(n)%T_m=T_a_m(n)
        ag%TV(n)%s_v_n(1)=sw_n(n)
        ag%TV(n)%s_v_n(2)=r_e*(sw_n(n)+1.0_PS)-1.0_PS
```

---

## Contained solver routines (FULL bodies)

### `cal_coef_svsteady_init` (`8406–8414`)
```fortran
    subroutine cal_coef_svsteady_init(g,n,a)
      type (Group)        :: g
      integer :: n!,i
      real(PS) :: a
      a=g%MS(1,n)%coef(1)*g%dt
    end subroutine cal_coef_svsteady_init
```

### `cal_air_temp` (`8416–8431`)
```fortran
    function cal_air_temp(Til,qr,qi) result(T)
      real(PS),intent(in) :: Til,qr,qi
      real(PS) :: T
      T=Til*(1.0_PS+(L_e*qr+L_s*qi)/(c_pa*253.0_PS))
      if(T>253.0_PS) then
         T=0.5*(Til+sqrt(Til**2+4.0_PS*Til/c_pa*(L_e*qr+L_s*qi)))
      endif
    end function cal_air_temp
```

### `func_liqvap_vec` (`8433–8530`)  — liquid mass formed/left from activation + condensational growth
```fortran
    subroutine func_liqvap_vec(used_Mr_vap,used_Mr_act,liq_left,noccnt &
                          ,mbx,Lbx,x)
      integer,intent(in) :: Lbx
      integer,intent(in),dimension(*) :: mbx
      real(PS),dimension(*),intent(in) :: x!,y
      real(PS),dimension(*),intent(inout) :: used_mr_vap,used_Mr_act,liq_left
      integer,dimension(*),intent(inout) :: noccnt
      real(PS) :: d_mean_mass
      integer :: m,n

      do m=1,Lbx
        n=mbx(m)
        used_Mr_vap(n)=0.0_PS
        used_Mr_act(n)=0.0_PS
        liq_left(n)=0.0_PS
        noccnt(n)=1
      enddo

      if(flagp_r>0) then
        do m=1,Lbx
          n=mbx(m)
          if(x(n)>0.0_PS.and.sw(n)>0.0_PS.and.used_Ma_act(n)>1.0e-25_PS) then
            if(sw_allact(n)>1.0e-25_PS) then
              akk(n)=min(1.0_PS,akk_lmt(n),x(n)/sw_allact(n))
            else
              akk(n)=min(1.0_PS,max(0.0_PS,akk_lmt(n)))
            endif
            used_Mr_act(n)=used_Mr_act(n)+akk(n)*used_Ma_act(n)
            noccnt(n)=0
          endif
        enddo
        do j=1,gr%N_BIN
          do m=1,Lbx
            n=mbx(m)
            if((mes_rc(n)==2.or.mes_rc(n)==4).and.&
               gr%MS(j,n)%con>=nlmt.and.gr%MS(j,n)%mass(rmt)>=mlmt ) then
              d_mean_mass=(gr%MS(j,n)%coef(1)*x(n)+gr%MS(j,n)%coef(2))*gr%dt
              if(d_mean_mass*gr%MS(j,n)%con&
                    <=gr%MS(j,n)%mass(rmat)-gr%MS(j,n)%mass(rmt)) then
                used_Mr_vap(n)=used_Mr_vap(n)-(gr%MS(j,n)%mass(rmt)-gr%MS(j,n)%mass(rmat))
                liq_left(n)=liq_left(n)+0.0_PS
              elseif(d_mean_mass<0.0_PS.and.&
                    (gr%MS(j,n)%r_act>&
                  coef3i4p1i3*((gr%MS(j,n)%mean_mass+d_mean_mass)/gr%MS(j,n)%den)**(1.0/3.0)&
                  )) then
                used_Mr_vap(n)=used_Mr_vap(n)-(gr%MS(j,n)%mass(rmt)-gr%MS(j,n)%mass(rmat))
                liq_left(n)=liq_left(n)+0.0_PS
              else
                used_Mr_vap(n)=used_Mr_vap(n)+d_mean_mass*gr%MS(j,n)%con
                liq_left(n)=liq_left(n)+gr%MS(j,n)%mass(rmt)-gr%MS(j,n)%mass(rmat)+d_mean_mass*gr%MS(j,n)%con
              end if
            endif
          enddo
        enddo
      endif
    end subroutine func_liqvap_vec
```
The three branches per bin: (1) full evaporation, (2) shrink below Köhler-critical `r_act` → return to haze (both remove all excess mass, leave `liq_left=0`), (3) normal growth/decay contributing `d_mean_mass*con`.

### `func_icevap_vec` (`8532–8990`) — ice mass, melting, deposition-freezing + DHF
Declaration and init (`8532–8574`):
```fortran
    subroutine func_icevap_vec(used_Mi_vap,used_Mi_vapliq,used_Mi_act,loss_Mi_mlt,ice_left &
                          ,noindep,nodhft &
                          ,mbx,Lbx,x,y)
      use class_Thermo_Var, only: get_sat_vapor_pres_lk
      use class_Ice_Shape, only: get_vip, cal_semiac_ip
      integer,intent(in) :: Lbx
      integer,dimension(*),intent(in) :: mbx
      real(PS),dimension(*),intent(in) :: x,y
      real(PS),dimension(*),intent(inout) :: used_Mi_vap,used_Mi_vapliq,used_Mi_act &
              ,loss_Mi_mlt,ice_left
      integer,dimension(*),intent(inout) :: noindep,nodhft
      real(PS) :: xi,xw,tmp,tmp_bf,d_mean_mass,m_w,m_icore,m_w_c,&
             v_cs_p,semi_aip_p,semi_cip_p,v_w,v_ip,v_space,m_shed,&
             dm_w,mass_ap,dia_icore,q
      real(PS) :: x0,x1,x2,gx0,gx1,gx2,dT_w
      real(8) :: TS_B1,TS_D2,w0,w1,w2,aL,bL,dL
      real(PS) :: flg1
      integer :: m,n

      do m=1,Lbx
        n=mbx(m)
        used_Mi_vap(n)=0.0_PS
        used_Mi_vapliq(n)=0.0_PS
        used_Mi_act(n)=0.0_PS
        loss_Mi_mlt(n)=0.0_PS
        ice_left(n)=0.0_PS
        noindep(n)=1
        nodhft(n)=1
      enddo
```

**`bin_loop1`** (`8577–8914`): for each ice bin, computes over-water/over-ice saturations `xw/xi`, then solves the ice **surface temperature** `tmp` by a locally-fit Lagrange parabola for `esat(Tsfc)/Tsfc`, iterated up to **5 trials** with shrinking windows `dT_w = 20,10,5,1,0.5` (quadratic root `TS_D2` discriminant). The depositional mass change splits on `phase2(j,n)`:

- `phase2==2` (over-ice growth) → `used_Mi_vap` (`8771–8807`)
- else (over-water) → `used_Mi_vapliq` (`8809–8854`)

with `d_mean_mass = (coef(1)*esat/y*(x+1) + coef(2)*esat(tmp)/tmp)*dt`. Melting mass `m_w`, shedding logic via `cal_semiac_ip/get_vip` and `q` polynomial (`8871–8909`), accumulating `loss_Mi_mlt`.

**KC04 deposition-freezing branch** (`8916–8954`) — the `flg1` (Meyers vs classical) selector added by kc04dep:
```fortran
!    Deposition freezing nucleation
!
         if(iflg_inuc>0.and.iflg_dep>0) then
           if(iflg_dep==1) then
             ! classical nucleation theory
             flg1=x(n)
           else
             ! Meyer's secheme
             flg1=-1.0
           endif
           do m=1,Lbx
             n=mbx(m)
             e_satw=get_sat_vapor_pres_lk(1,y(n),ag%estbar,ag%esitbar)
             e_sati=get_sat_vapor_pres_lk(2,min(T_0,y(n)),ag%estbar,ag%esitbar)
             r_e=e_satw/e_sati
             xi=r_e*(x(n)+1.0_PS)-1.0_PS
             if(( xi>0.0_PS.and.y(n) < T_0.and.flg1<0.0_PS) .and. &
                ( ga(2)%MS(1,n)%con>=nlmt.and.ga(2)%MS(1,n)%mass(1)>=mlmt )) then
               if(si_alldep(n)>1.0e-25_PS) then
                 akk_dep=min(1.0_PS,xi/si_alldep(n))
               else
                 akk_dep=0.0_PS
               endif
               used_Mi_act(n)=used_Mi_act(n)+akk_dep*used_Ma_dep(n)
               noindep(n)=0
             end if
           enddo
         endif
```

**KC04 DHF branch** (`8955–8987`) — the `nodhft` path (twin has neither):
```fortran
!    Deliquescence-heterogeneous freezing
!
         if(iflg_inuc>0.and.iflg_dhf>0) then
           do m=1,Lbx
             n=mbx(m)
             if(x(n)>s_c_dhfmin(n)) then
               if(ds_allDHF(n)>1.0e-25_PS) then
                 akk_DHF(n)=max(0.0_PS,min(1.0_PS,akk_lmt_DHF(n),&
                               (x(n)-s_c_dhfmin(n))/ds_allDHF(n)))
               else
                 akk_DHF(n)=0.0_PS
               endif
               if(akk_DHF(n)>1.0e-6) then
                 used_Mi_act(n)=used_Mi_act(n)+akk_DHF(n)*used_Ma_DHF(n)
                 nodhft(n)=0
               endif
             endif
           enddo
         endif
      endif
    end subroutine func_icevap_vec
```

### `zbrent_act_vec` (`8992–9275`) — Brent root-finder (CRITICAL)
```fortran
    subroutine zbrent_act_vec(sw_n,iphase,iswitch,mbx,Lbx,sw_o,T_a_o,sw_b,T_a)
      use class_Thermo_Var, only: get_sat_vapor_pres_lk
      integer,intent(in) :: iswitch
      integer,dimension(*),intent(in) :: iphase
      integer,intent(in) :: Lbx
      integer,dimension(*),intent(in) :: mbx
      real(PS),dimension(*),intent(in) :: sw_o,T_a_o
      real(PS),dimension(*),intent(in) :: sw_b,T_a
      real(PS),dimension(*),intent(inout) :: sw_n
      integer,dimension(LMAX) :: icond2
      integer :: Lbx2
      integer,dimension(LMAX) :: mbx2
      REAL(PS) :: tol=1.0e-6
      real(PS) :: eps
      PARAMETER (EPS=3.0e-8)
      REAL(PS),dimension(LMAX) :: a,b,c,d,e,fa,fb,fc
      REAL(PS) :: p,q,r,s,tol1,xm
      real(PS),dimension(LMAX) :: dum1,dum2,dum3
      integer,parameter :: nb_it=20
      integer :: itr,m,mm,n
      integer,parameter :: ITMAX_ini=30
      integer,parameter :: ITMAX=50
      integer :: ierror1(ag%L)

      if (Lbx == 0) return ! CHIARUI DEBUG

      ! ---- set initial brackets [b,a] per phase ----
      do m=1,Lbx
        n=mbx(m)
        if(iphase(n)==1) then
          b(n)=-0.5_PS
          a(n)=0.2_PS
        elseif(iphase(n)==2) then
          e_satw=get_sat_vapor_pres_lk(1,T_a_o(n),ag%estbar,ag%esitbar)
          e_sati=get_sat_vapor_pres_lk(2,min(T_0,T_a_o(n)),ag%estbar,ag%esitbar)
          r_e=e_satw/e_sati
          b(n)=(-0.5_PS+1.0_PS)/r_e-1.0_PS   ! super saturation over ice is -0.5
          a(n)=0.2_PS
        endif
      enddo

      do m=1,Lbx
        icond2(m)=1
        mbx2(m)=mbx(m)
      enddo
      Lbx2=Lbx

      do n=1,ag%L
        iterz(n)=0
        iterzi(n)=0
        ierror1(n)=0
      enddo

      ! ---- bracket-expansion phase (ITMAX_ini=30) ----
      do itr=1,ITMAX_ini
        call func_vec(dum1,dum2,dum3,fa,iswitch,mbx2,Lbx2,a,T_a)
        call func_vec(dum1,dum2,dum3,fb,iswitch,mbx2,Lbx2,b,T_a)
        do m=1,Lbx2
          n=mbx2(m)
          if(fa(n).lt.0.0_PS.and.fb(n).lt.0.0_PS) then
            if(iswitch==2) then
               a(n)=a(n)*2.0_PS
            else
               a(n)=a(n)+0.2_PS
            endif
            if(a(n)>1.0e+4) then
              ierror1(n)=1
              icond2(m)=0
            endif
          elseif(fa(n).gt.0.0_PS.and.fb(n).gt.0.0_PS) then
            b(n)=b(n)-0.2_PS
            if(b(n)<-2.0e+0) then
              ierror1(n)=2
              icond2(m)=0
            endif
          else
            icond2(m)=0
          end if
        enddo
        if(any(icond2(1:Lbx2)>0)) then
          iterzi(n)=iterzi(n)+1
          mm=0
          do m=1,Lbx2
            if(icond2(m)==1) then
              mm=mm+1
              mbx2(mm)=mbx2(m)
              icond2(mm)=1
            endif
          enddo
          Lbx2=mm
        else
          Lbx2=0
          exit
        endif
      enddo

      do m=1,Lbx2
        n=mbx2(m)
        ierror1(n)=3
      enddo

      if(any(ierror1(1:ag%L)>0)) then
        do m=1,Lbx
          n=mbx(m)
          if(ierror1(n)>0) then
             qv_max=ag%TV(n)%rv+qr_0(n)+qi_0(n)
             if ( debug ) then
                ... diagnostics ...
             end if
             sw_n(n)=sw_o(n)          ! bracketing failed -> fall back to old S
          endif
        enddo
      endif

      ! ---- main Brent iteration (ITMAX=50) ----
      do m=1,Lbx
        icond2(m)=1
        mbx2(m)=mbx(m)
      enddo
      Lbx2=Lbx
      do m=1,Lbx
        n=mbx(m)
        c(n)=b(n)
        fc(n)=fb(n)
      enddo

      do itr=1,ITMAX
        do m=1,Lbx2
          n=mbx2(m)
          if((fb(n).gt.0.0_PS.and.fc(n).gt.0.0_PS).or.(fb(n).lt.0.0_PS.and.fc(n).lt.0.0_PS))then
            c(n)=a(n)
            fc(n)=fa(n)
            d(n)=b(n)-a(n)
            e(n)=d(n)
          endif
          if(abs(fc(n)).lt.abs(fb(n))) then
            a(n)=b(n)
            b(n)=c(n)
            c(n)=a(n)
            fa(n)=fb(n)
            fb(n)=fc(n)
            fc(n)=fa(n)
          endif
          tol1=2.0_PS*EPS*abs(b(n))+0.5_PS*tol
          xm=0.5_PS*(c(n)-b(n))
          if(abs(xm).le.tol1 .or. fb(n).eq.0.)then
            sw_n(n)=b(n)
            icond2(m)=0
          else
            if(abs(e(n)).ge.tol1 .and. abs(fa(n)).gt.abs(fb(n))) then
              s=fb(n)/fa(n)
              if(a(n).eq.c(n)) then
                p=2.0_PS*xm*s
                q=1.0_PS-s
              else
                q=fa(n)/fc(n)
                r=fb(n)/fc(n)
                p=s*(2.*xm*q*(q-r)-(b(n)-a(n))*(r-1.0_PS))
                q=(q-1.0_PS)*(r-1.0_PS)*(s-1.0_PS)
              endif
              if(p.gt.0.) q=-q
              p=abs(p)
              if(2.0_PS*p .lt. min(3.0_PS*xm*q-abs(tol1*q),abs(e(n)*q))) then
                e(n)=d(n)
                d(n)=p/q
              else
                d(n)=xm
                e(n)=d(n)
              endif
            else
              d(n)=xm
              e(n)=d(n)
            endif
            a(n)=b(n)
            fa(n)=fb(n)
            if(abs(d(n)) .gt. tol1) then
              b(n)=b(n)+d(n)
            else
              b(n)=b(n)+sign(tol1,xm)
            endif
            iterz(n)=iterz(n)+1
            sw_n(n)=b(n)
          endif
        enddo
        if(any(icond2(1:Lbx2)>0)) then
          call func_vec(dum1,dum2,dum3,fb,iswitch,mbx2,Lbx2,b,T_a)
          mm=0
          do m=1,Lbx2
            if(icond2(m)==1) then
              mm=mm+1
              mbx2(mm)=mbx2(m)
              icond2(mm)=1
            endif
          enddo
          Lbx2=mm
        else
          Lbx2=0
          exit
        endif
      enddo

      if ( debug ) then
         do m=1,Lbx2
            n=mbx2(m)
            write(*,'("brent exceeding maximum iterations. ...")') ...
         enddo
      end if
     end subroutine zbrent_act_vec
```
Note the two iteration budgets: bracket-expansion `ITMAX_ini=30`, main Brent `ITMAX=50` (distinct from the outer routine's `ITMAX=200`). Bracket expansion: additive `a+=0.2` (or multiplicative `a*=2` when `iswitch==2`), `b-=0.2`, with failure bounds `a>1e4` / `b<-2`. `tol=1e-6`, `EPS=3e-8`, `tol1=2*EPS*|b|+0.5*tol`.

### `func_vec` (`9277–9387`) — residual function driving the solver
```fortran
     subroutine func_vec(x_n,y_n,xi_n,  fa,iswitch,mbx,Lbx,x,y)
      use class_Thermo_Var, only: get_sat_vapor_pres_lk
      real(PS),dimension(*),intent(in) :: x,y
      real(PS),dimension(*),intent(inout) :: x_n,xi_n,y_n,fa
      integer,intent(in) :: iswitch
      integer,intent(in) :: Lbx
      integer,dimension(*),intent(in) :: mbx

      if (Lbx == 0) return ! CHIARUI DEBUG

      call func_liqvap_vec(used_Mr_vap,used_Mr_act,liq_left,noccnt,mbx,Lbx,x)
      call func_icevap_vec(used_Mi_vap,used_Mi_vapliq,used_Mi_act,loss_Mi_mlt,ice_left &
                          ,noindep,nodhft,mbx,Lbx,x,y)
      do m=1,Lbx
        n=mbx(m)
        trans_Mi(n)=loss_Mi_mlt(n)-min(liq_left(n),gain_Mi_rim(n)+gain_Mi_frn(n))
        qr(n)=max(0.0_PS,qr_0(n)+(used_Mr_act(n)+used_Mr_vap(n)+trans_Mi(n))/ag%TV(n)%den)
        qi(n)=max(0.0_PS,qi_0(n)+(used_Mi_vap(n)+used_Mi_vapliq(n)+used_Mi_act(n)-trans_Mi(n))/ag%TV(n)%den)
        if(qr(n)+qi(n)>=qtp(n)) then
          em(n)=iswitch+10
        else
          em(n)=0
        endif
      enddo
      do m=1,Lbx
        n=mbx(m)
        y_n(n)=cal_air_temp(Til(n),qr(n),qi(n))
      enddo
      do m=1,Lbx
        n=mbx(m)
        e_satw=get_sat_vapor_pres_lk(1,y_n(n),ag%estbar,ag%esitbar)
        qv_n(n)=max(0.0_PS,qtp(n)-qr(n)-qi(n))
        x_n(n)=ag%TV(n)%P*qv_n(n)/(Rdvchiarui+qv_n(n))/e_satw-1.0_PS
        e_sati=get_sat_vapor_pres_lk(2,min(T_0,y_n(n)),ag%estbar,ag%esitbar)
        r_e=e_satw/e_sati
        xi_n(n)=r_e*(x_n(n)+1.0_PS)-1.0_PS
      enddo

      select case(iswitch)
      case(1)
        call func_liqvap_vec(used_Mr_vap,used_Mr_act,liq_left,noccnt,mbx,Lbx,x_n)
        call func_icevap_vec(used_Mi_vap,used_Mi_vapliq,used_Mi_act,loss_Mi_mlt,ice_left &
                            ,noindep,nodhft,mbx,Lbx,x_n,y_n)
        do m=1,Lbx
          n=mbx(m)
          trans_Mi(n)=loss_Mi_mlt(n)-min(liq_left(n),gain_Mi_rim(n)+gain_Mi_frn(n))
          qr2(n)=max(0.0_PS,qr_0(n)+(used_Mr_act(n)+used_Mr_vap(n)+trans_Mi(n))/ag%TV(n)%den)
          qi2(n)=max(0.0_PS,qi_0(n)+(used_Mi_vap(n)+used_Mi_vapliq(n)+used_Mi_act(n)-trans_Mi(n))/ag%TV(n)%den)
          fa(n)=(qr(n)+qi(n)-qr2(n)-qi2(n))/max(1.0e-9_RP,ag%TV(n)%rv)
        enddo
      case(2)
        do m=1,Lbx
          n=mbx(m)
          fa(n)=-x_n(n)+sw_allow
        enddo
      case(0)
        do m=1,Lbx
          n=mbx(m)
          fa(n)=(qr(n)+qi(n)-qr_b(n)-qi_b(n))/max(1.0e-9_RP,ag%TV(n)%rv)
        enddo
      case default
      end select
     end subroutine func_vec
  end subroutine cal_aptact_var8_kc04dep
```
Residual semantics: `iswitch=1` → fixed-point self-consistency residual (condensate at guessed S vs. condensate re-evaluated at the resulting S); `iswitch=2` → water-saturation-cap residual `sw_allow - x_n`; `iswitch=0` → oscillation-check residual against `qr_b/qi_b`.

---

## 2. KC04 / `act_type=1` specifics vs the `var8_vec` twin (`cal_aptact_var8_vec`, `3409–5773`)

Verified by structural diff. What kc04dep adds:

- **New arguments** (`5776–5781` vs twin `3409–3411`): `iflg_dhf`, `ap_sig_cp`, `ap_mean_cp`, `CRIC_RN_IMM`, `frac_dust`; `use ... get_hazerad_anal` and `get_inact`.
- **`nodhft` flag**: declared alongside `noindep,noccnt` (`5990`), initialized `nodhft(n)=1` (`6453`), reset per call in `func_icevap_vec` (`8573`), set `nodhft(n)=0` on successful DHF (`8982`), and gates the DHF transfer via `if(sum(nodhft(1:ag%L))==ag%L) goto 6666` (`7861`). The twin has no `nodhft` and no DHF machinery at all.
- **`akk_lmt` / `akk_lmt_DHF` limiting**: `akk_lmt_DHF` (`6020`, init `6439`) is kc04-only. In `func_liqvap_vec` the activation fraction is `akk(n)=min(1.0_PS,akk_lmt(n),x(n)/sw_allact(n))` (`8469`) with the `akk_lmt` CCN cap `min(CCNMAX,used_Na_act)-nr_0` (`7087–7088`). `akk_DHF` uses `min(1.0,akk_lmt_DHF,(x-s_c_dhfmin)/ds_allDHF)` (`8974`).
- **KC04 Köhler**: uses `get_critrad_anal` + `get_hazerad_anal` (Khvorostyanov & Curry analytical, `class_Mass_Bin.F90:2926/3098`); larger `max_r_n=10.0` (twin `max_r_n=1.0`).
- **Deposition branch in `func_icevap_vec`**: twin (`5360–5392`) uses `get_inact_tropic`, `cs1_i_depmode`, and `N_acti`/`used_Mi_act=max(0,cs1_i_depmode*xi,...)*N_acti` directly. kc04dep (`8919–8954`) instead introduces the classical-vs-Meyers selector `flg1` (`flg1=x(n)` if `iflg_dep==1`, else `-1.0`), only fires when `flg1<0.0` (Meyers path inside `func_icevap_vec`; classical CNT is precomputed in `Dep_IF1`), and consumes the precomputed `used_Ma_dep(n)` scaled by `akk_dep=min(1.0,xi/si_alldep(n))`. The condition also differs: kc04 uses `xi>0 .and. y(n)<T_0 .and. flg1<0` whereas twin uses `xi>0 .and. ((y<=T_0-5 .and. x>0) .or. (xi>=0.05 .and. y<T_0))`.
- **New process-type constants** (`6001–6004`): `pro_type_dhf_ice=3`, `pro_type_dhf_aer=6`, `pro_type_dep_ice=7`, `pro_type_dep_aer=7`; plus 40-point Gauss quadrature tables, `INMAX`, `max_r_n_DHF`, and the whole Savre & Ekman θ-PDF freezing physics (`DHF_IF1`, `Dep_IF1`).

---

## 3. Placing activated droplets into liquid bins

`add_simple_vec` (`mod_amps_core.F90:15452–15577`) — called at line 7676 for the "mean mass not constrained" case:
```fortran
  subroutine add_simple_vec(nbin,L,icond,g,Np,Mp,ratio_Mp,new_N,new_M,new_Q,Qp)
    ...
    ! locate target bin jbin by mean mass Mp/Np
    do n = 1, L
    do i = 1, nbin
      jbin(i,n)=-1
      if(icond(i,n)==1) then
        if(g%binb(g%N_binb)<Mp(i,n)/Np(i,n)) then
          jbin(i,n)=g%N_BIN
        elseif(g%binb(1)>Mp(i,n)/Np(i,n)) then
          jbin(i,n)=1
        else
          jbin(i,n)=0
        endif
      endif
    enddo
    enddo
    do j=1,g%N_BIN
       do n = 1, L
       do i = 1, nbin
        if(jbin(i,n)==0) then
          if(g%binb(j)<Mp(i,n)/Np(i,n).and.Mp(i,n)/Np(i,n)<=g%binb(j+1)) then
            jbin(i,n)=j
          end if
        endif
      enddo
      enddo
    enddo
    do n = 1, L
    do i= 1, nbin
      if(icond(i,n)==1) then
        j=jbin(i,n)
        new_M(j,n,1)=new_M(j,n,1)+Mp(i,n)
        Np(i,n)=max(Mp(i,n)/0.99_PS/g%binb(j+1),min(Mp(i,n)/1.01_PS/g%binb(j),Np(i,n)))
        new_N(j,n)=new_N(j,n)+Np(i,n)
      endif
    enddo
    enddo
    ! mass components: new_M(j,n,1+jj) += Mp*ratio_Mp(i,n,jj)
    ! non-mass: token-dependent (activated-IN concentration, ivcs/inex, iacr/iccr/... cubed)
```
Number is clamped so the placed mean mass stays inside `[g%binb(j), g%binb(j+1)]`. `add_samebin_vec` (`15579–15657`) is the sibling used elsewhere. In the production path the liquid-activation placement uses `add_simple_vec` (7676) followed by `cal_lincubprms_vec`/`cal_transbin_vec`; `add_samebin_vec` is not called from `cal_aptact_var8_kc04dep`.

**Köhler critical/haze radius usage** — `get_critrad_anal` (`class_Mass_Bin.F90:2926–2954`, KC 2014 eq 6.3.7):
```fortran
  Zd=r_n**beta*sqrt(sb/(3.0_PS*AA))
  Zdc=Zd
  if(Zd>1.0e-2) then
    Pp=(Zdc*Zdc*Zdc+sqrt(Zdc*Zdc*Zdc+0.25_PS)+0.5_PS)**(1.0/3.0)
    Pm=(Zdc*Zdc*Zdc-sqrt(Zdc*Zdc*Zdc+0.25_PS)+0.5_PS)**(1.0/3.0)
    r_cr=r_n*(Zdc+Pp+Pm)
  else
    r_cr=r_n*(1.0d+0+Zdc+2.0d+0*Zdc*Zdc*Zdc/3.0d+0)
  endif
```
`get_hazerad_anal` (`3098–3158`, KC 2014 eq 6.2.12/6.2.28) has three regimes: `S_w<0.97` closed form, small-`Zd` series, and the `S_w≈1` complex-root branch. `get_critrad_itr` (`2794–2924`) is the iterative fallback — present but **not** used on the kc04dep production path (the analytic variants are used at lines 6651/6655/6743).

---

## 4. `iflg_dhf` (micexfg 19) branches — verbatim

- **`DHF_IF1`** (precompute, `6707–6908`): `DHF_IF1: if(iflg_inuc>0.and.iflg_dhf>0) then ... endif DHF_IF1` — computes `N_DHF/M_DHF`, `s_c_dhf`, `s_c_dhfmin`, sets `noccn=2`. (Full physics quoted in §1d above.)
- **`DHF_IF2`** (transfer to ice, `7863–8061`): `DHF_IF2: if(iflg_inuc==1.and.iflg_dhf==1) then ... endif DHF_IF2` — reached after `if(sum(nodhft(1:ag%L))==ag%L) goto 6666`; applies `akk_DHF` and deposits frozen mass/number into ice bins under `pro_type_dhf_ice`/`pro_type_dhf_aer`. (Full block quoted in §1l.)
- The DHF sink inside `func_icevap_vec` (`8958–8986`) is quoted in full in §1/`func_icevap_vec` above (`if(iflg_inuc>0.and.iflg_dhf>0)` → sets `nodhft(n)=0`).

The `iflg_dhf` argument is `intent(in)` at line 5819. The twin `cal_aptact_var8_vec` does not take `iflg_dhf` and contains none of `DHF_IF1`/`DHF_IF2`/the DHF sink.

Scratch diff artifacts written to `/private/tmp/claude-503/-Users-jcanton-projects-scale-amps/cfc898c5-fd11-4ed1-82d3-c7210332e211/scratchpad/{twin.f90,kc04.f90}`.