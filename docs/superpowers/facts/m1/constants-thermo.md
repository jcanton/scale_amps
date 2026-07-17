# AMPS Constants + Thermodynamics — Verbatim Fortran Extraction

All units CGS unless noted. `PS`/`MP_KIND`/`RP` = SCALE's `RP` real kind; `DS` = double precision.

---

## 1. `mod_amps_const.F90` — FULL FILE (lines 1–62)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_const.F90:1-62`

```fortran
module mod_amps_const
  use acc_amps
  implicit none

  public

  ! reference pressure (g/cm/s^2)
  real(PS), parameter :: p00 = 1.0e6_PS

  ! temperature at triple point
  real(PS), parameter :: T_0 = 273.16_PS

  ! gravity (cm/s^2)
  real(PS), parameter :: gg = 980.0_PS

  ! gas constant for vapor (ergs/deg/g)
  real(PS), parameter :: M_w = 18.016_PS

  real(PS), parameter :: R_u = 8.31436e7_PS
  real(PS), parameter :: MR = M_w / R_u

  ! density of water (g/cm^3)
  real(PS), parameter :: den_w = 1.0_PS
  ! density of ice (bulk density) (g/cm^3)
  real(PS), parameter :: den_i = 0.91668_PS

  ! gas constant for dry air (ergs/deg/g)
  real(PS), parameter :: R_d = 287.04e4_PS
  ! gas constant for vapor (ergs/deg/g)
  real(PS), parameter :: R_v = 461.5e4_PS

  ! heat capacity of dry air
  real(PS), parameter :: C_pa = 1004.64e4_PS

  real(PS), parameter :: Racp = R_d / C_pa
  real(PS), parameter :: Rdvchiarui = R_d / R_v
  real(PS), parameter :: M_a = M_w / Rdvchiarui



  ! heat content of water
  real(PS), parameter :: c_w = 4.187e5_PS

  ! latent heat of condensation in ergs/g
  real(PS), parameter :: L_e = 2.5e10_PS
  ! latent heat of freeze in ergs/g
  real(PS), parameter :: L_f = 0.3337e10_PS
  ! latent heat of sublimation in ergs/g
  real(PS), parameter :: L_s = 2.8337e10_PS

  ! condensation coefficient  for liquid
  real(PS), parameter :: a_cliq = 0.036_PS


  ! thermal conductivity
  real(PS), parameter :: k_w = 0.58e5_PS


  real(PS), parameter :: undef = 999.9e30_PS
  integer,  parameter :: iundef = 999

end module mod_amps_const
```

---

## 2. `acc_amps.F90` — precomputed constants (lines 21–34; kinds at 10–19)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/acc_amps.F90:10-34`

```fortran
  integer, parameter  :: DS = DP
  ! just real(4)
  !integer, parameter  :: PS = SELECTED_REAL_KIND(6, 37)
  !integer, parameter  :: PS = SELECTED_REAL_KIND(15, 307)

  integer, parameter  :: PS = RP_SCALE
  integer, parameter  :: MP_KIND = PS
  integer, parameter  :: PS_KIND = PS

  integer, parameter  :: RP = PS

  !---------------------------------------------------------------------------------------------
  ! precalculated coefficients
  real(MP_KIND), parameter :: PI = 3.141592653589793238462643_PS
  real(MP_KIND), parameter :: sq_three = 1.7320508075688772935_PS
  ! 3 sqrt(3), 4pi/3, 2.0*pi
  real(MP_KIND),parameter :: coef3s=5.196152423_MP_KIND, coef4pi3=4.18879020478639_MP_KIND,coef2p=6.28318530717959_MP_KIND
  ! pi/6.0, sqrt(2*pi), 3*sqrt(3)
  real(MP_KIND),parameter :: coefpi6=0.523598776_MP_KIND,coefsq2p=2.506628274631_MP_KIND,coef3sq3=5.19615242270663_MP_KIND
  ! 4.0*pi,(3.0/4.0/pi)**(1.0/3.0)
  real(MP_KIND),parameter :: coef4p=1.25663706144E+01_MP_KIND,coef3i4p1i3=0.62035049089940_MP_KIND
  ! pi/6.0, sqrt(2*pi), 3*sqrt(3)
  real(DS),parameter :: coedpi6=0.523598775598299_DS,coedsq2p=2.506628274631_DS,coed3sq3=5.19615242270663_DS
  ! pi/180.0
  real(MP_KIND),parameter :: coefpi180=0.0174532925199433_MP_KIND
```

---

## 3. Thermodynamics functions

### 3a. e_sat lookup table generation — `QSPARM2` (Murphy & Koop 2005)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_utility.F90:3012-3033` — full subroutine. Index convention: `ESTBAR(K)`/`ESITBAR(K)` hold e_sat at T = 164+(K-1) K, i.e. K=1 → 164 K. Tables are in **hPa** (mb); the accessor multiplies by 10 to get CGS (g/s²/cm = dyn/cm² × 10⁻¹... precisely: hPa × 10 = dyn/cm²? note hPa = 1000 dyn/cm², and lookup multiplies by 10 — the MK formulas give Pa, and Pa × 10 = dyn/cm² (g/s²/cm), so tables are in **Pa**).

```fortran
      SUBROUTINE QSPARM2(ESTBAR,ESITBAR)
!        based on Murphy and Koop (2005)
      implicit none
      real(DS), DIMENSION(150) ::  ESTBAR
      real(DS), DIMENSION(111) ::  ESITBAR
      real(DS) :: T
      integer :: K, JD
      T=163.
      DO K=1,111
        T=T+1.
        ESITBAR(K)=exp(9.550426-5723.265/T+3.53068*log(T)-0.00728332*T)
      ENDDO

      T=163.
      DO JD=1,150
        T=T+1.
        ESTBAR(JD)=exp(54.842763-6763.22/T-4.210*log(T)+0.000367*T &
                  + tanh(0.0415*(T-218.8))*(53.878-1331.22/T &
                  - 9.44523*log(T)+0.014025*T))
      end do
      RETURN
      END subroutine QSPARM2
```

### 3b. `get_sat_vapor_pres` (analytic, Lowe & Ficke 1974 style)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:462-494`

```fortran
  function get_sat_vapor_pres(phase, d_T) result (e_sat)
    ! **********************************************************************
    ! Calculate the saturation vapor pressure over water and ice
    ! according to Lowe and Ficke (1974)
    ! **********************************************************************
    ! ambient temperature (K)
    real (PS), intent(in) :: d_T
    ! phase of water which the vapor pressure exist over
    !    water => 1   ice => 2 
    integer          :: phase
    ! saturation vapor pressure ( g/s^2/cm )
    real (PS)                   :: e_sat
    ! coefficient
    real (PS), dimension(7)     :: a

    if( phase == 1 ) then
       ! over liquid
       e_sat = 6.1070_PS*exp(17.15_PS*(d_T-273.16_PS)/(d_T-38.25_PS))
    else if( phase == 2) then
       ! over ice
       e_sat = 6.1064_PS*exp(21.88_PS*(d_T-273.16_PS)/(d_T-7.65_PS))
    end if

    e_sat = 1000.0_PS*e_sat ! g/s^2/cm
    
  end function get_sat_vapor_pres
```

(Comment lines with polynomial coefficient arrays `a(:)` at 481–482, 486–487 are commented out in source.)

### 3c. `get_sat_vapor_pres_lk` — lookup accessor

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:496-528`

```fortran
  function get_sat_vapor_pres_lk(phase, T, estbar, esitbar ) result(e_sat)
    !_______________________________________________________________________
    !     THIS ROUTINE COMPUTES SATURATION MIXING RATIO OVER ICE
    !     WATER BASED ON A TABLE LOOK UP PROCEEDURE DEVELOPED BY
    !     DERICKSON AND COTTON (1977) FOR THE LIQUID PHASE
    !_______________________________________________________________________
    implicit none
    integer,intent(in)  :: phase
    real(PS),intent(in) :: T
    real(DS),intent(in) :: estbar(150),esitbar(111)
    integer :: I,J
    real :: wt
    real(PS) :: e_sat
    ! g/s^2/cm
    if( phase==1 ) then
       I=MAX(1,MIN(int(T)-163,149))
       wt=MAX(MIN(T-real(I+163,PS_KIND),1.0_RP),0.0_RP)
       e_sat=(estbar(I)*(1.0_RP-wt)+estbar(I+1)*wt)*10.0_RP
    else
       ! phase == 2
       J=MAX(1,MIN(int(T)-163,110))
       wt=MAX(MIN(T-real(J+163,PS_KIND),1.0_RP),0.0_RP)
       e_sat=(esitbar(J)*(1.0_RP-wt)+esitbar(J+1)*wt)*10.0_RP
    end if
  end function get_sat_vapor_pres_lk
```

### 3d. `get_T_fesv_lk` — reverse lookup (T from e_sat)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:620-689`

```fortran
  function get_T_fesv_lk(phase, T, des,estbar, esitbar ) result(T_n)
    implicit none
    integer,intent(in)  :: phase
    real(PS),intent(in) :: T
    real(DS),intent(in) :: estbar(150),esitbar(111)
    integer :: I,J,k
    real(PS) :: wt,des,es,T_n
    es=des/10.0_PS
    ! g/s^2/cm
    if( phase==1 ) then
       I=MAX(1,MIN(int(T)-163,149))
       wt=MAX(MIN(T-real(I+163,PS_KIND),1.0_RP),0.0_RP)

       do k=max(1,I-5),min(149,I+5)
          if(estbar(k)<=es.and.estbar(k+1)>es) then
             wt=(es-estbar(k))/(estbar(k+1)-estbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
       do k=1,I-4
          if(estbar(k)<=es.and.estbar(k+1)>es) then
             wt=(es-estbar(k))/(estbar(k+1)-estbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
       do k=I+5,149
          if(estbar(k)<=es.and.estbar(k+1)>es) then
             wt=(es-estbar(k))/(estbar(k+1)-estbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
    else if( phase==2 ) then
       J=MAX(1,MIN(int(T)-163,110))
       wt=MAX(MIN(T-real(J+163,PS_KIND),1.0_RP),0.0_RP)
       do k=max(1,J-5),min(110,J+5)
          if(esitbar(k)<=es.and.esitbar(k+1)>es) then
             wt=(es-esitbar(k))/(esitbar(k+1)-esitbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
       do k=1,J-4
          if(esitbar(k)<=es.and.esitbar(k+1)>es) then
             wt=(es-esitbar(k))/(esitbar(k+1)-esitbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
       do k=J+5,110
          if(esitbar(k)<=es.and.esitbar(k+1)>es) then
             wt=(es-esitbar(k))/(esitbar(k+1)-esitbar(k))
             T_n=real(k+163,PS_KIND)+wt
             goto 10
          end if
       end do
    end if
10 continue
    
  end function get_T_fesv_lk
```

(Note: `des` is declared in the local `real(PS)` list, not with `intent(in)`; if no bracket matches, `T_n` is returned undefined.)

### 3e. `get_diffusivity` — D_v (Hall & Pruppacher 1976), cm²/s

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:381-392`

```fortran
  function get_diffusivity(P,T) result (D_v)
    ! **********************************************************************
    ! Calculate the diffusivity of water vapor in air for temperatures
    ! between -40 to 40 C according to Hall and Pruppacher (1976)
    ! **********************************************************************
    ! ambient pressure and temperature in g/s^2/cm and Kelvin
    real (PS), intent(in)   :: P, T
    ! diffusivity of water vapor in cm^2/sec
    real (PS)               :: D_v
    real (PS), parameter:: P_0 = 1013250.0_PS ! g/s^2/cm
    D_v = 0.211_PS * ((T/T_0)**1.94_PS) * (P_0/P)
  end function get_diffusivity
```

### 3f. `get_thermal_conductivity` — k_a (Beard & Pruppacher 1971a)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:450-460`

```fortran
  function get_thermal_conductivity(T) result (k_a)
    ! **********************************************************************
    ! Calculate the thermal conductivity of dry air
    ! according to Beard and Pruppacher (1971a)
    ! **********************************************************************
    ! ambient temperature in Kelvin
    real (PS), intent(in)   :: T
    ! thermal conductivity in 
    real (PS)               :: k_a
    k_a = (5.69 + 0.017*(T-273.16))*4.1868*1.0e+02 ! gcm^2/s^2 /cm/sec/deg
  end function get_thermal_conductivity
```

### 3g. `get_dynamic_viscosity` — d_vis, g/cm/s (Pruppacher & Klett)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:530-544`

```fortran
  function get_dynamic_viscosity( d_T) result (mu)
    ! calculate the dynamic viscosity, g/cm/s
    ! ambient temperature (K)
    real (PS), intent(in) :: d_T
    real (PS)             :: mu,Tc

    Tc=d_T-273.15_PS
    ! from Pruppacker's book 
    if(Tc>=0.0_PS) then
       mu=(1.718_PS+0.0049_PS*Tc)*1.0e-4_PS
    else
       mu=(1.718_PS+0.0049_PS*Tc-1.2e-5_PS*Tc*Tc)*1.0e-4_PS
    end if
  end function get_dynamic_viscosity
```

### 3h. `get_sfc_tension` — sig_wa, erg/cm² (P&K97 eq. 5-12)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:545-559`

```fortran
  function get_sfc_tension( d_T) result (sig_wa)
    ! calculate the surface tension of water: erg/cm^2
    real (PS), intent(in) :: d_T
    real (PS)             :: sig_wa,Tc
    ! (5-12) of Prupacker and Klett 97
    real(PS),parameter,dimension(1:7) :: &
        an=(/75.93,0.115,6.818e-2,6.511e-3,2.933e-4,6.283e-6,5.285e-8/)
    ! valid between 40 and -40 Celisus
    Tc=max(-45.0_RP,min(40.0_RP,d_T-273.15_RP))
    sig_wa=an(1)+an(2)*Tc+an(3)*Tc*Tc+an(4)*Tc*Tc*Tc+an(5)*Tc*Tc*Tc*Tc+&
           an(6)*Tc*Tc*Tc*Tc*Tc+an(7)*Tc*Tc*Tc*Tc*Tc*Tc
  end function get_sfc_tension
```

### 3i. `get_GTP` — growth (diffusivity) function

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:561-585`

```fortran
  function get_GTP( d_T, d_P, d_D_v, d_k_a, d_es, iphase) result (gtp)
    ! calculate the diffusivity function.
    ! ambient temperature (K), pressure (g/s^2/cm) 
    real (PS), intent(in)            :: d_T, d_P
    ! diffusivity of water vapor in air
    real (PS), intent(in)            :: d_D_v
    ! thermal conductivity of dry air
    real (PS), intent(in)            :: d_k_a
    ! saturation vapor pressure ove water or ice.
    real (PS), intent(in)            :: d_es
    ! phase of water, 1: liquid, 2: solid
    integer, intent(in)              :: iphase
    !
    ! diffusivity function
    real (PS)                        :: gtp
    !

    if( iphase == 1 ) then
       gtp = R_v*d_T/(d_es*d_D_v)+L_e*( L_e/(R_v*d_T) - 1.0_PS) /(d_k_a*d_T)
       gtp = 1.0_PS/gtp
    else if( iphase == 2 ) then
       gtp = R_v*d_T/(d_es*d_D_v)+L_s*( L_s/(R_v*d_T) - 1.0_PS) /(d_k_a*d_T)
       gtp = 1.0_PS/gtp
    endif
  end function get_GTP
```

### 3j. MGTP — modified diffusivity / conductivity (kinetic corrections) + `update_DvKa`

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:393-448` and `691-705`

```fortran
  function get_mod_diffusivity(iphase,radius,OD_v,T) result (out)
    ! **********************************************************************
    ! Calculate the modified diffusivity of water vapor in air
    ! **********************************************************************
    integer,intent(in) :: iphase
    real (PS), intent(in) :: radius
    ! ambient pressure and temperature in g/s^2/cm and Kelvin
    real (PS), intent(in)   :: T
    ! original diffusivity of water vapor in cm^2/sec
    real (PS), intent(in)               :: OD_v
    ! modified diffusivity of water vapor in cm^2/sec
    real (PS)                :: out
    ! condensation coefficient
    real(PS) :: a_c
    ! condensation coefficient  for liquid
    real (PS), parameter             :: a_cliq=1.0
!test    real (PS), parameter             :: a_cliq=0.036
    ! deposition coefficient  for ice
    real (PS), parameter             :: a_cice1=0.5
    ! deposition coefficient  for ice for T<-40C
    real (PS), parameter             :: a_cice2=0.006
    ! vapor jump length
    real (PS), parameter             :: del_v=1.0e-5

    if(iphase==1) then
       ! liquid
       a_c=a_cliq
    else
       !if(iphase==2) then
       ! ice
       a_c=a_cice1
    endif
    out=radius/(radius+del_v)+sqrt(2.0*PI/(R_v*T))*OD_v/(radius*a_c)
    out=OD_v/out
  end function get_mod_diffusivity

  function get_mod_thermal_cond(radius,OK_a,T,den) result (out)
    real (PS), intent(in) :: radius
    ! ambient pressure and temperature in g/s^2/cm and Kelvin
    real (PS), intent(in)   :: T,den
    ! original thermal diffusivity of water vapor in cm^2/sec
    real (PS), intent(in)               :: OK_a
    ! modified thermal diffusivity of water vapor in cm^2/sec
    real (PS)                :: out
    ! thermal accommodation coefficient
    real (PS), parameter             :: a_t=0.96
    ! vapor jump length
    real (PS), parameter             :: del_t=2.16e-5

    out=radius/(radius+del_t)+sqrt(2.0*PI*M_a/(R_u*T))*OK_a/(den*radius*a_t*c_pa)
    out=OK_a/out
  end function get_mod_thermal_cond
```

```fortran
  subroutine update_DvKa(th_var,radius)
    type (Thermo_Var)  :: th_var
    real(PS), intent(in)   :: radius ! precision test (MP_KIND)
    integer :: i
    real (PS) :: mD_v,mk_a
    ! calculate modified diffusivity and conducitivity
    mk_a = get_mod_thermal_cond(radius,th_var%K_a,th_var%T,th_var%den)
    
    do i=1, 2
       mD_v = get_mod_diffusivity(i,radius,th_var%D_v,th_var%T)
       th_var%MGTP(i) = get_GTP( th_var%T, th_var%P, mD_v, mk_a, &
         th_var%e_sat(i), i)
    end do

  end subroutine update_DvKa
```

### 3k. Helpers: `get_rv`, `renew_rv_var`, `get_fkn`

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:588-618, 707-730`

```fortran
  function get_rv( d_T, d_e, d_den) result (rv)
    ! calculate the diffusivity function.
    ! ambient temperature (K)
    real (PS), intent(in)            :: d_T
    ! vapor pressure ove water or ice.
    real (PS), intent(in)            :: d_e
    ! density of ambient air
    real (PS), intent(in)            :: d_den
    !
    ! mixing ratio of vapor
    real (PS)                        :: rv
    !

    rv = d_e/(R_v*d_T)/d_den
  end function get_rv

  subroutine renew_rv_var(th_var, rv)
    type (Thermo_Var)  :: th_var
    real(PS), intent(in)    :: rv
    integer   :: i
    !
    th_var%rv = rv
    th_var%e = (th_var%rv*th_var%den*R_v)*th_var%T

    ! for vectorization purpose.
    th_var%s_v(1) = th_var%e/th_var%e_sat(1) - 1.0_PS
    th_var%s_v_n(1)=th_var%s_v(1)
    th_var%s_v(2) = th_var%e/th_var%e_sat(2) - 1.0_PS
    th_var%s_v_n(2)=th_var%s_v(2)

  end subroutine renew_rv_var
```

```fortran
  function get_fkn(th_var,phase,r0) result(fkn)
    ! NOTE:
    ! this should be the same as get_mod_diffusivity.
    !
    type (Thermo_Var)  :: th_var
    real(PS), intent(in)   :: r0
    integer :: phase
    real(PS) :: beta
    ! deposition coefficient, thickness of a bc
    real(PS),parameter :: beta_w=0.036
    real(PS),parameter :: beta_i1=0.5,beta_i2=0.006
    real(PS),parameter :: delta=1.0e-5    
    real(PS) :: fkn

    if(phase==1) then
       beta=beta_w
    else
       beta=beta_i1
    end if
    ! +++ calculate kinetic effect +++
    fkn = r0/(r0+delta)+(th_var%D_v/beta/r0)*&
         sqrt(2.0_PS*PI/R_v/th_var%T)
    fkn = 1.0_PS/fkn
  end function get_fkn
```

---

## 4. Thermo_Var update sequence

### 4a. Canonical constructor `make_Thermo_Var3` (MKS inputs → CGS)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:249-308`

```fortran
  function make_Thermo_Var3 (drv, dden, pres, temp, dw,destbar,desitbar ) result (th_var)
    use mod_amps_utility, only: get_growth_mode_max
    real(MP_KIND), intent(in)   :: drv, dden, pres, temp, dw
    real(DS), intent(in)   :: destbar(150),desitbar(111)
    type (Thermo_Var)  :: th_var
    integer            :: i
    real(PS) :: check

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    ! Input variables are in MKS unit.
    ! They has to be changed to CGS unit.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    th_var%rv = drv
    th_var%den = dden*1.0e-3_PS
    th_var%P = pres*1.0e+1_PS
    th_var%T = temp
    th_var%W = dw*1.0e+2_PS

    th_var%D_v = get_diffusivity(th_var%P,th_var%T)
    th_var%k_a = get_thermal_conductivity(th_var%T)
    th_var%d_vis = get_dynamic_viscosity( th_var%T)
    th_var%sig_wa=get_sfc_tension(th_var%T)
    th_var%e=th_var%P*th_var%rv/(Rdvchiarui+th_var%rv)

    th_var%den_a=M_a*(th_var%P-th_var%e)/(R_u*th_var%T)

    do i=1,2
       th_var%e_sat(i) = get_sat_vapor_pres_lk(i, th_var%T, destbar, desitbar )
       th_var%s_v(i) = th_var%e/th_var%e_sat(i) - 1.0_PS
       th_var%GTP(i) = get_GTP( th_var%T, th_var%P, th_var%D_v, th_var%k_a, &
         th_var%e_sat(i), i)
       th_var%rv_sat(i)=Rdvchiarui*th_var%e_sat(i)/max((th_var%P-th_var%e_sat(i)),th_var%e_sat(i))

       if(th_var%rv_sat(i)<=th_var%rv) then
          th_var%s_v(i)=max(0.0_PS,th_var%s_v(i))
       else
          th_var%s_v(i)=min(0.0_PS,th_var%s_v(i))
       end if

       th_var%s_v_n(i)=th_var%s_v(i)
    end do
    th_var%dmassdt=0.0_PS
    th_var%dmassdt_v=0.0_PS

    th_var%e_sat(2)=min(th_var%e_sat(1),th_var%e_sat(2))
    th_var%rv_sat(2)=min(th_var%rv_sat(1),th_var%rv_sat(2))

    th_var%nuc_gmode=get_growth_mode_max(th_var%T,th_var%s_v(2))    

    th_var%e_sat_n=th_var%e_sat
    th_var%T_n=th_var%T
    th_var%T_m=th_var%T
  end function make_Thermo_Var3
```

### 4b. `ini_AirGroup` — the routine that actually fills Thermo_Var from (RV, DEN, PT, T, W)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_AirGroup.F90:156-262` (core fill loop, lines 191–249; inputs MKS: DEN kg/m³, PT Pa, W m/s):

```fortran
    do n=1,a%L
!       a%TV(n) = make_thermo_var3( RV(n),DEN(n),PT(n),T(n),W(n),a%estbar,a%esitbar )
      a%TV(n)%rv = rv(n)
      a%TV(n)%den = den(n)*1.0e-3_PS
      a%TV(n)%P = PT(n)*1.0e+1_PS
      a%TV(n)%T = T(n)
      a%TV(n)%W = W(n)*1.0e+2_PS

      a%TV(n)%D_v = get_diffusivity(a%TV(n)%P,a%TV(n)%T)
      a%TV(n)%k_a = get_thermal_conductivity(a%TV(n)%T)
      a%TV(n)%d_vis = get_dynamic_viscosity( a%TV(n)%T)
      a%TV(n)%sig_wa=get_sfc_tension(a%TV(n)%T)
      a%TV(n)%e=a%TV(n)%P*a%TV(n)%rv/(Rdvchiarui+a%TV(n)%rv)

      a%TV(n)%den_a=M_a*(a%TV(n)%P-a%TV(n)%e)/(R_u*a%TV(n)%T)

      a%TV(n)%dmassdt=0.0_PS

      a%TV(n)%e_sat(1) = get_sat_vapor_pres_lk(1, a%TV(n)%T, a%estbar, a%esitbar )
      a%TV(n)%e_sat(2) = get_sat_vapor_pres_lk(2, a%TV(n)%T, a%estbar, a%esitbar )

      a%TV(n)%s_v(1) = a%TV(n)%e/a%TV(n)%e_sat(1) - 1.0_PS
      a%TV(n)%GTP(1) = get_GTP( a%TV(n)%T, a%TV(n)%P, a%TV(n)%D_v, a%TV(n)%k_a, &
              a%TV(n)%e_sat(1), 1)
      a%TV(n)%rv_sat(1)=Rdvchiarui*a%TV(n)%e_sat(1)/max((a%TV(n)%P-a%TV(n)%e_sat(1)),a%TV(n)%e_sat(1))

      a%TV(n)%s_v(2) = a%TV(n)%e/a%TV(n)%e_sat(2) - 1.0_PS
      a%TV(n)%GTP(2) = get_GTP( a%TV(n)%T, a%TV(n)%P, a%TV(n)%D_v, a%TV(n)%k_a, &
              a%TV(n)%e_sat(2), 2)
      a%TV(n)%rv_sat(2)=Rdvchiarui*a%TV(n)%e_sat(2)/max((a%TV(n)%P-a%TV(n)%e_sat(2)),a%TV(n)%e_sat(2))

      if(a%TV(n)%rv_sat(1)<=a%TV(n)%rv) then
        a%TV(n)%s_v(1)=max(0.0_PS,a%TV(n)%s_v(1))
      else
        a%TV(n)%s_v(1)=min(0.0_PS,a%TV(n)%s_v(1))
      end if
      if(a%TV(n)%rv_sat(2)<=a%TV(n)%rv) then
        a%TV(n)%s_v(2)=max(0.0_PS,a%TV(n)%s_v(2))
      else
        a%TV(n)%s_v(2)=min(0.0_PS,a%TV(n)%s_v(2))
      end if

      a%TV(n)%dmassdt_v(1)=0.0_PS
      a%TV(n)%dmassdt_v(2)=0.0_PS
      a%TV(n)%dmassdt_v(3)=0.0_PS

      a%TV(n)%T_n=a%TV(n)%T
      a%TV(n)%T_m=a%TV(n)%T

      a%TV(n)%e_sat(2)=min(a%TV(n)%e_sat(1),a%TV(n)%e_sat(2))
      a%TV(n)%rv_sat(2)=min(a%TV(n)%rv_sat(1),a%TV(n)%rv_sat(2))

      a%TV(n)%s_v_n(1)=a%TV(n)%s_v(1)
      a%TV(n)%e_sat_n(1)=a%TV(n)%e_sat(1)
      a%TV(n)%s_v_n(2)=a%TV(n)%s_v(2)
      a%TV(n)%e_sat_n(2)=a%TV(n)%e_sat(2)
    enddo
```

followed by (class_AirGroup.F90:253-260):

```fortran
    call cal_growth_mode_inl_vec(igm,1,a%L,ihabit_gm_random &
                          ,a%TV(1:a%L)%T,a%TV(1:a%L)%s_v(2),rdsd)

    do n=1,a%L
      a%TV(n)%nuc_gmode=igm(n)
    enddo
```

### 4c. `update_AirGroup` — per-step refresh from (RV, T) keeping P, den

`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_AirGroup.F90:264-338`. Identical fill sequence except: only `rv` and `T` are reassigned (P, den, W retained), and supersaturations are clipped:

```fortran
      a%TV(n)%s_v(1) = min( a%TV(n)%e/a%TV(n)%e_sat(1) - 1.0_PS, 1000.0_PS )
      ...
      a%TV(n)%s_v(2) = min( a%TV(n)%e/a%TV(n)%e_sat(2) - 1.0_PS, 1000.0_RP )
```

(class_AirGroup.F90:293, 298; rest of loop is verbatim-identical to 4b including the rv_sat/s_v sign clipping, `e_sat(2)=min(e_sat(1),e_sat(2))`, `_n` copies, and the `cal_growth_mode_inl_vec` call at 331–336.)

---

## 5. `diag_t` — T recovered from theta_il + qr + qi

`/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_core.F90:12449-12550` (active code; `!tmp`/`!dbg` blocks omitted where noted):

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
  end subroutine diag_t
```

Notes: `qr_0`/`qi_0` are hydrometeor water mass per grid volume (bin mass minus aerosol mass, `imt/imat`, `rmt/rmat`) divided by CGS moist density → mixing ratio. All constants (`p00`, `Racp`, `L_e`, `L_s`, `c_pa`) are the CGS values from mod_amps_const. The two branches invert `theta_il = theta / (1 + (L_e qr + L_s qi)/(c_p max(T,253)))`: linear form when `max(T,253)=253`, quadratic (in T) form otherwise.

### 5b. `cal_til` / `cal_thetail` (parcel-model versions, MKS)

`/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_lib.F90:2217-2257`

```fortran
subroutine cal_til(til,NGRIDS,thetail,pt)
  implicit none
  integer,intent(in) :: NGRIDS
  real(MP_KIND), dimension(NGRIDS)               :: TIL,PT,thetail
  !

  integer :: i,ibr,ibi

  real(MP_KIND) :: pit,theta

  ! this is the easiest parcel model
  do i=1,NGRIDS
     pit=cp*(pt(i)/p0)**(r/cp)
     til(i)=pit/cp*thetail(i)
  enddo
end subroutine cal_til
subroutine cal_thetail(thetail,NGRIDS,qr,qi,pt,T)
  implicit none
  integer,intent(in) :: NGRIDS
  real(MP_KIND), dimension(NGRIDS)               :: TIL,PT, T,qr,qi,thetail
  !

  integer :: i,ibr,ibi
  ! latent heat of condensation in J/kg
  real(MP_KIND), parameter             :: L_e = 2.5e+6
  ! latent heat of sublimation in J/kg
  real(MP_KIND), parameter             :: L_s = 2.8337e+6
  ! gas constant of air
  real(MP_KIND), parameter             :: Ra = 287.0
  real(MP_KIND),parameter :: g=9.8
  real(MP_KIND),parameter :: cp=1004.0

  real(MP_KIND) :: pit,theta

  ! this is the easiest parcel model
  do i=1,NGRIDS
     pit=cp*(pt(i)/p0)**(Ra/cp)
     theta=T(i)*cp/pit
     thetail(i)=theta/(1.0_RP+(L_e*qr(i)+L_s*qi(i))/(cp*max(T(i),253.0_RP)))
  enddo
end subroutine cal_thetail
```

(Note `cal_thetail` uses local MKS constants — J/kg, `cp=1004.0`, `p0` from host module — whereas `diag_t` uses the CGS module constants; the closure formula is the same.)

---

### Extras relevant to porting (same file)

- `get_esi_fesw` / `get_T_fesw` (analytic inversions of the Lowe-Ficke esw formula): `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Thermo_Var.F90:731-759` — constants `a=6.1070, b=17.15, c=273.16, d=38.25`, `T=(d*log(esw/1000.0/a)/b-c)/(log(esw/1000.0/a)/b-1.0_PS)`.
- Lookup table index convention throughout: table index `K` ↔ `T = K+163` K (K=1 → 164 K; estbar spans 164–313 K, esitbar 164–274 K); stored in Pa, accessor `*10.0` → CGS (dyn/cm² = g/s²/cm).
- `make_Thermo_Var3_2` (from RH instead of rv): class_Thermo_Var.F90:309-369; `make_Thermo_Var2` (analytic e_sat): 205-248; `rv_sat = Rdvchiarui*e_sat/max(P-e_sat, e_sat)` everywhere.