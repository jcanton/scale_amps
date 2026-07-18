Below is the complete extraction. All code is quoted verbatim with `file:line` anchors. One important correction to the task's mental model up front, because it changes what "the paths" are:

- The repo does **not** use `iadvv=1`/`iadvv=2`. `sclsedprz_original` branches on `iadvv==12` (Semi-Lagrangian PPM), `iadvv==22` (Semi-Lagrangian mid-point rule), and an `else` (Upstream-linear / forward-Euler — this is the fall-through default). The Euler path lives in the `else` branch.
- `cal_flux_z` is imported into `sclsedprz_original` but is **never called there**. In the `iadvv==12` PPM path the flux is computed **inline** (only `cal_aLRa6_z` is called). `cal_flux_z` is called only by a *different* subroutine, `sclsedaer` (at `mod_amps_utility.F90:5217`). I quote `cal_flux_z` in full anyway since you asked, but note it is not on `sclsedprz_original`'s call path.

---

# AMPS Sedimentation + Terminal Velocity + Repair — Verbatim Extraction

Repo root: `/Users/jcanton/projects/scale_amps`; AMPS source: `contrib/AMPS/`.

---

## 1. `w_terminal_vel` + fall-speed helpers

### 1a. `w_terminal_vel` — `contrib/AMPS/mod_amps_core.F90:3258-3407`

This fills the `vt(2, mxnbin, L)` slots (index 1 = mass-weighted, 2 = concentration-weighted). Cloud droplets (`token==11`) use the stored bulk `%vtm`; rain/ice (`token==1|2`) build linear/cubic sub-bin distribution params via `cal_lincubprms_vec`, then call `cal_wterm_vel_v3_vec`.

```fortran
  subroutine w_terminal_vel(g,ag,vt)
    use class_Group, only: &
       cal_wterm_vel_v3_vec
    use mod_amps_utility, only: cal_lincubprms_vec
    !integer, intent(in)              :: level
    type (group), intent(in)       :: g
    type (airgroup),intent(in)    :: ag
    ! message from reality-check
!tmp    integer,pointer,dimension(:)  :: mes_rc
    !integer,dimension(*)   :: mes_rc
    ! weighted terminal velocity
    ! 1st argument
    !    bin number
    ! 2nd argument
    ! 1: mass weighted one
    ! 2: concentration weighted one
    ! 3: volume weighted one
    real(PS), dimension(2,mxnbin,*)  :: vt
    real(DS), dimension(2,g%n_bin,g%L)    :: vt_dum
    !real(ps), dimension(2)    :: vt_num,vt_num2
    ! parameter of linear distribution
    real(8), dimension(g%n_bin,g%L,4)   :: a2d
    real(8), dimension(g%n_bin,g%L)   :: vec1,vec2
    real(PS), dimension(g%n_bin,g%L,2)   :: binb3d

    ! mass bin boundaries
    !real(DS)                        :: left_bd, right_bd

    integer,dimension(g%n_bin,g%L)   :: error_number
    integer   :: i,j,n

    if( g%token == 11 ) then
      ! --- case of cloud droplets ---
      do n=1,g%l
        if( g%ms(1,n)%con < 1.0e-30_ps .or. g%ms(1,n)%mass(1) < 1.0e-30_ps ) then
          vt(1,1,n)=0.0_PS
          vt(2,1,n)=0.0_PS
        else
          vt(1,1,n)=real(g%ms(1,n)%vtm,PS_KIND)
          vt(2,1,n)=real(g%ms(1,n)%vtm,PS_KIND)
        end if
      enddo
    else if( g%token == 1 .or. g%token == 2 ) then
      ! --- case of rain and ice category ---
!      do in=1,g%n_bin*g%l
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
       do n = 1, g%l
       do i = 1, g%N_BIN
        vt(1,i,n)=0.0_PS
        vt(2,i,n)=0.0_PS

        vec1(i,n)=g%MS(i,n)%con
        vec2(i,n)=g%MS(i,n)%mass(1)
        binb3d(i,n,1)=g%binb(i)
        binb3d(i,n,2)=g%binb(i+1)
      enddo
      enddo

      ! 1. calculate the paramters of linear distribution
      call cal_lincubprms_vec(g%n_bin,g%n_bin,g%L,vec1,vec2,binb3d,&
              a2d,error_number,"w_terminal_vel1")


      ! 2. set up the left and right boundaries of mass bins
!      do in=1,g%n_bin*g%l
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%l
      do i = 1, g%N_BIN
        if(a2d(i,n,4)==-9.99d+100) then
          ! linear distribution
          if( a2d(i,n,1) > 0.0d+0 ) then
            ! --- case of non-negative bin ---
            vec1(i,n)=g%binb(i)   ! left_bd
            vec2(i,n)=g%binb(i+1) ! right_bd
          elseif( a2d(i,n,1) == -1.0d+0 ) then
            ! case of n(x'_1) < 0.0
            vec1(i,n) = a2d(i,n,2)
            vec2(i,n)= g%binb(i+1)
          else if( a2d(i,n,1) == -2.0d+0 ) then
            ! case of n(x'_2) < 0.0
            vec1(i,n) = g%binb(i)
            vec2(i,n) = a2d(i,n,2)
          end if

          ! re-write the parameters as n(m)=a0+a1*m
          a2d(i,n,1)=max(0.0d+0,a2d(i,n,1))-a2d(i,n,2)*a2d(i,n,3)
          a2d(i,n,2)=a2d(i,n,3)
          a2d(i,n,3)=0.0d+0
          a2d(i,n,4)=0.0d+0

        elseif(a2d(i,n,4)>-9.98d+100) then
          ! cubic distribution
          vec1(i,n)=g%binb(i)   ! left_bd
          vec2(i,n)=g%binb(i+1) ! right_bd
        endif
      enddo
      enddo

      ! 3. calculate weighted terminal velocity.
      call cal_wterm_vel_v3_vec(g, ag, 1, &
              a2d, vec1, vec2, error_number, vt_dum)

!      do j=1,2
!        do in=1,g%N_BIN*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
      do j=1,2
         vt(j,i,n)=vt_dum(j,i,n)
      enddo
      enddo
      enddo

      if(debug) then
!        do j=1,2
!          do in=1,g%N_BIN*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN
         do j=1,2
            if(vt(j,i,n)<0.0.or.vt(j,i,n)>1.5e+3_PS)then
              write(*,11) j,i,n,g%token,vt(j,i,n),g%MS(i,n)%vtm
11            format("terminal velocity can be wrong at weight,bin,grid,token,w_v,v",4i5,2es15.6)
            endif
         end do
         end do
         end do
      endif
!      do j=1,2
!        do in=1,g%N_BIN*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
      do j=1,2
          if(vt(j,i,n)<0.0.or.vt(j,i,n)>1.5e+3)then
            vt(j,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),PS_KIND)
          end if
      end do
      end do
      end do

    end if


  end subroutine w_terminal_vel
```

### 1b. Per-bin fall-speed helper `cal_wterm_vel_v3_vec` — `contrib/AMPS/class_Group.F90:8537-9010`

This is the Best-number / Reynolds (`NRe = a*X**b`, Böhm fit) helper. It integrates `Vt = A1*mass^B1` over each bin's sub-distribution with the piecewise `X_b` regime breakpoints, capping at `D_max`/`m_max` for liquid and a density-dependent `v_max` for solid.

```fortran
  subroutine cal_wterm_vel_v3_vec(g,ag, mode, &
       a2d, m1, m2, error_number, WVT)
    use class_Mass_Bin, only: &
       get_aspect_ratio_sol, &
       get_area_ratio_sol, &
       cal_best_number
    use class_Ice_Shape, only: &
       get_total_sfc_area, &
       get_effect_area, &
       get_circum_area, &
       get_perimeter, &
       get_vip
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate the terminal velocity weighted by concentration and mass.
    ! use the approximation of NRe= a*X**b fitted to results by Bohm equation.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    type (Group), intent(in)   :: g
    ! thermo variable object
    type (AirGroup), intent(in)  :: ag

    ! paarameter of linear distribution
    real(8),dimension(g%n_bin,g%L,4)     :: a2d
    real(DS), dimension(2,g%n_bin,*)  :: WVT

    real(8),dimension(g%n_bin,*),intent(in) :: m1,m2
    integer,dimension(g%n_bin,*),intent(in) :: error_number


    real(PS), dimension(4) :: a_rex,b_rex

    real(PS) :: alpha,char_len,OMEGA,P,&
         A_e,A_c,A_e_ic,A_c_ic,q_cs,q,X1,X2
    real(DS) :: A1,B1,D1,Q1,aM,aN
    real(DS),dimension(g%n_bin,g%L) :: aMT,aNT
    real(8),dimension(4)     :: a
    real(8) :: xA,xB,xC,xD,aN_abmx,aM_abmx

    ! parameters to define NRe=a*X**b
    ! for liquid
!!c    real(PS),parameter :: &
!!c         la1=4.072942E-02,lb1=9.953352E-01,&
!!c         la2=7.810803E-02,lb2=7.935192E-01,&
!!c         la3=3.443135E-01,lb3=6.126759E-01,&
!!c         la4=2.372899E+00,lb4=4.645688E-01,&
!!c         LXlim1=2.492336E+01,LXlim2=3.439475E+03,LXlim3=4.490854E+05

    ! does not include diameter > 5.0mm
    real(PS),parameter :: &
         la1=4.049863E-02,lb1=9.939921E-01,&
         la2=7.395895E-02,lb2=8.073323E-01,&
         la3=2.557850E-01,lb3=6.412501E-01,&
         la4=1.033221E+00,lb4=5.235578E-01,&
         LXlim1=2.517541E+01,LXlim2=1.700411E+03,LXlim3=1.356160E+05

    real(DS),parameter :: D_max=0.45,m_max=0.047712938426395
    real(DS) :: v_max

    ! for solid
    ! for ylim   1.000000E+00   2.000000E+01   2.000000E+02
    real(PS),parameter :: &
         sa1=4.742587E-02,sb1=9.910232E-01,&
         sa2=9.447403E-02,sb2=7.722376E-01,&
         sa3=2.549053E-01,sb3=6.299768E-01,&
         sa4=7.880099E-01,sb4=5.236878E-01,&
         SXlim1=2.167592E+01,SXlim2=1.027213E+03,SXlim3=3.934172E+04
    ! for ylim   1.000000E+00   1.000000E+01   2.000000E+02
!!c    real(PS),parameter :: &
!!c         sa1=4.688108E-02,sb1=9.868769E-01,&
!!c         sa2=8.315829E-02,sb2=8.049257E-01,&
!!c         sa3=2.125120E-01,sb3=6.480806E-01,&
!!c         sa4=8.059516E-01,sb4=5.225413E-01,&
!!c         SXlim1=2.221647E+01,SXlim2=3.838900E+02,SXlim3=3.876058E+04

    ! define terminal velocity regimes by masses
    real(PS),dimension(5) :: X_b
    ! masses and diameter devided by terminal velocity regimes.
    real(DS) :: m_sub1,m_sub2,D_1,D_2

    ! +++ mode of shape to calculate the ventilation coefficient
    !     and terminal velocity +++
    ! 1 : spheroid with a_len and c_len
    ! 2 : ice crystal model
    integer,intent(in)   :: mode

    integer :: i_den_gt0p5

    ! level of complexity
    !integer,intent(in)   :: level
    integer  :: i,j,n!,var_Status
    !integer :: em

    real(PS),parameter :: chiarui_1=0.0_PS, chiarui_2=1.0e+30

    if( g%token == 1 ) then
      ! --- in case of liquid phase ---
      a_rex=(/la1,la2,la3,la4/)
      b_rex=(/lb1,lb2,lb3,lb4/)
      X_b=(/chiarui_1,LXlim1,LXlim2,LXlim3,chiarui_2/)


!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        WVT(1,i,n)=0.0_DS
        WVT(2,i,n)=0.0_DS
        aNT(i,n)=0.0_DS
        aMT(i,n)=0.0_DS
      enddo
      enddo

      do j=1,4
!        do in=1,g%N_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN

      no_error_if1: if(error_number(i,n)==0) then

          a(1)=a2d(i,n,1)
          a(2)=a2d(i,n,2)
          a(3)=a2d(i,n,3)
          a(4)=a2d(i,n,4)

          ! calculate Best number for the bin limits
          alpha=g%MS(i,n)%c_len/g%MS(i,n)%a_len

          q = 1.0_PS
!tmp          q = 0.99999_PS

          char_len=2.0_PS*g%MS(i,n)%a_len

          ! ++ calculate Best number for the bin limits +++
          ! coef. from X
          call cal_best_number(D1,Q1,ag%TV(n),q,alpha)

          m_sub1=max(m1(i,n),real(X_b(j)/D1,DS))
          m_sub2=min(m2(i,n),real(X_b(j+1)/D1,DS))

          D_1=(m_sub1/coedpi6)**(1.0/3.0)
          D_2=(m_sub2/coedpi6)**(1.0/3.0)

          X1=D1*m1(i,n)
          X2=D1*m2(i,n)

          xA=m_sub2+m_sub1
          xB=m_sub2-m_sub1
          xC=m_sub2*m_sub1
          xD=m_sub2*m_sub2+m_sub1*m_sub1

          aN=a(1)*xB &
                +0.5d+0*a(2)*xA*xB &
                +a(3)*xB*(xD+xC)/3.0_RP &
                +0.25_RP*a(4)*xD*xA*xB

          aM=0.5d+0*a(1)*xA*xB &
                   +a(2)*xB*(xD+xC)/3.0d+0 &
                   +0.25d+0*a(3)*xD*xA*xB &
                   +0.2d+0*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)



          !
          ! Vt=d_vis/den_a * Nre/d
          !   =d_vis/den_a * (1/d_bar) * (mass_bar/m)**(1/3) * a*X^b
          !   =A1*mass^B1
          !
          A1=ag%TV(n)%d_vis**(1.0_RP-2.0_RP*b_rex(j))*ag%TV(n)%den_a**(b_rex(j)-1.0_RP)*&
             (8.0_PS*gg/PI)**b_rex(j)*(Q1*max(alpha,1.0_RP))**(-b_rex(j))*&
             g%MS(i,n)%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(j)

          B1=b_rex(j)-1.0_DS/3.0_DS

          if(X_b(j+1)<X1.or.X_b(j)>X2) then

          elseif(D_2<D_max) then
            !
            ! in case of all the drops are within D_max
            !
            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(2)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(3)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                     + a(4)/(B1+5.0_DS)*(m_sub2**(B1+5.0_DS)-m_sub1**(B1+5.0_DS)) &
                       ))

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS)) &
                     + a(2)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(4)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                       ))

            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM
          elseif(D_1<D_max) then
            !
            ! in case of D_2 boundary is above D_max
            !
            !   above D_max, the velocity becomes constant.
            !

            v_max=A1*m_max**B1

            xA=m_sub2+m_max
            xB=m_sub2-m_max
            xC=m_sub2*m_max
            xD=m_sub2*m_sub2+m_max*m_max

            aN_abmx=a(1)*xB &
                  +0.5d+0*a(2)*xA*xB &
                  +a(3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a(4)*xD*xA*xB

            aM_abmx=0.5d+0*a(1)*xA*xB &
                     +a(2)*xB*(xD+xC)/3.0d+0 &
                     +0.25d+0*a(3)*xD*xA*xB &
                     +0.2d+0*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)

            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(2)/(B1+3.0_DS)*(m_max**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(3)/(B1+4.0_DS)*(m_max**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                     + a(4)/(B1+5.0_DS)*(m_max**(B1+5.0_DS)-m_sub1**(B1+5.0_DS)) &
                      ) &
                  + &
                  v_max*aM_abmx)

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+1.0_DS)*(m_max**(B1+1.0_DS)-m_sub1**(B1+1.0_DS)) &
                     + a(2)/(B1+2.0_DS)*(m_max**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(3)/(B1+3.0_DS)*(m_max**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(4)/(B1+4.0_DS)*(m_max**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                       ) &
                  + &
                  v_max*aN_abmx)

            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM

          else
            !
            ! in case of all the drops are larger than D_max.
            !
            !   above D_max, the velocity becomes constant.
            !
            v_max=A1*m_max**B1

            xA=m_sub2+m_sub1
            xB=m_sub2-m_sub1
            xC=m_sub2*m_sub1
            xD=m_sub2*m_sub2+m_sub1*m_sub1

            aN_abmx=a(1)*xB &
                  +0.5d+0*a(2)*xA*xB &
                  +a(3)*xB*(xD+xC)/3.0d+0 &
                  +0.25d+0*a(4)*xD*xA*xB

            aM_abmx=0.5d+0*a(1)*xA*xB &
                     +a(2)*xB*(xD+xC)/3.0d+0 &
                     +0.25d+0*a(3)*xD*xA*xB &
                     +0.2d+0*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)

            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                       v_max*aM_abmx)

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                       v_max*aN_abmx)


            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM
          end if
         endif no_error_if1

        end do
        end do
      end do

!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        if(error_number(i,n)==0) then
          WVT(1,i,n)=min(max(0.0_DS,WVT(1,i,n)/aMT(i,n)),1.5e+3_DS)
          WVT(2,i,n)=min(max(0.0_DS,WVT(2,i,n)/aNT(i,n)),1.5e+3_DS)
        elseif(error_number(i,n)<10) then
          WVT(1,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
          WVT(2,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
        endif
      enddo
      enddo

    elseif( g%token == 2 ) then
      ! --- in case of liquid phase ---
      a_rex=(/sa1,sa2,sa3,sa4/)
      b_rex=(/sb1,sb2,sb3,sb4/)
      X_b=(/chiarui_1,SXlim1,SXlim2,SXlim3,chiarui_2/)

!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        WVT(1,i,n)=0.0_DS
        WVT(2,i,n)=0.0_DS
        aNT(i,n)=0.0_DS
        aMT(i,n)=0.0_DS
      enddo
      enddo

      do j=1,4
!        do in=1,g%N_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
         do n = 1, g%L
         do i = 1, g%N_BIN

      no_error_if2: if(error_number(i,n)==0) then

          a(1)=a2d(i,n,1)
          a(2)=a2d(i,n,2)
          a(3)=a2d(i,n,3)
          a(4)=a2d(i,n,4)

          ! +++ get aspect ratio +++
          alpha=get_aspect_ratio_sol(g%MS(i,n),g%IS(i,n))

          ! +++ calculate the total surface area +++
          OMEGA=get_total_sfc_area( g%IS(i,n), g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          ! +++ calculate the area ratio +++
          ! ++++++ calculate the effective cross section ++++++
          A_e = get_effect_area( g%IS(i,n),g%IS(i,n)%sh_type, &
                           g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          A_e_ic = get_effect_area( g%IS(i,n),1, &
                           g%MS(i,n)%a_len, g%MS(i,n)%c_len, mode)

          ! ++++++ calculate the circumscribed cross-sectional area ++++++
          A_c = get_circum_area( g%IS(i,n),g%IS(i,n)%sh_type, &
                           g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          A_c_ic = get_circum_area( g%IS(i,n),1, &
                           g%MS(i,n)%a_len, g%MS(i,n)%c_len, mode)

          ! +++ calculate the characteristic length based on
          !     perimeter and total surface area +++
          P = get_perimeter( g%IS(i,n), g%MS(i,n)%semi_a, g%MS(i,n)%semi_c, mode)

          ! +++ calculate characteristic length +++
          char_len = 2.0*sqrt(A_c/PI)

          ! +++ calculate area ratio +++
          q_cs=(g%MS(i,n)%mass(imc)+g%MS(i,n)%mass(ima)+g%MS(i,n)%mass(imr))&
                     /g%MS(i,n)%con/den_i/&
              get_vip(g%IS(i,n)%is_mod(2),g%IS(i,n)%phi_cs,g%IS(i,n)%semi_aip)

          q = get_area_ratio_sol(g%MS(i,n),g%IS(i,n),mode,A_e,A_c,A_e_ic,A_c_ic,q_cs)


          ! ++ calculate Best number for the bin limits +++
          ! coef. from X
          call cal_best_number(D1,Q1,ag%TV(n),q,alpha)

          X1=D1*m1(i,n)
          X2=D1*m2(i,n)


          if(X_b(j+1)<X1.or.X_b(j)>X2) then

          else

            m_sub1=max(m1(i,n),real(X_b(j)/D1,DS))
            m_sub2=min(m2(i,n),real(X_b(j+1)/D1,DS))

            xA=m_sub2+m_sub1
            xB=m_sub2-m_sub1
            xC=m_sub2*m_sub1
            xD=m_sub2*m_sub2+m_sub1*m_sub1

            aN=a(1)*xB &
                  +0.5_RP*a(2)*xA*xB &
                  +a(3)*xB*(xD+xC)/3.0_RP &
                  +0.25_RP*a(4)*xD*xA*xB

            aM=0.5_RP*a(1)*xA*xB &
                     +a(2)*xB*(xD+xC)/3.0_RP &
                     +0.25_RP*a(3)*xD*xA*xB &
                     +0.2_RP*a(4)*xB*(xA*xA*(xD-xC)+xC*xC)


            !
            ! Vt=d_vis/den_a * Nre/d
            !   =d_vis/den_a * (1/d_bar) * (mass_bar/m)**(1/3) * a*X^b
            !   =A1*mass^B1
            !
            A1=ag%TV(n)%d_vis**(1.0_RP-2.0_RP*b_rex(j))*ag%TV(n)%den_a**(b_rex(j)-1.0)*&
                 (8.0_DS*gg/PI)**b_rex(j)*(Q1*max(alpha,1.0_RP))**(-b_rex(j))*&
                 g%MS(i,n)%mean_mass**(1.0_RP/3.0_RP)/char_len*a_rex(j)

            B1=b_rex(j)-1.0_DS/3.0_DS

            ! 1. mass-weighted terminal velocity
            WVT(1,i,n)=WVT(1,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(2)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(3)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                     + a(4)/(B1+5.0_DS)*(m_sub2**(B1+5.0_DS)-m_sub1**(B1+5.0_DS)) &
                       ))

            ! 2. con-weighted terminal velocity
            WVT(2,i,n)=WVT(2,i,n)+max(0.0_DS,&
                  A1*( a(1)/(B1+1.0_DS)*(m_sub2**(B1+1.0_DS)-m_sub1**(B1+1.0_DS)) &
                     + a(2)/(B1+2.0_DS)*(m_sub2**(B1+2.0_DS)-m_sub1**(B1+2.0_DS)) &
                     + a(3)/(B1+3.0_DS)*(m_sub2**(B1+3.0_DS)-m_sub1**(B1+3.0_DS)) &
                     + a(4)/(B1+4.0_DS)*(m_sub2**(B1+4.0_DS)-m_sub1**(B1+4.0_DS)) &
                       ))

            aNT(i,n)=aNT(i,n)+aN
            aMT(i,n)=aMT(i,n)+aM
          endif

         endif no_error_if2
        end do
        end do
      enddo


!      do in=1,g%N_bin*g%L
!        n=(in-1)/g%N_BIN+1
!        i=in-(n-1)*g%N_BIN
      do n = 1, g%L
      do i = 1, g%N_BIN
        if(error_number(i,n)==0) then
          i_den_gt0p5=0.5*(1.0-sign(1.0_PS,0.5_PS-g%MS(i,n)%den))
          v_max=real(i_den_gt0p5,PS_KIND)*5.0e+3_DS+&
                  (1.0-real(i_den_gt0p5,PS_KIND))*2.0e+3_DS

          WVT(1,i,n)=min(max(0.0_DS,WVT(1,i,n)/aMT(i,n)),v_max)
          WVT(2,i,n)=min(max(0.0_DS,WVT(2,i,n)/aNT(i,n)),v_max)
        elseif(error_number(i,n)<10) then
          WVT(1,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
          WVT(2,i,n)=real(min(max(0.0_ps,g%MS(i,n)%vtm),1.5e+3_ps),8)
        endif
      enddo
      enddo

!dbg!CDIR NOVECTOR
!dbg      do in=1,g%N_bin*g%L
!dbg        n=(in-1)/g%N_BIN+1
!dbg        i=in-(n-1)*g%N_BIN
!dbg        if(error_number(i,n)==0) then
!dbg          if(WVT(1,i,n).gt.1.0e-2.and.WVT(2,i,n)<1.0e-5) then
!dbg            write(*,*) "wvt2 is 0",i,n,g%MS(i,n)%vtm,WVT(1:2,i,n) &
!dbg                ,aMT(i,n),aNT(i,n),a2d(i,n,1:4)
!dbg          endif
!dbg        endif
!dbg      enddo

    end  if


  end subroutine cal_wterm_vel_v3_vec
```

### 1c. `cal_best_number` — `contrib/AMPS/class_Mass_Bin.F90:3279-3310`

```fortran
  subroutine cal_best_number(D1,Q1,th_var,q,alpha)
    ! +++ return Davis number necesarry for calculating 
    !   terminal velocity

    type (Thermo_Var), intent(in)    :: th_var
    ! aspect ratio
    real(PS),intent(in)    :: alpha
    ! area ratio
    real(PS),intent(in) :: q
    real(DS),intent(inout) :: D1,Q1
    real(DS) :: D1_0,D1_1,Q1_0,Q1_1

    integer :: i_q_le1

    i_q_le1=0.5_RP*(1.0_RP+sign(1.0_PS,1.0_PS-q))

    D1_0=8.0_PS*gg*th_var%den/&
          (PI*(th_var%d_vis**2)*max(alpha, 1.0_PS)*(q**0.25))
    Q1_0=q**0.25

    D1_1=8.0_PS*gg*th_var%den/&
          (PI*(th_var%d_vis**2)*max(alpha, 1.0_PS)*q)
    Q1_1=q


    D1=real(i_q_le1,PS_KIND)*D1_0 +&
       (1.0_RP-real(i_q_le1,PS_KIND))*D1_1

    Q1=real(i_q_le1,PS_KIND)*Q1_0 +&
       (1.0_RP-real(i_q_le1,PS_KIND))*Q1_1

  end subroutine cal_best_number
```

---

## 2. `sclsedprz_original` — sedimentation driver

`contrib/AMPS/mod_amps_utility.F90:3668-4980`.

### 2a. Argument list + declarations — lines 3668-3773

```fortran
      subroutine sclsedprz_original(ql,qmlt,qc,qtp &
! <<< 2015/03 T. Hashino changed below
           ,ispray &
           ,n1,n2,n3,n4,npr,nbr,ncr,k1r,k2r,k1m,k2m &
           ,k1br,k2br &
           ,qrp,den,wb,wacc &
           ,spdsfc &
           ,theta,qv,thskin,pgnd &
           ,k1eta,z,zz,dzzm,dzvm,dz1 &         ! z center grid, zz half grid
           ,dtl,i,j,CPPM,CPPME,isnow,iadvv &
!           ,ictpc,ictsw,ictag,ictgr &
!           ,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl &
           ,level,nbhzcl,kmicv,imicv,jmicv &
           ,U &       ! [IN]    u velocity,  cell center
           ,V &       ! [IN]    v velocity,  cell center
           ,TEMP &    ! [IN]    temperature, cell center
           ,dens &    ! [IN]    dry+moist density, cell center
           ,momz &    ! [IN]    z momentum, cell face
           ,CV_WATER &! [IN]    specific heat water
           ,CV_ICE &  ! [IN]    specific heat ice
           ,den_t  &  ! [INOUT] density tendency, cell center
           ,momz_t &  ! [INOUT] z momentum tendency, cell face (added for scale momentum flux output)
           ,rhou_t &  ! [INOUT] rho u velocity tendency, cell center
           ,rhov_t &  ! [INOUT] rho v velocity tendency, cell center
           ,rhoe_t &  ! [INOUT] rho e velocity tendency, cell center
           ,sflx)     ! [INOUT] surface flux
        use scale_prc, only: &
           PRC_abort
        use par_amps
        use maxdims
        use mod_scladv_slppm, only: &
           cal_aLRa6_z, &
           cal_flux_z
        implicit none

      integer,intent(in) :: n1,n2,n3,n4 &
           ,npr,nbr,ncr &
           ,k1br(nbr,*),k2br(nbr,*) &
           ,i,j,k1eta &
           ,isnow,ispray,iadvv,level,nbhzcl
!           ,ictpc,ictsw,ictag,ictgr &
!           ,ircnfl,ipcnfl,iscnfl,iacnfl,igcnfl

      integer,intent(inout) :: k1r,k2r,k1m,k2m
      integer,intent(in) :: imicv(*),jmicv(*),kmicv(*)

      real(RP),intent(in) :: theta(*),qv(*),thskin &
           ,pgnd,den(*),wb(*),wacc(*) &
           ,z(*),zz(*),dzzm(*),dzvm(*) &
           ,spdsfc &
           ,CPPM(11,*),CPPME(11,n2mx,*)
      real(DS),intent(in) :: dtl
      real(RP),intent(inout) :: qrp(npr,nbr,ncr,*)
      real(RP),intent(inout) :: &
           ql(*),qmlt(*),qc(*),qtp(*)

      real(RP), intent(in)     :: dz1

      real(RP), intent(in)    :: momz(n1), U(n1), V(n1), TEMP(n1), dens(n1)
      real(RP), intent(in)    :: CV_WATER, CV_ICE
      real(RP), intent(inout) :: den_t(n1), momz_t(n1), rhou_t(n1), rhov_t(n1), rhoe_t(n1), sflx

! new space
      integer :: k,ka,kb,k1,k2 &
           ,ipr,ibr,icr,niter,iter &
           ,inc,ninc &
           ,isfcset,ibr2,ibr1,kk,nspray
      real(PS) :: aLR(n1,2),a6(n1),dzzz(n1) &
           ,vel(n1,2),vt(n1,2) &
           ,rold(n1),rold_1(n1) &
           ,vel_d(n1),flxsp(n1),fluxf(n1)
      real(PS) :: fluxf_scale(n1), flux(n1), eflx(n1)
      real(PS) :: fluxf_2, rold_2(n1), aLR_2(n1,2),a6_2(n1) ! for scale flux in PPM
      real(PS) :: mass_convert
      integer :: iwv,nwv
      character(len=10) :: isntyp

      real(PS) :: snow &
             ,diab &
             ,tair,tsoil,pignd &
             ,dtnew &
             ,cfl,wt &
             ,dz &
             ,am0,am0p1,am0m1,rm1,rp1,r0,r1,r2,r00,am00,dfmsdr0 &
             ,rinc, s2
      real(PS) :: cpdrd, gcp, coef, denl, rmax, rcld
      integer,dimension(n1) :: ierror1
      integer :: nprx(nprmx)
      integer :: n1mxd,n2mxd,n3mxd
      integer :: istp,ntn
      integer :: jmt_q,jcon_q

! local variables for semi-lagrangian mid-point rule
      real(RP) :: rfdz2(n1), dvel_d(n1), semilag_dist(n1), Z_src
      integer  :: k_dst, k_src(n1)
```

### 2b. Setup / velocity prep common to all schemes — lines 3764-4035

Constants, `nprx`, level→`nwv/jmt_q/jcon_q`, spray flux, then per-`(icr,ibr)` the bin terminal velocity is pulled from `qrp` (last `nwv` parameter slots) and clamped:

```fortran
      cpdrd=cpd/r
      gcp=g/cpd
      coef=4./3.*3.14
      denl=1.e3
      rmax=500.
      rcld=10.

      fluxf_scale = 0.0_RP
      flux = 0.0_RP
      eflx = 0.0_RP
```

Level branch (3824-3841):

```fortran
      if(level.eq.3) then
        nwv=1
        jmt_q=1
        jcon_q=2
      elseif(level.eq.4) then
        nwv=2
        jmt_q=1
        jcon_q=2
      elseif(level.eq.5) then
        nwv=2
        if(isnow.eq.1) then
          jmt_q=imt_q
          jcon_q=icon_q
        else
          jmt_q=rmt_q
          jcon_q=rcon_q
        endif
      end if

      dzzz(k1eta)=dz1
      do k=k1eta+1,n1
        dzzz(k)=1.0/dzzm(k)
      end do
```

Terminal velocity for the current bin (3947-3956):

```fortran
          do iwv=1,nwv
            do k=k1,k2
              vt(k,iwv)=qrp(npr-nwv+iwv,ibr,icr,k)
              if(qrp(1,ibr,icr,k).ne.0.0.and.vt(k,iwv).le.0.) then
                vt(k,iwv)=sign(min(abs(vt(k,iwv)*wacc(k)),50.0_RP),wacc(k))
              else
                vt(k,iwv)=0.
              endif
            enddo
          enddo
```

### 2c. CFL substep count logic — lines 3967-3990

`niter = 1 + cfl/0.95`, with `ka/kb` widened by `niter`:

```fortran
          cfl=0.
          do k=max(k1,k1eta+1),k2
            cfl=max(real(cfl,DS), &
                max(abs(vt(k,1)),abs(vt(k-1,1)) &
                   ,abs(vt(k,2)),abs(vt(k-1,2)))*dtl &
             /(zz(k)-zz(k-1)))
          enddo
!
!----------------------------------------------------------------------------------------
!     Find how many iterations (small timesteps)  necessary for stability
!----------------------------------------------------------------------------------------
!
          niter=1.+cfl/0.95
...
          ka=max(k1eta,k1-niter+1)
          kb=min(n1-1,k2+niter-1)
```

Face-velocity fill (`velfill` + averaging to "w" points) — lines 3995-4035 (`vel(k,iwv)=0.5*(vt(k,iwv)+vt(k+1,iwv))`, with one-sided ends).

The `iwv` selection per parameter `ipr` (mass- vs con-weighted) — lines 4063-4101; the mass-weighted vs con-weighted `vel_d` is chosen here.

### 2d. `advscheme_if` opener — line 4103

```fortran
            advscheme_if: if(iadvv.eq.12) then
```

### 2e. Euler (Upstream-linear) path — the DEFAULT `else` branch — lines 4611-4827

This is the forward-upstream scheme. It has a multi-iteration sub-branch (`niter>1`, lines 4621-4725) and a single-iteration sub-branch (`else`, lines 4726-4826). Full verbatim:

```fortran
            else
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Upstream linear
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Multiple iterations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
              dtnew=dtl/float(niter)
              if(niter.gt.1) then
                !write(*,*) "niter is more than 1",niter,dtnew
                do iter=1,niter

!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                 Forward - Upstream precip setteling over niter iterations
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!
                  if (isnow == 1) then
                    if (ipr <= 8) then
                      mass_convert = 1.0D0
                    else
                      mass_convert = 0.001D0
                    endif
                  else
                    if (ipr <= 3) then
                      mass_convert = 1.0D0
                    else
                      mass_convert = 0.001D0
                    endif
                  endif
                  do k=ka-1,kb
                    if(vel_d(k).gt.0.0_RP) then
                      fluxf(k)=vel_d(k)*rold(k)*den(k)*mass_convert

                      if (isnow == 1 .and. (ipr==imt_q)) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(imw_q,ibr,icr,k) - qrp(imat_q,ibr,icr,k))*den(k)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_ICE &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow == 1 .and. (ipr==imw_q)) then
                        fluxf_2 = vel_d(k)*rold(k)*den(k)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow /= 1 .and. ipr==rmt_q) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(rmat_q,ibr,icr,k))*den(k)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      endif
                      LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                      call PRC_abort
                    else
                      fluxf(k)=vel_d(k)*rold(k+1)*den(k+1)*mass_convert

                      if (isnow == 1 .and. (ipr==imt_q)) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(imw_q,ibr,icr,k+1) - qrp(imat_q,ibr,icr,k+1))*den(k+1)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_ICE &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow == 1 .and. ipr==imw_q) then
                        fluxf_2 = vel_d(k)*rold(k+1)*den(k+1)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      elseif (isnow /= 1 .and. ipr==rmt_q) then
                        fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(rmat_q,ibr,icr,k+1))*den(k+1)/float(niter)
                        fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                        eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                                + fluxf_2/dzvm(k)*g
                      endif
                    endif
                  enddo
!
!
!
!
!------------------------------------------------------------------------
!                 Precip Rate and accumulation Rate for iteration case
!------------------------------------------------------------------------
!
                  if(ka.eq.k1eta) then
                    if(ipr.eq.1) then
                        ! scale surface flux
                        sflx = sflx + fluxf(k1eta-1)/float(niter)
                    endif
!
!     surface spray flux
!
                    if(ipr.eq.1.and.icr.eq.1.and.ispray.eq.1) then
                      fluxf(k1eta-1)=fluxf(k1eta-1)+flxsp(ibr)
                    endif
                  endif
!
!     predict settling now!
!
                  do k=ka,kb
                    if(k.eq.k1eta) then
                      rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                        /((z(k)-z(k-1))*den(k)*mass_convert)
                      !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      !  /((z(k)-z(k-1))*den(k))
                    else
                      rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                        /((zz(k)-zz(k-1))*den(k)*mass_convert)
                      !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtnew &
                      !  /((zz(k)-zz(k-1))*den(k))
                    endif
                  enddo
                enddo
              else
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Single Iteration
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!                 Forward - Upstream vertical advection just 1 iteration
!------------------------------------------------------------------------
!------------------------------------------------------------------------
                if (isnow == 1) then
                  if (ipr <= 8) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                else
                  if (ipr <= 3) then
                    mass_convert = 1.0D0
                  else
                    mass_convert = 0.001D0
                  endif
                endif
                do k=ka-1,kb
                  if(vel_d(k).gt.0.0_RP) then
                    fluxf(k)=vel_d(k)*rold(k)*den(k)*mass_convert

                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(imw_q,ibr,icr,k) - qrp(imat_q,ibr,icr,k))*den(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_ICE &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow == 1 .and. (ipr==imw_q)) then
                      fluxf_2 = vel_d(k)*rold(k)*den(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k) - qrp(rmat_q,ibr,icr,k))*den(k)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    endif
                    LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                    call PRC_abort
                  else
                    fluxf(k)=vel_d(k)*rold(k+1)*den(k+1)*mass_convert

                    if (isnow == 1 .and. (ipr==imt_q)) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(imw_q,ibr,icr,k+1) - qrp(imat_q,ibr,icr,k+1))*den(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_ICE &
                              + fluxf_2/dzvm(k)*g           ! potential energy
                    elseif (isnow == 1 .and. (ipr==imw_q)) then
                      fluxf_2 = vel_d(k)*rold(k+1)*den(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g               ! potential energy
                    elseif (isnow /= 1 .and. ipr==rmt_q) then
                      fluxf_2 = vel_d(k)*max(0.0_RP,rold(k+1) - qrp(rmat_q,ibr,icr,k+1))*den(k+1)
                      fluxf_scale(k) = fluxf_scale(k) + fluxf_2
                      eflx(k) = eflx(k) + fluxf_2*TEMP(k+1)*CV_WATER &
                              + fluxf_2/dzvm(k)*g           ! potential energy
                    endif
                  endif
                enddo
                if(ka.eq.k1eta) then
!
!------------------------------------------------------------------------
!                 Precip Rate and accumulation Rate for no iteration case
!------------------------------------------------------------------------
!
                  if(ipr.eq.1) then
                      ! scale surface flux
                      sflx = sflx + fluxf(k1eta-1)
                  endif
!
!     surface spray flux
!
                  if(ipr.eq.1.and.icr.eq.1.and.ispray.eq.1) then
                    fluxf(k1eta-1)=fluxf(k1eta-1)+flxsp(ibr)
                  endif
                endif
!
!     predict settling now!
!
                do k=ka,kb
                  if(k.eq.k1eta) then
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                      /((z(k)-z(k-1))*den(k)*mass_convert)
                    !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                    !  /((z(k)-z(k-1))*den(k))
                  else
                    rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                      /((zz(k)-zz(k-1))*den(k)*mass_convert)
                    !rold(k)=rold(k)-(fluxf(k)-fluxf(k-1))*dtl &
                    !  /((zz(k)-zz(k-1))*den(k))
                  endif
                enddo
              endif
            endif advscheme_if
```

### 2f. PPM path (`iadvv==12`) structure + inline flux — lines 4103-4302

The PPM path calls `cal_aLRa6_z` (NOT `cal_flux_z`) to get left/right interface values (`aLR`) and the parabola coefficient (`a6`), then computes the fractional-cell PPM flux inline. A second "scale" distribution (`rold_2 = total − aerosol − meltwater`) is built and integrated in parallel to fill `fluxf_scale`/`eflx`. Structure:

```fortran
            advscheme_if: if(iadvv.eq.12) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            Semi-Lagrangian PPM
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              dtnew=dtl/float(niter)
...
              do iter=1,niter

!     semi-Lagrangian PPM method (modo Ong Chia Rui)
!     1. calculation of points at left and right interfaces and coefficient.
                call cal_aLRa6_z(i,j,max(k1eta,ka-1),min(n1,kb+1) &
                       ,k1eta,n1,n1mxd,n2mxd,n3mxd &
                       ,n1,istp &
                       ,CPPM,CPPME &
                       ,rold(1),aLR,a6)

                if (isnow == 1 .and. (ipr==imt_q)) then
                  do k=1,n1
                    rold_2(k) = max(0.0_RP,rold(k) - qrp(imw_q,ibr,icr,k) - qrp(imat_q,ibr,icr,k))
                  enddo
                  call cal_aLRa6_z(i,j,max(k1eta,ka-1),min(n1,kb+1) &
                         ,k1eta,n1,n1mxd,n2mxd,n3mxd &
                         ,n1,istp &
                         ,CPPM,CPPME &
                         ,rold_2(1),aLR_2,a6_2)
                elseif (isnow == 1 .and. (ipr==imw_q)) then
                  rold_2 = rold
                  aLR_2 = aLR
                  a6_2 = a6
                elseif (isnow /= 1 .and. (ipr==rmt_q)) then
                  do k=1,n1
                    rold_2(k) = max(0.0_RP,rold(k) - qrp(rmat_q,ibr,icr,k))
                  enddo
                  call cal_aLRa6_z(i,j,max(k1eta,ka-1),min(n1,kb+1) &
                         ,k1eta,n1,n1mxd,n2mxd,n3mxd &
                         ,n1,istp &
                         ,CPPM,CPPME &
                         ,rold_2(1),aLR_2,a6_2)
                endif
```

The inline flux for the (positive-velocity, aborted) and (negative-velocity, used) cases — lines 4172-4211:

```fortran
                  if(vel_d(k).gt.0.0_RP) then
                    s2=(zz(k)-(zz(k)-vel_d(k)*dtnew))/dzzz(k)
                    fluxf(k)= &
                      den(k)*(aLR(k,2)-0.5*s2*(aLR(k,2)-aLR(k,1) &
                         -(1.0-2.0*s2/3.0)*a6(k))) &
                     *vel_d(k)*mass_convert
                    ...
                    LOG_ERROR("sclsedprz_original",*) "TVEL>0"
                    call PRC_abort
                  else
                    s2=((zz(k)-vel_d(k)*dtnew)-zz(k))/dzzz(k+1)
                    fluxf(k)= &
                      den(k+1)*(aLR(k+1,1)+0.5*s2*(aLR(k+1,2)-aLR(k+1,1) &
                          +(1.0-2.0*s2/3.0)*a6(k+1))) &
                     *vel_d(k)*mass_convert
                    ...
                  endif
```

(The `iadvv==22` Semi-Lagrangian mid-point path is lines 4304-4609; it uses a Taylor-expanded backward-trajectory `semilag_dist` and integer-cell sweeping to build `fluxf`. Same `fluxf_scale`/`eflx` accumulation structure.)

### 2g. Per-parameter deposition back into `ql/qmlt/qc/qtp` and write-back to `qrp` — lines 4835-4930

After each scheme updates `rold`, the difference vs the old `qrp` is deposited into the appropriate bulk arrays (level 3/4 into `ql`+`qtp`; level≥5 splits total-minus-aerosol into `qc`/`ql`, meltwater into `qmlt`, with an `ierror1`/`any()` negativity guard on `qtp`), then `qrp(ipr,ibr,icr,k)=rold(k)`.

### 2h. **M2 replication target — the momentum + energy tendency accumulation** — lines 4932-4980

After all `icr/ibr/ipr` loops, the accumulated `fluxf_scale` (mass flux of the "condensate minus aerosol/meltwater" part) and `eflx` (internal + potential energy flux) drive the SCALE tendencies:

```fortran
!
!  momentum flux for scale
!
      ! z center grid, zz half grid
      do k = k1eta, n1 ! k1eta = 2
         if(k .eq. k1eta) then
            den_t(k) = den_t(k)-(fluxf_scale(k)-fluxf_scale(k-1))/(z(k)-z(k-1))
         else
            den_t(k) = den_t(k)-(fluxf_scale(k)-fluxf_scale(k-1))/(zz(k)-zz(k-1))
         endif
      enddo
      flux(n1) = 0.0_RP
      ! AMPS z-grid index: 1~n1, 1 is underground, n1=KE-KS+2. SCALE z-grid index: KE-1.
      ! please make reference to scale_atmos_phy_mp_common.f90/ATMOS_PHY_MP_precipitation_mom
      do k = 2, n1-1
        flux(k) = (fluxf_scale(k) + fluxf_scale(k-1))*momz(k)/(dens(k) + dens(k+1))
      enddo
      !flux(n1-1) = 0.0_RP

      ! z momentum
      do k = 2, n1-1
        momz_t(k) = momz_t(k) - (flux(k+1) - flux(k))*dzvm(k)
      enddo
      momz_t(n1) = 0.0_RP

      ! x momentum
      do k = 1, n1-1
        flux(k) = fluxf_scale(k)*U(k+1)
      enddo
      flux(n1) = 0.0_RP
      do k = 2, n1
        rhou_t(k) = rhou_t(k) - (flux(k) - flux(k-1))*dzzm(k)
      enddo

      ! y momentum
      do k = 1, n1-1
        flux(k) = fluxf_scale(k)*V(k+1)
      enddo
      flux(n1) = 0.0_RP
      do k = 2, n1
        rhov_t(k) = rhov_t(k) - (flux(k) - flux(k-1))*dzzm(k)
      enddo

      ! internal energy
      do k = 2, n1
        rhoe_t(k) = rhoe_t(k) - (eflx(k) - eflx(k-1))*dzzm(k)
      enddo
      return
      end subroutine sclsedprz_original
```

Note the grid staggering: `den_t`/`rhou_t`/`rhov_t`/`rhoe_t` are divergence of a face flux onto centers (`dzzm`), `momz_t` uses `dzvm`, the surface flux `sflx` is accumulated inside the loops (`sflx = sflx + fluxf(k1eta-1)[/niter]` for `ipr==1`), and `SFLX` is the surface value `fluxf(k1eta-1)`.

### 2i. PPM helper `cal_flux_z` (imported but not called by this driver) — `contrib/AMPS/mod_scladv_slppm.F90:423-593`

```fortran
subroutine cal_flux_z(k,k1eta,n1,n1mx,lbin,zz,dz,dt,den,a,aLR,a6,zarr,w,flux,from)
  use scale_prc, only: &
     PRC_abort
  implicit none
!!c  integer,intent(in) :: k,k1eta,n1,n1mx
!!c  real,dimension(n1) :: a,a6,dz,den
!!c  real,dimension(n1mx) :: z
!!c  real,dimension(n1,2) :: aLR
!!c  real,intent(in)    :: w,zarr
!!c  real,intent(inout) :: flux
  integer,intent(in) :: k,k1eta,n1,n1mx,lbin
  real(PS) :: a(lbin),a6(n1),dz(n1)
  real(RP) :: den(lbin),zz(n1mx),zarr
  real(PS),dimension(n1,2) :: aLR
  real(PS),intent(in)    :: w
  real(PS),intent(inout) :: flux
  character*(*) :: from

  real(PS) :: dt,s1,s2,yisum,zdep
  integer :: kdep,kck
  integer :: ka,n
  
  
  ! find departure point
  ! search for departure point
  kck=0
  zdep=zarr-w*dt

  if(zdep<zz(k1eta)-dz(k1eta)) then
     zdep=zz(k1eta)-dz(k1eta)
     kck=1
  else if(zdep>zz(n1-1)) then
     ! case of periodic boundary
     zdep=zz(n1-1)
     kck=1
  end if
 
  ! find the interface departure point
!!c  kdep=int((zdep-zz(1))/dz(2))+2
  if(w<0.0) then
    do n=k+1,n1
       if(zz(n)>zdep.and.zz(n-1)<=zdep) then
          kdep=n
          exit
       end if
    end do
  else
    do n=k,k1eta,-1
!!c       if(zz(n)>zdep.and.zz(n-1)<=zdep) then
       if(zz(n)>=zdep.and.zz(n-1)<zdep) then
          kdep=n
          exit
       end if
    end do
 end if
  
  
!!c  ka=kdep-1
  ka=kdep
  
  flux=0.0
  if(w>0.0) then
     ! ka is the cell number for which the interpolation function is used.
     ka=kdep
     ! calculate integer flux
     yisum=0.0
...
     do n=ka+1,k
        flux=flux+a(n)*den(n)*dz(n)
!!c        flux=flux+a(n)*dz(n)
        yisum=yisum+dz(n)
     end do

     ! calculate fractional flux
     s1=zz(kdep)-zdep
     s2=s1/dz(ka)
     if(yisum+s1>0.0_PS) then
        flux=(flux+&
!!c             (aLR(ka,2)-0.5*s2*(aLR(ka,2)-aLR(ka,1)-(1.0-2.0*s2/3.0)*a6(ka)))*s1&
             den(ka)*(aLR(ka,2)-0.5_PS*s2*(aLR(ka,2)-aLR(ka,1)-(1.0_PS-2.0_PS*s2/3.0_PS)*a6(ka)))*s1&
             )/(yisum+s1)
     end if

     if ( debug ) then
        if(flux>1.0e+10_PS.or.flux<-1.0e+10_PS) then
           write(*,*) from
           write(*,'("flux large in slppm 1",2I5,20ES15.6)') kdep,ka,flux,w,den(ka),aLR(ka,1),aLR(ka,2),a6(ka),s1,yisum,dz(ka),s2,zdep
        end if
     end if

  elseif(w<0.0_PS) then
     ka=kdep
     ! calculate integer flux
     yisum=0.0_PS
...
     do n=k+1,ka-1
        flux=flux+a(n)*den(n)*dz(n)
!!c        flux=flux+a(n)*dz(n)
        yisum=yisum+dz(n)
     end do
     ! calculate fractional flux
     if(kdep-1>=1) then
        s1=zdep-zz(kdep-1)
     else 
        LOG_ERROR("cal_flux_z",*) "something wrong in vdy_scladv_slppm flux_z",kdep
        call PRC_abort
     end if
     s2=s1/dz(ka)
     if(yisum+s1>0.0) then
        flux=(flux+&
!!c             (aLR(ka,1)+0.5*s2*(aLR(ka,2)-aLR(ka,1)+(1.0-2.0*s2/3.0)*a6(ka)))*s1&
             den(ka)*(aLR(ka,1)+0.5*s2*(aLR(ka,2)-aLR(ka,1)+(1.0-2.0*s2/3.0)*a6(ka)))*s1&
             )/(yisum+s1)
     end if
     if ( debug ) then
        if(flux>1.0e+10.or.flux<-1.0e+10) then
           write(*,*) from
           write(*,'("flux large in slppm 2",2I5,20ES15.6)') kdep,ka,flux,w,den(ka),aLR(ka,1),aLR(ka,2),a6(ka),s1,yisum,dz(ka),s2,zdep
        end if
     end if
  end if
  flux=flux*w

  if ( debug ) then
     if(flux>1.0e+10.or.flux<-1.0e+10) then
        write(*,*) from
        write(*,'("flux large in slppm 3",2I5,20ES15.6)') kdep,ka,flux,w,den(ka),aLR(ka,1),aLR(ka,2),a6(ka),s1,yisum,dz(ka),s2,zdep
     end if
  end if
  
end subroutine cal_flux_z
```

### 2j. `cal_aLRa6_z` (PPM interface + parabola coefficients) — `contrib/AMPS/mod_scladv_slppm.F90:272-421`

```fortran
subroutine cal_aLRa6_z(i,j,ns,ne,k1eta,n1,n1mx,n2mx,n3mx,lbin,istp,CPPM,CPPME,a,aLR,a6)
  implicit none
  ! ns and ne is the starting and ending grid cell for aLR and a6
  integer,intent(in) :: i,j,ns,ne
  integer,intent(in) :: k1eta,n1,n1mx,n2mx,n3mx,lbin,istp
  real(PS) :: a(lbin),a6(n1)
  real(PS) :: aLR(n1,2)
  real(RP) :: CPPM(11,n1mx),CPPME(11)
  real(PS),dimension(n1) :: da,da2
  real(PS) :: ita,ita_w, aLd, aRd
  real(PS),parameter  :: ita_1=20.0,ita_2=0.05,eps=0.01
  integer :: n,m,k,l
  
  
!!c  open(1,file='check_cppm_eta.dat') 
!!c    write(1,'(11ES15.6)') (CPPME(m),m=1,11)
!!c  close(1)
!!c  stop

  ! calculate the slopes
  if(ns<=k1eta+1) then
     da(k1eta)=CPPME(1)*a(k1eta-1)+CPPME(2)*a(k1eta)+CPPME(3)*a(k1eta+1)
     if((a(k1eta+1)-a(k1eta))*(a(k1eta)-a(k1eta-1))>0.0) then
        da(k1eta)=min(abs(da(k1eta)),2.0*abs(a(k1eta)-a(k1eta-1)),2.0*abs(a(k1eta+1)-a(k1eta)))*sign(1.0_PS,da(k1eta))
     else
        da(k1eta)=0.0
     end if
     da2(k1eta)=CPPME(4)*a(k1eta-1)+CPPME(5)*a(k1eta)+CPPME(6)*a(k1eta+1)
  end if

  do k=max(ns-1,k1eta+1),min(n1-1,ne+1)
     da(k)=CPPM(1,k)*a(k-1)+CPPM(2,k)*a(k)+CPPM(3,k)*a(k+1)
     
     if((a(k+1)-a(k))*(a(k)-a(k-1))>0.0) then
        da(k)=min(abs(da(k)),2.0*abs(a(k)-a(k-1)),2.0*abs(a(k+1)-a(k)))*sign(1.0_PS,da(k))
     else
        da(k)=0.0
     end if
     
     da2(k)=CPPM(4,k)*a(k-1)+CPPM(5,k)*a(k)+CPPM(6,k)*a(k+1)
  end do

  if(ns==k1eta) then
     da(k1eta-1)=0.0
     da2(k1eta-1)=0.0
  endif
  if(ne>=n1-1) then
     da(n1)=0.0
     da2(n1)=0.0
  endif 
  
  ! interpolation at the left and right interfaces of cell i
  if(ns<=k1eta+1) then
     aLR(k1eta,2)=CPPME(7)*a(k1eta)+CPPME(8)*a(k1eta+1)+CPPME(9)*da(k1eta)+CPPME(10)*da(k1eta+1)
     aLR(k1eta+1,1)=aLR(k1eta,2)
  end if
  do k=max(ns-1,k1eta+1),min(n1-2,ne)
     aLR(k,2)=CPPM(7,k)*a(k)+CPPM(8,k)*a(k+1)+CPPM(9,k)*da(k)+CPPM(10,k)*da(k+1)
     aLR(k+1,1)=aLR(k,2)
  end do
  if(ns==k1eta) then
     aLR(k1eta,1)=aLR(k1eta,2)
!!c     aLR(k1eta,1)=0.0
  endif
  if(ne>=n1-1) then
     aLR(n1-1,2)=a(n1)
  endif

  if(istp==1) then
     if(ns==k1eta) then
        if(-da2(k1eta+1)*da2(k1eta-1)>0.0.or.abs(a(k1eta+1)-a(k1eta-1))-eps*min(abs(a(k1eta+1)),abs(a(k1eta-1)))>0.0) then
           ita_w=CPPME(11)*(da2(k1eta+1)-da2(k1eta-1))/(a(k1eta+1)-a(k1eta-1))
        else
           ita_w=0.0
        end if
        ita=max(0.0_PS,min(ita_1*(ita_w-ita_2),1.0_PS))
        
        aLd=a(k1eta-1)+0.5*da(k1eta-1)
        aRd=a(k1eta+1)-0.5*da(k1eta+1)
        
        aLR(k1eta,1)=aLR(k1eta,1)*(1.0-ita)+aLd*ita
        aLR(k1eta,2)=aLR(k1eta,2)*(1.0-ita)+aRd*ita
     end if
     
     ! treatment of discontinuity in the cell
     do k=max(ns,k1eta+1),min(n1-1,ne)
        if(-da2(k+1)*da2(k-1)>0.0.or.abs(a(k+1)-a(k-1))-eps*min(abs(a(k+1)),abs(a(k-1)))>0.0) then
           ita_w=CPPM(11,k)*(da2(k+1)-da2(k-1))/(a(k+1)-a(k-1))

!!c           ita_w=-((da2(i+1)-da2(i-1))/(0.5*dx(i+1)+dx(i)+0.5*dx(i-1)))*&
!!c                (( (0.5*dx(i)+0.5*dx(i-1))**3.0+(0.5*dx(i+1)+0.5*dx(i))**3.0)/(a(i+1)-a(i-1)))
        else
           ita_w=0.0
        end if
        ita=max(0.0_PS,min(ita_1*(ita_w-ita_2),1.0_PS))
        
        aLd=a(k-1)+0.5*da(k-1)
        aRd=a(k+1)-0.5*da(k+1)
        
        aLR(k,1)=aLR(k,1)*(1.0-ita)+aLd*ita
        aLR(k,2)=aLR(k,2)*(1.0-ita)+aRd*ita
     end do
  end if
  
  ! monotonicity algorithm
  do k=ns,min(n1-1,ne)


     if((aLR(k,2)-a(k))*(a(k)-aLR(k,1))<0.0) then
       aLR(k,1)=a(k)
       aLR(k,2)=a(k)          
     else 
...
       if((aLR(k,2)-aLR(k,1))*(a(k)-0.5*(aLR(k,1)+aLR(k,2)))>&
          (aLR(k,2)-aLR(k,1))*(aLR(k,2)-aLR(k,1))/6.0) then


         aLR(k,1)=3.0*a(k)-2.0*aLR(k,2)
       else if((aLR(k,2)-aLR(k,1))*(a(k)-0.5*(aLR(k,1)+aLR(k,2)))<&
          -(aLR(k,2)-aLR(k,1))*(aLR(k,2)-aLR(k,1))/6.0) then
         aLR(k,2)=3.0*a(k)-2.0*aLR(k,1)
       endif
     end if
     
  end do
! this is for just upstream scheme
!  if(ns==k1eta) then
!     aLR(k1eta,1)=a(k1eta)
!     aLR(k1eta,2)=a(k1eta)
!  endif
  if(ne==n1) then
     aLR(n1,1)=a(n1)
     aLR(n1,2)=a(n1)
  endif 
  
  ! calculation of coefficient
  do k=ns,ne
     a6(k)=6.0*(a(k)-0.5*(aLR(k,1)+aLR(k,2)))
  end do
  
end subroutine cal_aLRa6_z
```

### 2k. `cal_ppmcoef` (grid-dependent PPM coefficient table) — `contrib/AMPS/mod_scladv_slppm.F90:14-97`

```fortran
subroutine cal_ppmcoef(ns,n2,n2max,xx,CPPM)
  implicit none
  integer,intent(in) :: ns,n2,n2max
  real(MP_KIND), intent(in), dimension(n2) :: xx ! n2max -> n2
  real(MP_KIND), intent(out), dimension(11,n2max) :: CPPM ! how n2max -> n2

  real(MP_KIND),dimension(n2) :: dx
  integer :: i,j
  
  ! calculate the difference of grid interfaces
  do i=2,n2-1
     dx(i)=xx(i)-xx(i-1)
  end do
  dx(1)=dx(2)
  dx(n2)=dx(n2-1)

  ! coefficients for calculation of slope
  do i=2,n2-1
     !! for da1, these are unitless
     ! for a(i-1)
     CPPM(1,i)=-(dx(i)+2.0*dx(i+1))/(dx(i-1)+dx(i))
     ! for a(i)
     CPPM(2,i)=-(2.0*dx(i-1)+dx(i))/(dx(i+1)+dx(i))+(dx(i)+2.0*dx(i+1))/(dx(i-1)+dx(i))
     ! for a(i+1)
     CPPM(3,i)=(2.0*dx(i-1)+dx(i))/(dx(i+1)+dx(i))
     do j=1,3
        CPPM(j,i)=CPPM(j,i)*dx(i)/(dx(i-1)+dx(i)+dx(i+1))
     end do
     
     !! for da2, these have unit of L^(-2)
     ! for a(i-1)
     CPPM(4,i)=1.0/(dx(i)+dx(i-1))/(dx(i-1)+dx(i)+dx(i+1))
     ! for a(i)
     CPPM(5,i)=(-1.0/(dx(i+1)+dx(i))-1.0/(dx(i)+dx(i-1)))/(dx(i-1)+dx(i)+dx(i+1))
     ! for a(i+1)
     CPPM(6,i)=1.0/(dx(i+1)+dx(i))/(dx(i-1)+dx(i)+dx(i+1))
     
  end do
  
  
  ! coefficients for interpolation of interface values
  do i=ns,n2-2
     ! for a(i), dimensionless
     CPPM(7,i)=1.0-dx(i)/(dx(i)+dx(i+1))+&
          ((dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))-(dx(i+2)+dx(i+1))/(2.0*dx(i+1)+dx(i)))*&
          (-1.0)*2.0*dx(i+1)*dx(i)/(dx(i)+dx(i+1))/&
          (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))
     
     ! for a(i+1), dimensionless
     CPPM(8,i)=dx(i)/(dx(i)+dx(i+1))+&
          ( (dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))-(dx(i+2)+dx(i+1))/(2.0*dx(i+1)+dx(i)) )*&
          (1.0)*2.0*dx(i+1)*dx(i)/(dx(i)+dx(i+1))/&
          (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))
     
     
     ! for da(i), dimensionless
     CPPM(9,i)=dx(i+1)*(dx(i+1)+dx(i+2))/(dx(i)+2.0*dx(i+1))/&
            (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))            
     ! for da(i+1), dimensionless
     CPPM(10,i)=(-1.0)*dx(i)*(dx(i-1)+dx(i))/(2.0*dx(i)+dx(i+1))/&
          (dx(i-1)+dx(i)+dx(i+1)+dx(i+2))
     
...

     ! dimension is L^2
     CPPM(11,i)=-1.0/(0.5*dx(i+1)+dx(i)+0.5*dx(i-1))*&
                ( (0.5*dx(i)+0.5*dx(i-1))**3.0+(0.5*dx(i+1)+0.5*dx(i))**3.0)

  end do
  CPPM(7,n2-1)=CPPM(7,n2-2)
  CPPM(8,n2-1)=CPPM(8,n2-2)
  CPPM(9,n2-1)=CPPM(9,n2-2)
  CPPM(10,n2-1)=CPPM(10,n2-2)
  CPPM(11,n2-1)=CPPM(11,n2-2)
  
  
end subroutine cal_ppmcoef
```

---

## 3. `shift_bin_vec` — post-sedimentation mean-mass bin reshift

`contrib/AMPS/mod_amps_utility.F90:6584-7720`. Full verbatim. It rebuilds each spectrum after sedimentation has shifted mean masses (`mmass_p = total/con`) relative to the fixed bin boundaries: it computes shifted boundaries `binb3d = binbr ± (mmass_p − mmass)`, fits sub-bin linear/cubic params via `cal_lincubprms_vec`, then redistributes back to the original bins via `cal_transbin_vec`. `iphase==1` is liquid, `iphase==2` is ice.

```fortran
      subroutine shift_bin_vec(n1,npr,nbr,ncr,k1r,k2r,binbr &
          ,qrpv,mmass,den,iphase &
          ,kmic,imic,jmic)
!     this subroutine has to be consistent with how to shift bins in
!     class_group
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      use maxdims
      use com_amps, only: fid_alog
!      use com_amps
      use par_amps
      implicit none
      integer :: n1,npr,nbr,ncr,k1r,k2r,iphase
      real(RP) :: qrpv(npr,nbr,ncr,*),mmass(nbr,ncr,*),den(*)
      real(PS),dimension(*),intent(in) :: binbr
      integer :: kmic(*),imic(*),jmic(*)
      integer :: naxis,nmass,ncon
      real(RP) :: den_i
      parameter(den_i= 0.91668)


      real(8),dimension(mxnbin,n1)               :: new_N
      real(8),dimension(mxnbin,n1,mxnmasscomp+1) :: new_M
      real(8),dimension(mxnbin,n1,mxnnonmc)      :: new_Q

      real(RP),dimension(mxnbin,n1,mxnnonmc) :: Qp


      ! total concentration in the shifted bin
      real(8), dimension(mxnbin+1,n1)           :: Np
      ! total mass in the shifted bin
      real(8), dimension(mxnbin+1,n1)           :: Mp

      ! averaged mass tendency after the time step for the bin
      real(8), dimension(mxnbin,n1)             :: new_mtend

      ! parameter of distribution in each bin
      real(8), dimension(mxnbin+1,n1,4)           :: a2d

!     shifted bounds of a bin, number of parameters for subdistribution
      real(RP),dimension(mxnbin+1,n1,2) :: binb3d
      ! mass tendency of shifted bin
      real(PS), dimension(mxnbin+1,n1)                     :: mtend
      real(RP),dimension(mxnbin+1,n1,mxnmasscomp) :: ratio_Mp

      ! axis ratio of an ice crystal in a shifted bin
      ! 1: c/a, 2: d/a, 3: r/a, 4: e/a
      real(PS),dimension(mxnbin+1,n1,mxnaxis-1) :: axr_p
      ! bulk sphere density of dry ice particle, and bulk crystal density in the shifted bin.
      real(PS), dimension(mxnbin+1,n1)           :: den_ip_p,den_ic_p
      ! habit in shifted bin
      integer, dimension(mxnbin+1,n1)          :: habit_p
      ! aspect ratio of circumscribing cylinder (not used)
      real(PS)              :: spx_p

      ! aspect ratio of ice particle
      real(PS), dimension(mxnbin+1,n1)              :: asr_p
      ! type in shifted bin
      integer, dimension(mxnbin+1,n1)          :: type_p


      ! ratio of ag^3 to a^3
      real(PS), dimension(mxnbin+1,n1)           :: rag_p,rcg_p
      ! number of extra ice crystals
      real(RP), dimension(mxnbin+1,n1)           :: n_exice_p

      ! activated IN fraction for contact parameter diagnosis
      real(PS), dimension(mxnbin+1,n1)           :: actINF_p


      integer,dimension(mxnbin+1,n1)           :: error_number

      integer,dimension(mxnbin+1,n1)           :: ierror1

      real(PS),dimension(mxnbin,n1) :: mmass_p
      real(PS),dimension(mxnbin,n1) :: rmod
      real(PS) :: v_sp
      real(8),dimension(n1) :: ap_dN,ap_dM
!       4*pi/3
      real(PS) :: coef3
      parameter(coef3=4.18879020478639)
!     error message
      integer :: em,icem,ncem
      integer :: ia,im,iq,k,ibr,icr,ik,j,jk
      real(PS),dimension(n1) :: den1
!     minimum fraction of aerosols in hydrometeors
      real(RP) :: min_fapt_r,min_faps_r,min_fapt_s,min_faps_s
      parameter(min_fapt_r=1.0e-18_RP,min_faps_r=1.0e-6_RP &
               ,min_fapt_s=1.0e-18_RP,min_faps_s=0.0_RP)

!     possible ratio of mean mass to bin boundaries
      real(RP) :: brat1,brat2
      parameter(brat1=1.001_RP,brat2=0.999_RP)

!     factor
      real(RP) :: frac,fct1,fct2
      parameter(frac=0.1_RP,fct1=0.1_RP,fct2=0.99_RP)

      real(PS) :: dum,max_rad,rat

!     minimum possible volume of ice crystalc
      real(RP),parameter :: V_csmin=1.1847688e-11_RP

      real(PS) :: sum_m1(max(nprmx,npimx,npamx),n1) &
             ,sum_m2(max(nprmx,npimx,npamx),n1)

      integer,dimension(n1) :: icond1
      integer,dimension(mxnbin*n1) :: icond2,icond3

      real(PS),parameter :: amin_mmass=1.0e-25
      ! bin boundary modification factor.
      real(PS),parameter    :: bbmf=0.2
!      write(fid_alog,*) "I am in shifted bin",npr,nbr,ncr,iphase

!!!      write(fid_alog,*) "k1r,k2r",k1r,k2r

!!!      return

      if(k1r>k2r) return

!   initialization

      icem=0
      ncem=0


      do k=1,n1
        den1(k)=den(k)*1.0e-3_PS

        !  these are not used.
        ap_dN(k)=0.0d+0
        ap_dM(k)=0.0d+0

        if(k>=k1r.and.k<=k2r) then
          icond1(k)=1
        else
          icond1(k)=0
        endif
      enddo

      if(iphase==1) then
!
!   in case of liquid spectrum
!
!

        do icr=1,ncr
!
!     initialization
!
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr
            new_N(ibr,k)=0.0_DS
            new_mtend(ibr,k)=0.0_DS
            rmod(ibr,k)=1.0
          enddo
          do j=1,nvar_mcp_liq+1
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_M(ibr,k,j)=0.0_DS
            end do
          enddo
          do j=1,nvar_nonmcp_liq
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_Q(ibr,k,j)=0.0_DS
            end do
          enddo

          do jk=1,(npr-2)*n1
            k=(jk-1)/(npr-2)+1
            j=jk-(k-1)*(npr-2)
            sum_m1(j,k)=0.0_PS
            sum_m2(j,k)=0.0_PS
          end do

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m1(j,k)=sum_m1(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

          !
          ! initialize
          !
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond1(k)==1.and.&
               qrpv(rmt_q,ibr,icr,k)>1.0e-30_RP.and. &
               qrpv(rcon_q,ibr,icr,k)>1.0e-30_RP) then
              icond3(ik)=1
              mmass_p(ibr,k)=qrpv(rmt_q,ibr,icr,k)/qrpv(rcon_q,ibr,icr,k)
            else
              icond3(ik)=0
              mmass_p(ibr,k)=0.0_RP
            endif
          enddo

!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then


!     calculate mass components ratio
!         for aerosols
              ratio_Mp(ibr,k,rmat_m)=min( &
                 max(qrpv(rmat_q,ibr,icr,k)/qrpv(rmt_q,ibr,icr,k) &
                        ,min_fapt_r),1.0_RP)

              ratio_Mp(ibr,k,rmas_m)=min( &
                 max(qrpv(rmas_q,ibr,icr,k)/qrpv(rmt_q,ibr,icr,k) &
                        ,min_faps_r),ratio_Mp(ibr,k,rmat_m))

              actINF_p(ibr,k)=0.0_RP

              binb3d(ibr,k,1)=max(0.0_RP,binbr(ibr) &
                                 +mmass_p(ibr,k)-mmass(ibr,icr,k))
              binb3d(ibr,k,2)=max(0.0_RP,binbr(ibr+1) &
                                 +mmass_p(ibr,k)-mmass(ibr,icr,k))


              Np(ibr,k)=qrpv(rcon_q,ibr,icr,k)*den1(k)
              Mp(ibr,k)=qrpv(rmt_q,ibr,icr,k)*den1(k)
              mtend(ibr,k)=0.0_RP
            else
              ratio_Mp(ibr,k,rmat_m)=0.0_RP
              ratio_Mp(ibr,k,rmas_m)=0.0_RP
              actINF_p(ibr,k)=0.0_RP

              binb3d(ibr,k,1)=binbr(ibr)
              binb3d(ibr,k,2)=binbr(ibr+1)
!org              Np(ibr,k)=0.0_RP
!org              Mp(ibr,k)=0.0_RP
              Np(ibr,k)=max(0.0_RP,qrpv(rcon_q,ibr,icr,k)*den1(k))
              Mp(ibr,k)=max(0.0_RP,qrpv(rmt_q,ibr,icr,k)*den1(k))
              mtend(ibr,k)=0.0_RP
            endif
          enddo

          call cal_lincubprms_vec(mxnbin+1,nbr,n1,Np,Mp,binb3d  &
                              ,a2d,error_number,"shift_liq")

...  (debug block elided-in-place; original has commented writes) ...

          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr
            if(1<=error_number(ibr,k).and.error_number(ibr,k)<=4) then
              icem=icem+1
              ncem=ncem+1
            elseif(error_number(ibr,k)==0) then
              ncem=ncem+1
            endif
          enddo

          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.&
               ( (1<=error_number(ibr,k).and.error_number(ibr,k)<=4).or.&
                 mmass(ibr,icr,k)<amin_mmass) ) then
              icond2(ik)=1
            else
              icond2(ik)=0
            endif
          enddo

!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond2(ik)==1) then

              new_M(ibr,k,rmt)=new_M(ibr,k,rmt) &
                                +qrpv(rmt_q,ibr,icr,k)*den1(k)


              if(mmass_p(ibr,k)>binbr(ibr+1) &
                                .or.mmass_p(ibr,k)<binbr(ibr))then
                rmod(ibr,k)=qrpv(rmt_q,ibr,icr,k) &
                      /max(brat1*binbr(ibr) &
                      ,min(brat2*binbr(ibr+1),mmass_p(ibr,k))) &
                      /qrpv(rcon_q,ibr,icr,k)
              else
                rmod(ibr,k)=1.0
              end if

              new_N(ibr,k)=new_N(ibr,k) &
                             +qrpv(rcon_q,ibr,icr,k)*den1(k)*rmod(ibr,k)


              ! this make these cases to skip in trans_bin
              error_number(ibr,k)=10
            endif
          enddo

!     for mass components
          do iq=i1_mcp_liq,i2_mcp_liq
            im=rmt+1+iq-i1_mcp_liq
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then

                new_M(ibr,k,im)=new_M(ibr,k,im) &
                             +qrpv(iq,ibr,icr,k)*den1(k)
              endif
            enddo
          enddo

!     for concentration component
          do iq=i1_ccp_liq,i2_ccp_liq
            im=iq-i1_ccp_liq+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                         +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            end do
          end do

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! calculation of transferred concentration and mass into original bins
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          call cal_transbin_vec(iphase &
                           ,n1,nvar_mcp_liq &
                           ,nbr,nbr &
                           ,binbr &
                           ,error_number &
                           ,a2d,binb3d,mtend &
                           ,new_N,new_M,new_Q &
                           ,new_mtend &
                           ,ratio_Mp,den_ip_p,axr_p,spx_p &
                           ,habit_p,den_ic_p &
                           ,rag_p,rcg_p,n_exice_p &
                           ,actINF_p &
                           ,0)
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     put new variables back into the originals
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond1(k)==1) then
              qrpv(rmt_q,ibr,icr,k)=max(new_M(ibr,k,rmt)/den1(k),0.0_DS)
              qrpv(rcon_q,ibr,icr,k)=max(new_N(ibr,k)/den1(k),0.0_DS)
            endif
          enddo

! mass component
          do iq=i1_mcp_liq,i2_mcp_liq
            im=rmt+1+iq-i1_mcp_liq
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)=max(new_M(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

!   for concentration component
          do iq=i1_ccp_liq,i2_ccp_liq
            im=iq-i1_ccp_liq+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m2(j,k)=sum_m2(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

          if(debug) then
            do k=1,n1
              if(sum_m1(1,k)>1.0e-30) then
                if(abs((sum_m2(1,k)-sum_m1(1,k))/sum_m1(1,k))>1.0e-4) then
                  write(fid_alog,*) "bf liq:sum", &
                    kmic(k),imic(k),jmic(k),sum_m1(1:npr-2,k)
                  write(fid_alog,*) "af liq:sum", &
                    kmic(k),imic(k),jmic(k),sum_m2(1:npr-2,k)
                  write(fid_alog,*) "np",(np(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "mp",(mp(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "error_number",(error_number(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "mmass",(mmass(ibr,icr,k),ibr=1,nbr)
                  write(fid_alog,*) "mmass_p",(mmass_p(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "rmod",(rmod(ibr,k),ibr=1,nbr)
                  write(fid_alog,*) "a2d(1)",(a2d(ibr,k,1),ibr=1,nbr)
                  write(fid_alog,*) "a2d(2)",(a2d(ibr,k,2),ibr=1,nbr)
                  write(fid_alog,*) "a2d(3)",(a2d(ibr,k,3),ibr=1,nbr)
                end if
              endif
            enddo
          endif

!     calculate mean mass for cloud droplet bin
!CDIR NODEP
          do k=1,n1
            mmass(1,icr,k)= &
                qrpv(rmt_q,1,icr,k)/max(1.0e-30_RP,qrpv(rcon_q,1,icr,k))
          enddo
        end do

        if(debug) then
          if(ncem>0) then
            write(fid_alog,301) 100.0_PS*real(icem,PS_KIND)/real(ncem,PS_KIND)
301   format("error % in liq shift_bin:",ES15.6)
          end if
        end if



      elseif(iphase==2) then
!
!  in case of ice spectrum
!

        do icr=1,ncr
!
!     initialization
!
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr
            new_N(ibr,k)=0.0_PS
            new_mtend(ibr,k)=0.0_PS
            rmod(ibr,k)=1.0
          enddo
          do j=1,nvar_mcp_ice+1  ! including total mass
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_M(ibr,k,j)=0.0_PS
            end do
          enddo
          do j=1,nvar_nonmcp_ice
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              new_Q(ibr,k,j)=0.0_PS
            end do
          enddo

          do jk=1,(npr-2)*n1
            k=(jk-1)/(npr-2)+1
            j=jk-(k-1)*(npr-2)
            sum_m1(j,k)=0.0
            sum_m2(j,k)=0.0
          end do

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m1(j,k)=sum_m1(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

          !
          ! initialize
          !
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            mmass_p(ibr,k)=0.0
            error_number(ibr,k)=0

            if(icond1(k)==1.and.&
               qrpv(imt_q,ibr,icr,k)>1.0e-30.and. &
               qrpv(icon_q,ibr,icr,k)>1.0e-30) then

              icond3(ik)=1
              mmass_p(ibr,k)=qrpv(imt_q,ibr,icr,k)/qrpv(icon_q,ibr,icr,k)
            else
              icond3(ik)=0
              mmass_p(ibr,k)=0.0

            endif
          enddo

!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then
!
!     calculate mass components ratio
!
!         for ice crystal mass
              ratio_Mp(ibr,k,imc_m)=min( &
                      max(max( &
                     qrpv(icon_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k)*4.763209003e-12 &
                    ,qrpv(imc_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k)) &
                      ,0.0_RP),1.0_RP)
!         for rimed mass
              ratio_Mp(ibr,k,imr_m)=min( &
                      max(qrpv(imr_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),1.0_RP)

!         for melt water mass
              ratio_Mp(ibr,k,imw_m)=min( &
                      max(qrpv(imw_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),1.0_RP)

!         for aggregation mass
              ratio_Mp(ibr,k,ima_m)=min( &
                      max(qrpv(ima_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),1.0_RP)

!         for freezing nucleation mass
              ratio_Mp(ibr,k,imf_m)=min( &
                      max(qrpv(imf_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,0.0_RP),ratio_Mp(ibr,k,imc_m))

!         normalization
              dum = &
                ratio_Mp(ibr,k,imr_m)+ratio_Mp(ibr,k,ima_m)+ratio_Mp(ibr,k,imw_m)

              if(dum>1.0e-20_RP) then
                ratio_Mp(ibr,k,imr_m)= &
                      (1.0_RP-ratio_Mp(ibr,k,imc_m))*ratio_Mp(ibr,k,imr_m)/dum

                ratio_Mp(ibr,k,ima_m)= &
                      (1.0_RP-ratio_Mp(ibr,k,imc_m))*ratio_Mp(ibr,k,ima_m)/dum

                ratio_Mp(ibr,k,imw_m)= &
                      (1.0_RP-ratio_Mp(ibr,k,imc_m))*ratio_Mp(ibr,k,imw_m)/dum
              else
                ratio_Mp(ibr,k,imc_m)=1.0_RP
                ratio_Mp(ibr,k,imr_m)=0.0_RP
                ratio_Mp(ibr,k,ima_m)=0.0_RP
                ratio_Mp(ibr,k,imw_m)=0.0_RP
              end if

              ratio_Mp(ibr,k,imf_m)=min(ratio_Mp(ibr,k,imf_m) &
                                       ,ratio_Mp(ibr,k,imc_m))

!         for aerosol total mass
              ratio_Mp(ibr,k,imat_m)=min( &
                      max(qrpv(imat_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,min_fapt_s),1.0_RP)

!         for soluble aerosol mass
              ratio_Mp(ibr,k,imas_m)=min( &
                      max(qrpv(imas_q,ibr,icr,k)/qrpv(imt_q,ibr,icr,k) &
                     ,min_faps_s),ratio_Mp(ibr,k,imat_m))

!         for circumscribing volume
              Qp(ibr,k,ivcs)= &
                      max(qrpv(ivcs_q,ibr,icr,k)/qrpv(icon_q,ibr,icr,k) &
                                 ,mmass_p(ibr,k)/den_i)

            else
!         for crystal mass
              ratio_Mp(ibr,k,imc_m)=0.0_RP
!         for rimed mass
              ratio_Mp(ibr,k,imr_m)=0.0_RP
!         for melt water mass
              ratio_Mp(ibr,k,imw_m)=0.0_RP
!         for aggregation mass
              ratio_Mp(ibr,k,ima_m)=0.0_RP
!         for freezing nucleation mass
              ratio_Mp(ibr,k,imf_m)=0.0_RP

!         for aerosol total mass
              ratio_Mp(ibr,k,imat_m)=0.0_RP
!         for soluble aerosol mass
              ratio_Mp(ibr,k,imas_m)=0.0_RP

!         for circumscribing volume
              Qp(ibr,k,ivcs)=0.0_RP

            endif
          enddo


!         for axis component
          do iq=i1_acp_ice,i2_acp_ice
            im=nvar_vcp_ice+iq-i1_acp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then
                Qp(ibr,k,im)=(max(qrpv(iq,ibr,icr,k)/ &
                        qrpv(icon_q,ibr,icr,k),0.0_RP))**(1.0_RP/3.0_RP)
              else
                Qp(ibr,k,im)=0.0_RP
              endif
            enddo
          enddo

!         for concentration component
          do iq=i1_ccp_ice,i2_ccp_ice
            im=nvar_vcp_ice+nvar_acp_ice+iq-i1_ccp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then
                Qp(ibr,k,im)=qrpv(iq,ibr,icr,k)/ &
                       qrpv(icon_q,ibr,icr,k)
              else
                Qp(ibr,k,im)=0.0_RP
              endif
            enddo
          enddo

          ierror1=0
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.mmass(ibr,icr,k)>=amin_mmass) then

!         Realitiy check
              if(Qp(ibr,k,ivcs)<=1.1847688E-11) then
                Qp(ibr,k,ivcs)=1.1847688E-11
              end if
              if(Qp(ibr,k,iacr)<1.0e-4) then
                Qp(ibr,k,iacr)=1.0e-4
              end if
              if(Qp(ibr,k,iccr)<1.0e-4) then
                Qp(ibr,k,iccr)=1.0e-4
              end if

              Qp(ibr,k,idcr)=min(Qp(ibr,k,iacr),max(0.0_RP,Qp(ibr,k,idcr)))
              if(Qp(ibr,k,iacr)<=1.0e-4_RP) then
                Qp(ibr,k,idcr)=0.0_RP
              end if
              Qp(ibr,k,iag)=max(0.0_RP,Qp(ibr,k,iag))
              Qp(ibr,k,icg)=max(0.0_RP,Qp(ibr,k,icg))
              Qp(ibr,k,inex)=min(1.0_RP,max(0.0_RP,Qp(ibr,k,inex)))

              ! c and d length ratio to the a length
              ia=1
              axr_p(ibr,k,ia)=Qp(ibr,k,iacr+ia)/Qp(ibr,k,iacr)
              ia=2
              axr_p(ibr,k,ia)=Qp(ibr,k,iacr+ia)/Qp(ibr,k,iacr)

!     pristine crystal
!         ratio for center of gravities
              rag_p(ibr,k)=(Qp(ibr,k,iag)/Qp(ibr,k,iacr))**3
              rcg_p(ibr,k)=(Qp(ibr,k,icg)/Qp(ibr,k,iccr))**3

!         number of extra crystals
              n_exice_p(ibr,k)=Qp(ibr,k,inex)
!         activated IN concentration fraction
              actINF_p(ibr,k)=0.0_RP

!                          hexagonal crystals
              v_sp=coef3*((1.0+axr_p(ibr,k,1)**2)**1.5_RP)  &
                                     *Qp(ibr,k,iacr)**3
              habit_p(ibr,k)=1
              spx_p=axr_p(ibr,k,1)

              den_ic_p(ibr,k)=max(1.0e-4_RP,&
                       mmass_p(ibr,k)*ratio_Mp(ibr,k,imc_m)/v_sp)

              axr_p(ibr,k,3)=0.0_RP
              axr_p(ibr,k,4)=0.0_RP

              if(den_ic_p(ibr,k)>den_i) then
                v_sp=mmass_p(ibr,k)*ratio_Mp(ibr,k,imc_m)/den_i

                den_ic_p(ibr,k)=den_i

                Qp(ibr,k,iacr)=(v_sp/ &
                     (coef3*((1.0_RP+spx_p**2)**1.5_RP)))**(1.0_RP/3.0_RP)
                Qp(ibr,k,iccr)=Qp(ibr,k,iacr)*spx_p
              end if

              Qp(ibr,k,ivcs)=max(Qp(ibr,k,ivcs),V_csmin)

              den_ip_p(ibr,k)=mmass_p(ibr,k)/Qp(ibr,k,ivcs)

              binb3d(ibr,k,1)=max(0.0_RP,binbr(ibr) &
                             +mmass_p(ibr,k)-mmass(ibr,icr,k))
              binb3d(ibr,k,2)=max(0.0_RP,binbr(ibr+1) &
                             +mmass_p(ibr,k)-mmass(ibr,icr,k))

              Np(ibr,k)=qrpv(icon_q,ibr,icr,k)*den1(k)
              Mp(ibr,k)=qrpv(imt_q,ibr,icr,k)*den1(k)
              mtend(ibr,k)=0.0_RP

            else

              ! c and d length ratio to the a length
              ia=1
              axr_p(ibr,k,ia)=1.0_RP
              ia=2
              axr_p(ibr,k,ia)=1.0_RP

!     pristine crystal
!         ratio for center of gravities
              rag_p(ibr,k)=0.0_RP
              rcg_p(ibr,k)=0.0_RP

!         number of extra crystals
              n_exice_p(ibr,k)=0.0_RP
!         activated IN concentration fraction
              actINF_p(ibr,k)=0.0_RP

!                          hexagonal crystals
              habit_p(ibr,k)=0
              spx_p=0.0

              den_ic_p(ibr,k)=1.0_RP

              axr_p(ibr,k,3)=0.0_RP
              axr_p(ibr,k,4)=0.0_RP

              den_ip_p(ibr,k)=1.0_RP


              binb3d(ibr,k,1)=binbr(ibr)
              binb3d(ibr,k,2)=binbr(ibr+1)
              Np(ibr,k)=max(0.0_RP,qrpv(icon_q,ibr,icr,k)*den1(k))
              Mp(ibr,k)=max(0.0_RP,qrpv(imt_q,ibr,icr,k)*den1(k))
              mtend(ibr,k)=0.0_RP

            endif

            if(den_ic_p(ibr,k)<1.0e-30_RP) then
              ierror1(ibr,k)=1
            endif
          enddo
          if(any(ierror1>0)) then
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              if(debug .and. ierror1(ibr,k)==1) then
                write(fid_alog,*) "Shift_bin_vec: Error1:",k,ibr,icond3(ik),den_ic_p(ibr,k),np(ibr,k),mp(ibr,k),&
                      mmass(ibr,icr,k),mmass_p(ibr,k),qrpv(imt_q,ibr,icr,k),qrpv(icon_q,ibr,icr,k)
              endif
            enddo
            stop
          endif

          call cal_lincubprms_vec(mxnbin+1,nbr,n1,Np,Mp,binb3d  &
                              ,a2d,error_number,"shift_ice")


          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(1<=error_number(ibr,k).and.error_number(ibr,k)<=4) then
              icem=icem+1
              ncem=ncem+1
            elseif(error_number(ibr,k)==0) then
              ncem=ncem+1
            endif
          enddo

          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond3(ik)==1.and.&
               ( (1<=error_number(ibr,k).and.error_number(ibr,k)<=4).or.&
                 mmass(ibr,icr,k)<amin_mmass) ) then
              icond2(ik)=1
            else
              icond2(ik)=0
            endif
          enddo

!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond2(ik)==1) then

              new_M(ibr,k,imt)=new_M(ibr,k,imt) &
                           +qrpv(imt_q,ibr,icr,k)*den1(k)

              if(mmass_p(ibr,k)>binbr(ibr+1) &
                 .or.mmass_p(ibr,k)<binbr(ibr))then

                rmod(ibr,k)=qrpv(imt_q,ibr,icr,k) &
                     /max(brat1*binbr(ibr) &
                     ,min(brat2*binbr(ibr+1),mmass_p(ibr,k))) &
                     /qrpv(icon_q,ibr,icr,k)
              else
                rmod(ibr,k)=1.0_RP
              end if

              new_N(ibr,k)=new_N(ibr,k) &
                       +qrpv(icon_q,ibr,icr,k)*den1(k)*rmod(ibr,k)

              ! this make these cases to skip in trans_bin
              error_number(ibr,k)=10
            endif
          enddo

!     for volume component
          do iq=i1_vcp_ice,i2_vcp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                      +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            end do
          enddo

!     for axis component
          do iq=i1_acp_ice,i2_acp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                      +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            enddo
          enddo

!     for concentration component
          do iq=i1_ccp_ice,i2_ccp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_Q(ibr,k,im)=new_Q(ibr,k,im) &
                         +qrpv(iq,ibr,icr,k)*den1(k)*rmod(ibr,k)
              endif
            end do
          end do

          do iq=i1_mcp_ice,i2_mcp_ice
            im=imt+1+iq-i1_mcp_ice
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond2(ik)==1) then
                new_M(ibr,k,im)=new_M(ibr,k,im) &
                                  +qrpv(iq,ibr,icr,k)*den1(k)
              endif
            enddo
          enddo

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! calculation of transferred concentration and mass into original bins
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          call cal_transbin_vec(iphase &
                           ,n1,nvar_mcp_ice &
                           ,nbr,nbr &
                           ,binbr &
                           ,error_number &
                           ,a2d,binb3d,mtend &
                           ,new_N,new_M,new_Q &
                           ,new_mtend &
                           ,ratio_Mp,den_ip_p,axr_p,spx_p &
                           ,habit_p,den_ic_p &
                           ,rag_p,rcg_p,n_exice_p &
                           ,actINF_p &
                           ,0)
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            ierror1(ibr,k)=0
            if(error_number(ibr,k)/=10) then
              if(den_ic_p(ibr,k)<1.0e-30) then
                ierror1(ibr,k)=1
              endif
            endif
          enddo
          if(any(ierror1>0)) then
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr
              if(debug .and. ierror1(ibr,k)==1) then
                write(fid_alog,*) "Shift_bin_vec: Error3:",k,ibr,icond2(ik),icond3(ik),mmass(ibr,icr,k),mmass_p(ibr,k),&
                           den_ic_p(ibr,k),np(ibr,k),mp(ibr,k),&
                           error_number(ibr,k)
              endif
            enddo
            stop
          endif

!     put new variables back into the originals
!CDIR NODEP
          do ik=1,nbr*n1
            k=(ik-1)/nbr+1
            ibr=ik-(k-1)*nbr

            if(icond1(k)==1) then
! total mass
              qrpv(imt_q,ibr,icr,k)=max(new_M(ibr,k,imt)/den1(k),0.0_DS)
! concentration
              qrpv(icon_q,ibr,icr,k)=max(new_N(ibr,k)/den1(k),0.0_DS)

            endif
          enddo

! mass component
          do iq=i1_mcp_ice,i2_mcp_ice
            im=imt+1+iq-i1_mcp_ice
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)=max(new_M(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo


!  for volume component
          do iq=i1_vcp_ice,i2_vcp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            end do
          end do

!   for axis component
          do iq=i1_acp_ice,i2_acp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

!   for concentration component
          do iq=i1_ccp_ice,i2_ccp_ice
            im=iq-i1_vcp_ice+1
!CDIR NODEP
            do ik=1,nbr*n1
              k=(ik-1)/nbr+1
              ibr=ik-(k-1)*nbr

              if(icond1(k)==1) then
                qrpv(iq,ibr,icr,k)= &
                      max(new_Q(ibr,k,im)/den1(k),0.0_DS)
              endif
            enddo
          enddo

          do ibr=1,nbr
!CDIR NODEP
            do jk=1,(npr-2)*n1
              k=(jk-1)/(npr-2)+1
              j=jk-(k-1)*(npr-2)
              sum_m2(j,k)=sum_m2(j,k)+qrpv(j,ibr,icr,k)*den1(k)
            end do
          end do

        end do

      end if
      end subroutine shift_bin_vec
```

(Two large fully-commented `!tmp`/`!!c` diagnostic blocks at 6838-6841, 6976, 7196-7197, 7271-7281, 7537-7569, 7697-7717 are elided in the quote above where marked "…"; they are pure comments with no executable effect.)

---

## 4. Repair + budget closure — `mod_amps_check.F90`

### 4a. Repair driver structure — `contrib/AMPS/mod_amps_check.F90:22-1181`

Signature + key declarations + iteration limits:

```fortran
  subroutine repair(level,ag,gr,gs,ncat_a,ga,flagp_r,flagp_s,flagp_a,&
                    mes_rc, & !act_type,
                    irep_vap,irep_col,&
                    iupdate_gr,iupdate_gs,iupdate_ga,&
                    qtp,ID,JD,KD,from,itcall)
    ...
    type (AirGroup),intent(inout)    :: ag
    type (Group), intent(inout)         :: gr, gs
    integer,intent(in) :: ncat_a
    type (Group), dimension(ncat_a)         :: ga
    ...
    ! accumulated (multiplied over the iteration) modification constant
    real (PS), dimension(mxnbin,mxntend,LMAX)    :: acc_mod_r
    real (PS), dimension(mxnbin,mxntend,LMAX)    :: acc_mod_s
    real (PS), dimension(mxnbina,mxntend,ncamx,LMAX)    :: acc_mod_a
    ! mark of modification
    integer,dimension(LMAX)   :: mark_mod
    integer,dimension(LMAX)   :: mark_all
    integer,dimension(LMAX)   :: mark_mass_r, mark_mass_s
    integer,dimension(LMAX)   :: mark_con_r, mark_con_s
    integer,dimension(LMAX)   :: mark_v_r, mark_v_s
    integer,dimension(ncamx,LMAX)    :: mark_mass_a,mark_con_a,mark_v_a
    ! maximum number of iteration
    integer                   :: NITER, NITER_TOTAL
    ...
    logical,parameter :: ifix_mass_col=.true.
    logical,parameter :: ifix_con_col=.true.
    logical,parameter :: ifix_vol_col=.true.

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! set up the maximum number of iteration
    NITER = 300
    NITER_TOTAL = 4
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```

The driver runs two independent repair phases:

**Phase 1 — vapor** (`if(irep_vap)`, lines 329-431): initialize `acc_mod_*=1`, `mark_*=0`, then the inner convergence loop:

```fortran
      do itr1 = 1, NITER
        ! calculate the budget of each group and modify the mass tendency
        call cal_mass_budget_vapor( ag, gr, gs, ncat_a,ga, mes_rc, &
               acc_mod_r, acc_mod_s, acc_mod_a, &
               mark_mod, &
               mark_mass_r, mark_mass_s, mark_mass_a, flagp_r, flagp_s, flagp_a,&
               ID,JD,KD)

        if( any(mark_mod(1:ag%L)>0) )  then
        else
          exit
        end if
      end do
      if( any(mark_mod(1:ag%L) > 0) ) then    ! non-convergence => abort
        ...
             call PRC_abort
      end if
```
then `mod_other_tendency_vap` is applied to whichever groups were marked (`any(mark_mass_r/s/a /= 0)`).

**Phase 2 — collision** (`if(irep_col)`, lines 445-772): the outer `NITER_TOTAL` loop wrapping three sub-blocks (mass, con, vol), each initializing `acc_mod`/`mark`, running its own `do itr1=1,NITER` any()-guarded convergence loop, then calling `mod_other_tendency` to propagate the accumulated modification into the other tendencies:

```fortran
      do itr2 = 1, NITER_TOTAL

      if(ifix_mass_col) then
        ...  ! reinit acc_mod_r/s/a=1, mark_mass_*=0
        do itr1 = 1, NITER
          call cal_mass_budget_col( ag, gr, gs, ncat_a,ga, &
                  acc_mod_r, acc_mod_s, acc_mod_a, &
                  mark_mod, &
                  mark_mass_r, mark_mass_s, mark_mass_a, flagp_r, flagp_s, flagp_a,&
                  iupdate_gr,iupdate_gs,iupdate_ga)
          if( any(mark_mod(1:ag%L) > 0 )) then
          else
            exit
          end if
        end do
        if( any(mark_mod(1:ag%L)>0) )then   ! non-convergence => abort
          ... call PRC_abort
        end if
        if( any(mark_mass_r(1:gr%L) /= 0 )) call mod_other_tendency(gr, acc_mod_r, 1, iupdate_gr)
        if( any(mark_mass_s(1:gs%L) /= 0 )) call mod_other_tendency(gs, acc_mod_s, 1, iupdate_gs)
        if(level>=4) then
          do ica=1,ncat_a
            if( any(mark_mass_a(ica,1:ga(ica)%L) /= 0) ) call mod_other_tendency_ap( ga, acc_mod_a, ica,1, iupdate_ga)
          end do
        end if
      endif

      if(ifix_con_col) then
        ...  ! reinit, do itr1=1,NITER: call cal_con_budget_col(...)  any()-guarded
        if( any(mark_con_r(1:gr%L) /= 0 )) call mod_other_tendency(gr, acc_mod_r, 2, iupdate_gr)
        if( any(mark_con_s(1:gs%L) /= 0 )) call mod_other_tendency(gs, acc_mod_s, 2, iupdate_gs)
        if(level>=4) then ... mod_other_tendency_ap( ga, acc_mod_a, ica,2, iupdate_ga) ... end if
      endif

      if(ifix_vol_col) then
        do nvol=1,1
          ...  ! reinit, do itr1=1,NITER: call cal_vol_budget_col(...)  any()-guarded
          if( any(mark_v_r(1:gr%L) /= 0 )) call mod_other_tendency(gr, acc_mod_r, 2+nvol, iupdate_gr)
          if( any(mark_v_s(1:gs%L) /= 0 )) call mod_other_tendency(gs, acc_mod_s, 2+nvol, iupdate_gs)
        end do
      endif

        mark_all=mark_mass_r+mark_mass_s+ &
                 mark_mass_a(1,:)+mark_mass_a(2,:) +&
                 mark_con_a(1,:)+mark_con_a(2,:)   +&
                 mark_con_r + mark_con_s +&
                 mark_v_a(1,:) + mark_v_a(2,:) + &
                 mark_v_r + mark_v_s
        if(sum(mark_all)==0) then
          exit
        end if
      end do
    endif
```

So the convergence structure is: inner `do itr1=1,NITER(=300)` per budget, exiting early on `.not.any(mark_mod>0)`; outer `do itr2=1,NITER_TOTAL(=4)` over the coupled mass/con/vol budgets, exiting when `sum(mark_all)==0`. The `any()`/`sum()` global reductions over `1:ag%L` are the host-loop reduction targets for M2.

### 4b. Representative `liqbin_loop` repair block — `cal_mass_budget_col`, `contrib/AMPS/mod_amps_check.F90:1261-1466`

Relevant local decls (from the subroutine header, 1186-1244): `total_rloss, total_rgain` scalars, `modc_r(LMAX)` modification constant, `real(PS),parameter :: total_limit=1.0e-30_PS`, `icond1/icond2(LMAX)`. The mass-conservation rescale computes total gain vs total loss per column, forms `modc_r = min(1, gain/max(loss,total_limit))`, and multiplies each losing process tendency (and its `acc_mod`) by `modc_r`:

```fortran
      liqbin_loop1: do i=1,gr%n_bin
        do n=1,gr%L

          total_rloss = &
               ! - collision advection
               max(-gr%MS(i,n)%dmassdt(rmt,2),0.0_ds)*real(iupdate_gr(2),PS_KIND) + &
               ! - riming process by a solid hyrometeor
               max(-gr%MS(i,n)%dmassdt(rmt,4), 0.0_ds)*real(iupdate_gr(4),PS_KIND) + &
               ! - collision breakup
               max(-gr%MS(i,n)%dmassdt(rmt,5),0.0_ds)*real(iupdate_gr(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(-gr%MS(i,n)%dmassdt(rmt,6),0.0_ds)*real(iupdate_gr(6),PS_KIND) + &
               ! - loss of cloud droplets by contact nucleation
               max(-gr%MS(i,n)%dmassdt(rmt,8),0.0_ds)*real(iupdate_gr(8),PS_KIND) + &
               ! - auto conversion
               max(-gr%MS(i,n)%dmassdt(rmt,9), 0.0_ds)*real(iupdate_gr(9),PS_KIND) + &
               ! - ice nucleation (immersion)
               max(-gr%MS(i,n)%dmassdt(rmt,11), 0.0_ds)*real(iupdate_gr(11),PS_KIND) + &
               ! - ice nucleation (homogeneous freezing)
               max(-gr%MS(i,n)%dmassdt(rmt,12), 0.0_ds)*real(iupdate_gr(12),PS_KIND)

          total_rgain = &
               ! - existing mass
               gr%MS(i,n)%mass(rmt)/gr%dt + &
               ! - collision advection
               max(gr%MS(i,n)%dmassdt(rmt,2),0.0_ds)*real(iupdate_gr(2),PS_KIND) + &
               ! - collision breakup
               max(gr%MS(i,n)%dmassdt(rmt,5),0.0_ds)*real(iupdate_gr(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(gr%MS(i,n)%dmassdt(rmt,6),0.0_ds)*real(iupdate_gr(6),PS_KIND) + &
               ! - auto conversion
               max(gr%MS(i,n)%dmassdt(rmt,9), 0.0_ds)*real(iupdate_gr(9),PS_KIND) + &
               ! - melting-shedding
               max(gr%MS(i,n)%dmassdt(rmt,10),0.0_ds)*real(iupdate_gr(10),PS_KIND)

          icond1(n)=0
          modc_r(n)=1.0_PS

          if( total_rgain < 0.99999_ps*total_rloss.and.total_rloss>total_limit) then

            icond1(n)=1

            modc_r(n) = min(1.0_ps,&
                      (total_rgain)/max(total_rloss,total_limit))

            modc_r(n)=max(0.0_PS,modc_r(n))

            mark_mod(n) = mark_mod(n) + 1
            mark_acc_r(n) = mark_acc_r(n) + 1

          else if( total_rgain == 0.0_ps .and. total_rloss == 0.0_ps ) then
            gr%MS(i,n)%mark = 4
          end if
          if(modc_r(n)>1.0e+4) then
             LOG_ERROR("cal_mass_budget_col",*) "something wrong",n,gr%MS(i,n)%dmassdt(rmt,1:gr%n_tendpros)
             call prc_abort
          endif
        enddo

        if(all(icond1(1:gr%L)==0)) cycle

        do j=1,gr%n_tendpros
          if(iupdate_gr(j)==0) cycle

          do n=1,gr%L
            icond2(n)=0
            if( gr%MS(i,n)%dmassdt(rmt,j) < 0.0_ps ) then
              icond2(n)=1
            endif
          enddo

          if(j == 2 .or. j == 5 .or. j == 6.or.j==9) then
            ! collision-coalescence
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if(icond2(n)==1) then
                gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
          elseif(j==4.and.flagp_s>0) then
            ! riming process
            do n=1,gr%L
              if( icond2(n)==1) then
                gr%MS(i,n)%dmassdt(rmt,j)=gr%MS(i,n)%dmassdt(rmt,j)*modc_r(n)
                acc_mod_r(i,j,n)=acc_mod_r(i,j,n)*modc_r(n)
              endif
              total_rimr_tend(n)=0.0_ps
              total_rims_tend(n)=0.0_ps
            enddo
            do k=1,gr%N_bin
              do n=1,gr%L
                total_rimr_tend(n)=total_rimr_tend(n)-gr%MS(k,n)%dmassdt(rmt,j)
              enddo
            enddo
            do k=1,gs%N_bin
              do n=1,gs%L
                total_rims_tend(n)=total_rims_tend(n)+gs%MS(k,n)%dmassdt(imt,j)
              enddo
            enddo
            do n=1,gr%L
              mod_dum(n)=1.0_PS
              if( icond2(n)==1) then
                mod_dum(n)=max(0.0_PS,min(1.0_PS,total_rimr_tend(n)/max(total_limit,total_rims_tend(n))))
              endif
              mark_acc_s(n)=mark_acc_s(n)+icond1(n)*icond2(n)
            enddo

            do n = 1, gs%L
            do k = 1, gs%n_bin
              gs%MS(k,n)%dmassdt(imt,j)=gs%MS(k,n)%dmassdt(imt,j)*mod_dum(n)
              acc_mod_s(k,j,n)=acc_mod_s(k,j,n)*mod_dum(n)
            end do
            end do

          else if((j == 11 .or. j == 12 ).and. flagp_s > 0 ) then
            ! - immersion freezing or homogeneous freezing
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if( icond2(n)==1) then
                gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n = 1, gs%L
            do k = 1, gs%n_bin
              if( icond2(n)==1 ) then
                gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_r(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n=1,gs%L
              mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
            enddo
          else if( j == 8 ) then
            ! contact ice nucleation
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if( icond2(n)==1) then
                gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n = 1, gs%L
            do k = 1, gs%n_bin
              if( icond2(n)==1) then
                gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_r(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n=1,gs%L
              mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
            enddo

            do n=1,ag%L
              if( icond2(n)==1) then
                ga(2)%ms(1,n)%dmassdt(amt,j) = ga(2)%ms(1,n)%dmassdt(amt,j) * modc_r(n)
                acc_mod_a(1,j,2,n) = acc_mod_a(1,j,2,n) * modc_r(n)
              endif
              mark_acc_a(2,n)=mark_acc_a(2,n) + icond1(n)*icond2(n)
            enddo
          else
            do n=1,ag%L
              if( icond2(n)==1) then
                gr%MS(i,n)%dmassdt(rmt,j) = gr%MS(i,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(i,j,n) = acc_mod_r(i,j,n) * modc_r(n)
              end if
            enddo
          endif
        end do
      enddo liqbin_loop1
```

The analogous concentration-rescale liqbin loop lives in `cal_con_budget_col` (`liqbin_loop1` at `mod_amps_check.F90:3207-3451`) and the circumscribing-volume rescale in `cal_vol_budget_col` (3925-4274); both follow the identical gain/loss → `modc` → per-process multiply pattern, keyed on concentration and `Vcs` respectively.

---

### File/line index

| Item | File | Lines |
|---|---|---|
| `w_terminal_vel` | `contrib/AMPS/mod_amps_core.F90` | 3258-3407 |
| `cal_wterm_vel_v3_vec` | `contrib/AMPS/class_Group.F90` | 8537-9010 |
| `cal_best_number` | `contrib/AMPS/class_Mass_Bin.F90` | 3279-3310 |
| `sclsedprz_original` | `contrib/AMPS/mod_amps_utility.F90` | 3668-4980 |
| — Euler/upstream default (`else`) | same | 4611-4827 |
| — PPM (`iadvv==12`) | same | 4103-4302 |
| — mid-point (`iadvv==22`) | same | 4304-4609 |
| — CFL `niter` | same | 3967-3990 |
| — tendency accumulation (M2) | same | 4932-4980 |
| `cal_ppmcoef` | `contrib/AMPS/mod_scladv_slppm.F90` | 14-97 |
| `cal_aLRa6_z` | same | 272-421 |
| `cal_flux_z` | same | 423-593 |
| `shift_bin_vec` | `contrib/AMPS/mod_amps_utility.F90` | 6584-7720 |
| `repair` driver | `contrib/AMPS/mod_amps_check.F90` | 22-1181 (limits at 112-113) |
| `cal_mass_budget_col` liqbin_loop1 | same | 1261-1466 |
| `cal_con_budget_col` liqbin_loop1 | same | 3207-3451 |
| `cal_vol_budget_col` | same | 3925-4274 |

Two behavioral notes for M2 replication: (1) `sclsedprz_original` has no `iadvv=1`/`iadvv=2`; the Euler scheme is the fall-through `else` and PPM is `iadvv==12`. (2) `cal_flux_z` is not on `sclsedprz_original`'s path (PPM flux is inlined at 4172-4211); `cal_flux_z`'s only caller is `sclsedaer` at `mod_amps_utility.F90:5217`.