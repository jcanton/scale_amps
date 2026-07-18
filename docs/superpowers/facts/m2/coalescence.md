# AMPS Collision–Coalescence Engine — Verbatim Extraction (RAIN-RAIN warm case + collisional breakup)

Source file: `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_core.F90` (28357 lines) except `cal_needgive` which lives in `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Group.F90`.

All quoted code is copied verbatim with `file:line` anchors. Fortran fixed offsets given as line ranges.

---

## 0. Call graph / M2-relevant summary

`coalescence` (`mod_amps_core.F90:1061`) is the generic two-group collection driver. For the warm rain-rain case `g_1%token==1 .and. g_2%token==1`. Per column `n` it:

1. builds the O(nbins²) collection kernel `KC(i,j,n)` and coalescence efficiency `E_coal(i,j,n)` via `cal_collision_kernel_func` (`:1602`, def `:15818`);
2. accumulates collected number `N_col(i,j,n)` and used amounts in a bin-pair double loop (`:1634`);
3. iteratively fixes over-depletion (`iter_loop1`, `:1715`);
4. in `collector_loop1` (`:1971`) redistributes each collector bin `i`'s mass to sub-bins (fortunate/unfortunate/high categories), computes transferred mass, calls `cal_ratio_mass_col_vec` (`:2483`, def `:17586`), performs the bin-boundary search + linear-remap (`cal_lincubprms_vec`/`cal_linprms_vec_s`, `:2539–2608`) and scatters into target original bins via `cal_transbin_vec` (`:2746`);
5. when `ibreak==1` (micexfg 18, ON for cloudlab) and `pro_type==2` adds Low-List fragments via `add_fragments_col_vec` (`:2890`, def `:15659`), using the pre-tabulated breakup distribution `bu_fd`/`bu_tmass` filled by `cal_breakup_dis_LL` (`:12019`).

The **`ibreak` flag path** threads through: passed into `coalescence` (`:1098`) → `cal_collision_kernel_func` (`:1602`,`:15818`) where at `:16341` it sets `E_coal(i,j)=1.0` for pairs with no tabulated breakup mass → the `ibreak==1` blocks in the `N_col` accounting loops (`:1668`, `:1754`, `:1839`) that add breakup-consumed number `used_N_b` → `add_fragments_col_vec` at `:2890`. Breakup number is always `N_bk = N_col(i,j,n)*max(0, 1 - E_coal(i,j,n))`.

`add_simple_vec` (`:15452`) and `add_samebin_vec` (`:15579`) are the alternate scatter helpers (bin-boundary search + scatter, and same-bin accumulate respectively); in this code path the main scatter is done by `cal_transbin_vec`, but `add_simple_vec` shows the canonical `Mp/Np → jbin` boundary search that M2 must emit.

---

## 1. `coalescence` — generic two-group collection engine

### 1.1 Argument list + key declarations (`:1061–1400`)

Header and dummy arguments (`:1061–1116`):

```fortran
  subroutine coalescence(g_1, g_2, ag, level,col_level,mes_rc,ibreak &
              ,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd,ID,JD,KD &
              ,adrpdrp,drpdrp &
              ,ahexdrp,hexdrp &
              ,abbcdrp,bbcdrp &
              ,acoldrp,coldrp &
              ,agp1drp,gp1drp &
              ,agp4drp,gp4drp &
              ,agp8drp,gp8drp &
! <<< 2014/10 T. Hashino added for KiD
              ,dM_auto,dM_accr)
    use scale_prc, only: &
       PRC_abort
    use class_Group, only: &
       col_lut_aux
    use mod_amps_utility, only: &
       cal_lincubprms_vec, &
       cal_linprms_vec_s, &
       cal_transbin_vec
    type (Group), intent(inout)   :: g_1, g_2
    type (AirGroup), intent(in) :: ag
    integer, intent(in)           :: level, col_level
    integer,intent(in),dimension(*)  :: mes_rc
    ! flag for collisional breakup calculation
    ! 1: implement calculation.
    integer,intent(in)  :: ibreak
    integer,intent(in) :: imin_bk,imax_bk,jmin_bk,jmax_bk
    real(PS),intent(in) :: bu_fd(2,*),bu_tmass(*)  ! 2014/10 T. Hashino modify for KID
    !
    type(col_lut_aux),intent(in) :: adrpdrp
    real(PS),intent(in),dimension(adrpdrp%nr,*) :: drpdrp
    type(col_lut_aux),intent(in) :: ahexdrp,abbcdrp,acoldrp
    real(PS),intent(in),dimension(ahexdrp%nr,*) :: hexdrp
    real(PS),intent(in),dimension(abbcdrp%nr,*) :: bbcdrp
    real(PS),intent(in),dimension(acoldrp%nr,*) :: coldrp
    type(col_lut_aux),intent(in) :: agp1drp,agp4drp,agp8drp
    real(PS),intent(in),dimension(agp1drp%nr,*) :: gp1drp
    real(PS),intent(in),dimension(agp4drp%nr,*) :: gp4drp
    real(PS),intent(in),dimension(agp8drp%nr,*) :: gp8drp
    !
    integer :: ID(*),JD(*),KD(*)
    real(PS),intent(inout),dimension(*) :: dM_auto,dM_accr
```

The `ibreak` documentation is at `:1096–1098`: "flag for collisional breakup calculation / 1: implement calculation."

The 7 collision-efficiency lookup tables arrive as pairs `(a<name>, <name>)` where `a<name>` is a `col_lut_aux` descriptor (see §3): `drpdrp` (drop-drop), `hexdrp`, `bbcdrp`, `coldrp`, `gp1drp`, `gp4drp`, `gp8drp`.

The central O(nbins²) work arrays (all `mxnbin`-dimensioned, per column `n≤LMAX`):

```fortran
    ! collection Kernel
    real(PS), dimension(mxnbin,mxnbin,LMAX)          :: KC          ! :1214
    ! coalescence efficiency
    real(PS), dimension(mxnbin,mxnbin,LMAX)          :: E_coal      ! :1218
    ! number of collitional breakup
    real(PS)      :: N_bk                                           ! :1227
```

```fortran
    ! total number of collected hydrometeors in j bin by i bin, N_col(i,j)
    real(PS),dimension(mxnbin,mxnbin,LMAX)   ::  N_col              ! :1362
    ! ratio of collection to total concentration of collector, N_col/N_i
    real(PS),dimension(mxnbin,LMAX)   ::  col_ratio                 ! :1366
    ! left over
    real(PS), dimension(mxnbin,LMAX)   ::  left_N, left_M           ! :1371
    ! total mass that each sub bin collects from jth bin
    real(PS), dimension(mxnbin+1,mxnbin,LMAX)   ::  dM_col          ! :1376
    ! number of drops that each sub bin collects from jth bin (per parent)
    real(PS), dimension(mxnbin+1,mxnbin,LMAX)   ::  N_drops         ! :1381
    ! number of drops that sub bin i collects from all bins (per parent)
    real(PS), dimension(mxnbin+1,LMAX)   ::  dN_ice                 ! :1386
    ! checker of collision, 0: no collision, 1: collide
    integer, dimension(mxnbin+1,mxnbin,LMAX)       :: checker       ! :1391
```

Parameters (`:1397–1400`):

```fortran
    ! bin boundary modification factor.
    real(PS),parameter    :: bbmf=0.2D0
    ! minimum collector mass to be considered.
    real(PS),parameter :: amin_mass=4.188790d-12
```

`used_N_2`, `used_M_2` (double precision `:1248–1249`), `used_N_b` (`:1262`), and the new-bin accumulators `new_N_1(mxnbin,LMAX)` (`:1129`), `new_M_1(mxnbin,g_1%L,1+mxnmasscomp)` (`:1130`), `new_Q_1` (`:1175`) are the outputs.

### 1.2 Kernel build + the O(nbins²) bin-pair double loop (`:1581–1664`)

The pair activation mask `icycle_ijn(i,j,n)` (0 = active) is built first (`:1581–1595`):

```fortran
       do n = 1, g_1%L
       do i = 1, g_1%N_BIN
       do j = 1, g_2%N_BIN
         if( g_1%MS(i,n)%con > 1.0e-30_PS .and. &
             g_1%MS(i,n)%mass(1) > 1.0e-30_PS .and. &
             g_1%MS(i,n)%mean_mass > 1.0e-15_PS .and. &
             g_2%MS(j,n)%con > 1.0e-30_PS .and. &
             g_2%MS(j,n)%mass(1) > 1.0e-30_PS .and. &
             g_2%MS(j,n)%mean_mass > 1.0e-15_PS ) then
           icycle_ijn(i,j,n)=0
         endif
       enddo
       enddo
       enddo
```

Kernel evaluation, once per column (`:1600–1621`):

```fortran
       do n=1,g_1%L
          if(icycle_n(n)==0) then
            call cal_collision_kernel_func(g_1,g_2,ag%TV(n),col_level,ibreak,&
               imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,&
               n,KC(:,:,n),E_coal(:,:,n), &
               adrpdrp,drpdrp,&
               ahexdrp,hexdrp,&
               abbcdrp,bbcdrp,&
               acoldrp,coldrp,&
               agp1drp,gp1drp,&
               agp4drp,gp4drp,&
               agp8drp,gp8drp)
          endif
       enddo
```

**The core O(nbins²) collection kernel** — this is the loop nest M2 must emit as SPLIT operators (`:1634–1664`). Note the collector index `i` runs *descending* `g_1%N_BIN → 1`:

```fortran
       do i=g_1%N_BIN, 1, -1
          do n = 1, g_2%L
          do j = 1, g_2%N_BIN

           if(icycle_ijn(i,j,n)==0) then

             ! +++ prevent errors from numerical diffusion +++
             if( g_1%MS(i,n)%con < 1.0e-30_PS .or. g_2%MS(j,n)%con < 1.0e-30_PS ) then
                N_col(i,j,n) = 0.0_PS
             else if( KC(i,j,n) > 0.0_PS ) then
                ! --- only consider the case of i catching j +++
                N_col(i,j,n) =  g_1%MS(i,n)%con * KC(i,j,n)
             else if( KC(i,j,n) < 0.0_PS .and. g_1%token /= g_2%token ) then
                ! --- case of riming ...
                N_col(i,j,n) =  g_1%MS(i,n)%con * abs(KC(i,j,n))
             end if

             if( N_col(i,j,n) > g_2%MS(j,n)%con ) then
                N_col(i,j,n) = g_2%MS(j,n)%con
             end if
             used_N_2(j,n)=used_N_2(j,n)+N_col(i,j,n)*E_coal(i,j,n)
             used_M_2(j,n)=used_M_2(j,n)+N_col(i,j,n)*E_coal(i,j,n)*g_2%MS(j,n)%mean_mass
           endif
         enddo
         enddo
       enddo
```

For warm rain (`g_1%token==g_2%token==1`), only `KC>0` is used; the `abs(KC)` riming branch is inactive.

### 1.3 `ibreak` breakup-consumption accounting (`:1668–1707`)

```fortran
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

### 1.4 Over-depletion fix `iter_loop1` (`:1713–1890`)

The full iteration is long; the essential structure per collectee bin `i` (`:1713–1799` number-limit half, and `:1801–1884` mass-limit half). Number-limit half:

```fortran
       do i=1,g_2%N_BIN
         iter_loop1: do ii=1,iter
           mod_ratio(1:g_2%L)=1.0_PS
           icond2_n(1:g_2%L)=0
           do n=1,g_2%L
             if( g_2%MS(i,n)%con<used_N_2(i,n) ) then
               mod_ratio(n)=g_2%MS(i,n)%con/used_N_2(i,n)
               icond2_n(n)=1
             endif
           enddo
           idone_n=.true.
           if(any(icond2_n(1:g_2%L)>0)) then
             idone_n=.false.
             do n = 1, g_1%L
             do j = 1, g_1%N_BIN
               N_col(j,i,n)=N_col(j,i,n)*min(1.0_PS,mod_ratio(n))
             enddo
             enddo
             ...
             do n=1,g_2%L
               if(icond2_n(n)==1) then
                 used_N_2(i,n)=g_2%MS(i,n)%con
                 used_marker(i,n)=1
               endif
             enddo
             if(ibreak==1) then
                ... rescale N_col(i,j,n) for j/=i, recompute used_N_b ...
             endif
           endif
```

The mass-limit half (`:1801–1884`) is the mirror using `g_2%MS(i,n)%mass(1)` and `used_M_2`, setting `used_M_2(i,n)=g_2%MS(i,n)%mass(1)` and `used_marker(i,n)=1`. Convergence at `:1886`:

```fortran
           if(idone_n.and.idone_m) then
             exit
           endif
         enddo iter_loop1
       enddo
```

### 1.5 Left-over computation (`:1913–1933`)

For like-token (rain-rain) collisions the residual (uncollided) part is retained:

```fortran
       if( g_1%token == g_2%token ) then
          do n = 1, g_1%L
          do i = 1, g_1%N_BIN
           left_N(i,n)=max(g_2%MS(i,n)%con-used_N_2(i,n),0.0_DS)
           left_M(i,n)=max(g_2%MS(i,n)%mass(1)-used_M_2(i,n),0.0_DS)
         end do
         end do
       else
          do n = 1, g_1%L
          do i = 1, g_1%N_BIN
           left_N(i,n) = g_1%MS(i,n)%con
           left_M(i,n) = g_1%MS(i,n)%mass(1)
         end do
         end do
       end if
```

### 1.6 `collector_loop1` — scatter-to-target-bin (the hardest M2 kernel)

The collector loop runs descending, `i = g_1%N_bin → icolbin_min` (`:1971`), where `icolbin_min` is the smallest collector bin whose upper boundary reaches `amin_mass` (`:1518–1523`):

```fortran
       icolbin_min=1
       do i=1,g_1%N_bin
         if(g_1%binb(i+1)<amin_mass) then
           icolbin_min=i
         endif
       enddo
```

```fortran
       collector_loop1: do i=g_1%N_bin,icolbin_min,-1
```

Per-collector work arrays are zeroed (`:1982–2016`).

**Categorization of collectee bins into collection-rate groups** (`:2028–2054`). `col_ratio(j,n)=E_coal*N_col/N_collector`; >10 or [1,10] → "high" list `jbin_h`; (0,1) or [1,10] → "low/medium" list `jbin_lm`:

```fortran
         do j= jmax, 1, -1
           do n=1,g_1%L
             if(icond1(i,n)==1) then
               col_ratio(j,n) = E_coal(i,j,n)*N_col(i,j,n)/g_1%MS(i,n)%con
               if( N_col(i,j,n) > 0.0_PS ) then
                 if( col_ratio(j,n) < 0.0_PS ) then
                 else if( col_ratio(j,n) > 10.0_PS ) then
                   n_h(n) = n_h(n) + 1
                   jbin_h(n_h(n),n) = j
                 else if( col_ratio(j,n) > 0.0_PS .and. col_ratio(j,n) < 1.0_PS ) then
                   n_lm(n) = n_lm(n) + 1
                   jbin_lm(n_lm(n),n) = j
                 else if( col_ratio(j,n) >= 1.0_PS .and. col_ratio(j,n) <= 10.0_PS ) then
                   n_h(n) = n_h(n) + 1
                   jbin_h(n_h(n),n) = j
                   n_lm(n) = n_lm(n) + 1
                   jbin_lm(n_lm(n),n) = j
                 endif
               end if
             end if
           end do
         end do
         n_lm_max=maxval(n_lm(1:g_1%L))
         n_h_max=maxval(n_h(1:g_1%L))
```

**Fortunate / unfortunate fraction determination** (`:2060–2087`): a "fortunate" fraction `p_F` (particles that collide) per low/medium collectee, accumulated in `sum_ndrop_B`, seeds sub-bins `Np(jj,n)`, `Mp(jj,n)`:

```fortran
         sum_ndrop_B(1:g_1%L) = 0.0_PS
         icond2_n(1:g_1%L)=0
         do jj= 1, n_lm_max
           do n=1,g_1%L
             if(icond1(i,n)==1.and.jj<=n_lm(n).and.icond2_n(n)==0) then
               j = jbin_lm(jj,n)
               if( col_ratio(j,n) > 0.0_PS .and. col_ratio(j,n) < 1.0_PS ) then
                 p_F = col_ratio(j,n)
                 sum_ndrop_B(n) = sum_ndrop_B(n) + p_F
               else if( col_ratio(j,n) >= 1.0_PS .and. col_ratio(j,n) <= 10.0_PS ) then
                 p_F = col_ratio(j,n) - aint(col_ratio(j,n))
                 sum_ndrop_B(n) = sum_ndrop_B(n) + p_F
               else
                 p_F = -1.0_PS
               end if
               if(p_F>0.0_PS) then
                 if( sum_ndrop_B(n) >= 1.0_PS ) then
                   sum_ndrop_B(n) = sum_ndrop_B(n) - p_F ! subtract back as sum (R) > 1
                   icond2_n(n)=1 ! these g2 bins are excluded as sum (R) > 1
                 else
                   Np(jj,n)=p_F*left_N(i,n) ! fraction of g1 particles that will collide
                   Mp(jj,n) = Np(jj,n)*g_1%MS(i,n)%mean_mass ! total original mass
                   n_f(n) = n_f(n) + 1
                 end if
               end if
             end if
           end do
         end do
```

The unfortunate group (remaining `1-sum_ndrop_B` fraction) is placed in the last sub-bin `n_all(n)` (`:2089` onward): `Np(n_all(n),n)=(1-sum_ndrop_B(n))*left_N(i,n)`; `Mp(...)=Np*mean_mass`.

Then the mass-component seed `Mc(j,k,n)` is built from `g_1%MS(i,n)%mass` ratios (`:2185`-region, offsets `189–239` of the earlier map): `Mc(j,k,n) = Mp(j,n)*(g_1%MS(i,n)%mass(1+k)/g_1%MS(i,n)%mass(1))`.

The **collection-of-mass accumulation** into `dM_col(k,j,n)` and per-parent counts `N_drops(k,j,n)` happens in three passes (high `jbin_h`, fortunate, low/medium), e.g. the "H" pass (map offsets 259–306):

```fortran
               Mp(k,n) = Mp(k,n) + dM_dum
               dM_col(k,j,n) = dM_col(k,j,n) + dM_dum
               N_drops(k,j,n) = N_drops(k,j,n) + ndrop_B
```

and the fortunate pass (`:2277`-region) with `dM_dum = N_col(i,j,n)*E_coal(i,j,n)*g_2%MS(j,n)%mean_mass`.

**No-collision rollback** using `checker` (map offset 480–502, i.e. ~`:2434`):

```fortran
            do n = 1, g_2%L
            do k = 1, g_2%N_BIN
            do j = 1, n_all_max
               if( N_drops(j,k,n) > 0.0_PS .and. checker(j,k,n) == 0) then
                 dM_dum = N_drops(j,k,n)*&
                       g_1%MS(i,n)%con*g_2%MS(k,n)%mean_mass*(Np(j,n)/max(left_N(i,n),1e-30_PS))
                 Mp(j,n)=Mp(j,n)-dM_dum
                 dmass_max(j,n)=dmass_max(j,n)-N_drops(j,k,n)*g_2%binb(k+1)
                 dmass_min(j,n)=dmass_min(j,n)-N_drops(j,k,n)*g_2%binb(k)
                 used_N_2(k,n)=used_N_2(k,n)-N_drops(j,k,n)*g_1%MS(i,n)%con
               endif
```

**Mass-component ratios** `ratio_Mp` via `cal_ratio_mass_col_vec` (`:2483`):

```fortran
         ratio_Mp(:,:,:) = 0.0_PS
         call cal_ratio_mass_col_vec( g_1, g_2, ag, level, icond3, &
                 n_all, Mc, dN_ice, ratio_M_2, dM_col, Mp, ratio_Mp)
```

**Bin-boundary construction + linear-remap** (`:2537–2608`). The shifted bin boundaries are the original collector boundaries offset by the accumulated `dmass_min/dmass_max`, then a cubic/linear-parameter fit is done; on error the boundaries are widened by `bbmf` and a linear fit is retried:

```fortran
         do n = 1, g_1%L
         do j = 1, n_all_max
           binb3d(j,n,1)=g_1%binb(i)  +dmass_min(j,n)
           binb3d(j,n,2)=g_1%binb(i+1)+dmass_max(j,n)
           Npd(j,n)=Np(j,n)
           Mpd(j,n)=Mp(j,n)
         enddo
         enddo
         call cal_lincubprms_vec(mxnbin+1,n_all_max,g_1%L,Npd,Mpd,binb3d  &
             ,a2d,error_number,"coal_1")

         do n = 1, g_1%L
         do j = 1, n_all_max
           if(1<=error_number(j,n).and.error_number(j,n)<=4) then
             dum1 = Mpd(j,n)/Npd(j,n)
             i_d_ge_b=int(0.5_PS*(1.0_PS+sign(1.0_PS,dum1-binb3d(j,n,2))))
             binb3d(j,n,2)=real(i_d_ge_b,PS_KIND)*dum1*(1.0+bbmf) + &
                           (1.0-real(i_d_ge_b,PS_KIND))*binb3d(j,n,2)
             i_b_ge_d=int(0.5_PS*(1.0_PS+sign(1.0_PS,binb3d(j,n,1)-dum1)))
             binb3d(j,n,1)=real(i_b_ge_d,PS_KIND)*max(g_1%binb(1),dum1*(1.0-bbmf)) + &
                           (1.0-real(i_b_ge_d,PS_KIND))*binb3d(j,n,1)
             a(1) = a2d(j,n,1); a(2) = a2d(j,n,2); a(3) = a2d(j,n,3); a(4) = a2d(j,n,4)
             call cal_linprms_vec_s(Npd(j,n),Mpd(j,n), &
                      binb3d(j,n,1),binb3d(j,n,2),a(:), &
                      error_number(j,n))
           endif
         enddo
         enddo
```

(The `i_d_ge_b`/`i_b_ge_d` idiom is a branchless `sign()` clamp: `int(0.5*(1+sign(1,x)))` is 1 if `x>=0` else 0.)

**Scatter into target original bins** via `cal_transbin_vec` (`:2746–2758`):

```fortran
         call cal_transbin_vec(g_1%token, &
                           g_1%L,g_1%N_masscom, &
                           g_1%N_bin,n_all_max, &
                           g_1%binb, &
                           error_number, &
                           a2d,binb3d,mtend, &
                           new_N_1,new_M_1,new_Q_1, &
                           new_mtend, &
                           ratio_Mp,den_ip_p,axr_p,spx_p, &
                           habit_p,den_ic_p, &
                           rag_p,rcg_p,n_exice_p, &
                           actINF_p, &
                           0)
```

For warm rain there is a second `cal_transbin_vec` call for the autoconversion diagnostic (`:2847`, guarded by `if(g_1%token==1.and.g_2%token==1)` at `:2770`), producing `dM_auto`/`dM_accr` (`:2861–2874`).

**Collisional-breakup fragment scatter** at the end of the collector loop (`:2888–2895`):

```fortran
         ! calculation of concentration and mass transfer
         ! by collisional breakup
         if(pro_type==2.and.ibreak==1) then
           call add_fragments_col_vec(new_N_1,new_M_1, &
                    i,g_1,g_2,icond1,used_marker,E_coal,N_col, &
                    imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd)
         endif
       enddo collector_loop1
```

**Post-loop left-over re-add** for like-token (`:2899–2925`):

```fortran
       if( g_1%token == g_2%token ) then
          do n = 1, g_1%L
          do i = 1, g_1%N_BIN
             if( used_marker(i,n) == 1 ) then
                new_N_1(i,n) = new_N_1(i,n) + left_N(i,n)
                new_M_1(i,n,1) = new_M_1(i,n,1) + left_M(i,n)
             endif
          enddo
          enddo
          do j=1,g_1%N_masscom
             do n = 1, g_1%L
             do i = 1, g_1%N_BIN
                if( used_marker(i,n) == 1 ) then
                   new_M_1(i,n,1+j) = new_M_1(i,n,1+j) &
                                    + left_M(i,n)*ratio_M_2(i,j,n)
                endif
             enddo
             enddo
          enddo
       end if
```

Followed by `call check_csvolume_vec(g_1,new_N_1,new_M_1,new_Q_1)` (`:2931`).

---

## 2. Collision-kernel evaluation `cal_collision_kernel_func` (`:15818–16361`)

Header (`:15818–15827`) and the warm-rain-relevant declarations:

```fortran
  subroutine cal_collision_kernel_func(g_1,g_2,th_var,col_level,ibreak,&
               imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,&
               n,KC,E_coal,&
               adrpdrp,drpdrp,&
               ahexdrp,hexdrp,&
               abbcdrp,bbcdrp,&
               acoldrp,coldrp,&
               agp1drp,gp1drp,&
               agp4drp,gp4drp,&
               agp8drp,gp8drp)
```

```fortran
    real (PS), intent(inout), dimension(mxnbin,*)     :: KC,E_coal    ! :15852
    real(PS),dimension(g_1%n_bin,g_2%n_bin)  :: E_c                   ! collision efficiency :15868
    real(PS),dimension(g_1%n_bin,g_2%n_bin)  :: A_c                   ! sweep cross-section  :15870
    real(PS) :: Nre_rp,Nre_sp                                          ! :15884
    integer :: i1,j1                                                   ! :15885
    real(PS) :: x1,y1,wx,wy                                            ! :15886
    real(PS),parameter :: ec_min=0.0D0                                ! :15887
    real(PS) :: NreL_p,rrat,rrat_p                                    ! :15888
```

Initialization + per-pair activation mask `icond1(i,j)` (`:15894–15907`):

```fortran
    do j = 1, g_2%N_BIN
    do i = 1, g_1%N_BIN
     KC(i,j)=0.0_PS
      E_c(i,j)=0.0_PS
      E_coal(i,j)=0.0_PS
      em(i,j)=0
      if((g_1%MS(i,n)%con>1.0e-30_PS .and. g_2%MS(j,n)%con>1.0e-30_PS ).and.&
         (g_1%MS(i,n)%mass(1)>1.0e-30_PS .and. g_2%MS(j,n)%mass(1)>1.0e-30_PS)) then
        icond1(i,j)=0
      else
        icond1(i,j)=1
      endif
    enddo
    enddo
```

**Warm-rain sweep cross-section** `A_c` (`:15913–15942`, token 1-1):

```fortran
    if(g_1%token==1.and.g_2%token==1 ) then
      if (col_level==0) then
         do j = 1, g_2%N_BIN
         do i = 1, g_1%N_BIN
          if(icond1(i,j)==0) then
            if(g_2%MS(j,n)%len > 0.0001_PS .or. g_1%MS(i,n)%len > 0.0001_PS) then
              A_c(i,j)=0.0_PS
            else
              A_c(i,j)=0.25_PS*PI*(g_1%MS(i,n)%len+g_2%MS(j,n)%len)**2.0
            endif
          endif
        enddo
        enddo
      else
         do j = 1, g_2%N_BIN
         do i = 1, g_1%N_BIN
          if(icond1(i,j)==0) then
            A_c(i,j)=0.25_PS*PI*(g_1%MS(i,n)%len+g_2%MS(j,n)%len)**2.0
          endif
        enddo
        enddo
      endif
```

Warm-rain **collision efficiency** `E_c` from the `drpdrp` table — see §3 for the verbatim bilinear gather (`:16061–16107`).

**Coalescence efficiency** (`:16277–16285`):

```fortran
    do j = 1, g_2%N_BIN
    do i = 1, g_1%N_BIN
      if(icond1(i,j)==0) then
        call cal_Coalescence_Efficiency(g_1,i,g_2,j,n,th_var,E_coal(i,j),&
             D_L,D_S,S_T,S_C,DS_S,CKE)
      endif
    enddo
    enddo
```

**Kernel assembly** `KC = E_c·(v_i−v_j)·A_c·N_j·dt` (`:16294–16310`):

```fortran
    do j = 1, g_2%N_BIN
    do i = 1, g_1%N_BIN
      if(icond1(i,j)==0) then
        KC(i,j)=E_c(i,j)*( g_1%MS(i,n)%vtm - g_2%MS(j,n)%vtm )*A_c(i,j)* &
                g_2%MS(j,n)%con*g_1%dt
      endif
    enddo
    enddo
```

(Commented alternatives for Golovin / constant kernels are present at `:16299–16307`.)

**`ibreak` path inside the kernel** (`:16341–16359`) — pairs outside the breakup mass table, or with no tabulated breakup mass, are forced `E_coal=1` so no fragments are generated:

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

### 2.1 `cal_Coalescence_Efficiency` (`:11796–12017`) — warm-rain branch

Header + Low-List (1982) constants:

```fortran
  subroutine cal_Coalescence_Efficiency(g_1,i,g_2,j,n,th_var,E_coal,&
       D_L,D_S,S_T,S_C,DS_S,CKE)
    type (Group), intent(in)               :: g_1, g_2
    integer, intent(in)                    :: i, j,n
    type (Thermo_Var) :: th_var
    real(PS)                    :: E_coal
    real(PS)                    :: E_int,E_stick
    real(PS)                    :: D_L, D_S, V_L, V_S
    real(PS)                    :: E_T, CKE
    real(PS)                    :: S_T, S_C, dS_S
    real(PS), parameter         :: den_w = 1000.0D0
    real(PS), parameter         :: a = 0.778D0
    real(PS), parameter         :: b = 2.61d+6
    real(PS),parameter :: D_0=0.01D0
```

The liquid-liquid (token 1-1) case — the warm-rain path (`:11842–11880`):

```fortran
    if( g_1%token == 1 .and. g_2%token == 1 ) then
       ! define large diameter and small diameter (m)
       D_L = max( g_1%MS(i,n)%len, g_2%MS(j,n)%len ) * 1.0e-2
       D_S = min( g_1%MS(i,n)%len, g_2%MS(j,n)%len ) * 1.0e-2

       if(min(D_L,D_S)*100.0<D_0) then
         E_coal=1.0_PS
       else
         V_L = max( g_1%MS(i,n)%vtm, g_2%MS(j,n)%vtm ) * 1.0e-2
         V_S = min( g_1%MS(i,n)%vtm, g_2%MS(j,n)%vtm ) * 1.0e-2

         S_T = PI*th_var%sig_wa*1.0e-3_DS*( D_L**2.0 + D_S**2.0)
         S_C = PI*th_var%sig_wa*1.0e-3_DS*( D_L**3.0 + D_S**3.0)**(2.0/3.0)
         dS_S = S_T - S_C
         CKE = (den_w*PI/12.0_PS)*((V_L-V_S)**2.0)*&
              ((D_L*D_S)**3.0)/(D_L**3.0+D_S**3.0)
         E_T = CKE + dS_S

         if( E_T < 5.0e-6_PS ) then
            E_coal = a * ((1.0_PS + D_S/D_L)**(-2.0)) * &
               exp( - b*th_var%sig_wa*1.0e-3_DS*(E_T**2.0)/S_c)
         else
            E_coal = 0.0_PS
         end if
       end if
```

(The remaining branches handle liquid-ice / ice-ice aggregation with `stick_lmt`, `E_int`, `E_stick`; not on the warm path.)

### 2.2 Low-List breakup distribution `cal_breakup_dis_LL` (`:12019–12447`)

This is the routine that pre-fills `bu_fd`/`bu_tmass` (micexfg 18). Full verbatim; header + argument declarations:

```fortran
  subroutine cal_breakup_dis_LL(g_1,i,j,imin_bk,jmin_bk,bu_tmass,bu_fd,&
                                xD_L,xD_S,xS_T,xS_C,xCKE)
    use mod_amps_utility, only: &
       getznorm2
    type (Group), intent(in)               :: g_1
    integer,intent(in)   :: i,j
    integer,intent(in) :: imin_bk,jmin_bk
    real(PS),intent(inout) :: bu_fd(2,*),bu_tmass(*)
    real(PS)                    :: xD_L,xD_S
    real(PS),intent(in)                    :: xCKE
    real(PS),intent(in)         :: xS_T,xS_C
    real(DS)                    :: D_L,D_S
    real(PS) :: S_C
    real(DS)                    :: W1,W2
    real(DS)                   :: CKE
    real(DS)         :: S_T
    real(DS),parameter          :: CKE0=8.93d-7
    real(DS),parameter          :: W0=0.86D0
    real(DS)                    :: R_f,R_s,R_d
```

CGS conversion + coalesced-drop size (`:12102–12118`):

```fortran
    sq_twod=sqrt(2.0_DS)
    D_L=xD_L; D_S=xD_S; S_C=xS_C; CKE=xCKE; S_T=xS_T; S_C=xS_C
    D_L=D_L*100.0_DS
    D_S=D_S*100.0_DS
    D_coal=(D_L**3.0+D_S**3.0)**(1.0/3.0)
    m_coal=coedpi6*D_coal**3.0
    if(D_coal<=D_0) return
```

Weber numbers + breakup-type fractions (filament R_f, sheet R_s, disk R_d) (`:12122–12145`):

```fortran
    W1=CKE/S_C
    W2=CKE/S_T
    if(CKE>=CKE0) then
       R_f=1.11e-4_DS*CKE**(-0.654)
    else
       R_f=1.0_DS
    end if
    if(W2>=W0) then
       R_s=0.685_DS*(1.0_DS-dexp(-1.63_DS*(W2-W0)))
    else
       R_s=0.0_DS
    end if
    if(R_s+R_f>1.0_DS) then
       R_s=1.0_DS-R_f
       R_d=0.0_DS
    else
       R_d=dmax1(1.0_DS-R_f-R_s,0.0_DS)
    end if
```

Mean fragment counts F_f/F_s/F_d (`:12150–12165`):

```fortran
    F_f=(-2.25e+4_DS*(D_L-0.403_DS)*(D_L-0.403_DS)-37.9_DS)*D_S**2.5+&
         9.67_DS*(D_L+0.170_DS)*(D_L+0.170_DS)+4.95_DS
    F_f=dmax1(2.0_DS,dmin1(F_f,app*D_S**bpp+2.0_DS))
    F_s=dmax1(5.0_DS*&
         (2.0_DS*getznorm2(sq_twod*(S_T-2.53e-6_DS)/1.85e-6_DS)-1.0_DS)&
         +6.0_DS,2.0_DS)
    F_d=dmax1(297.5_DS+23.76_DS*dlog(CKE),2.0_DS)
```

Parent-distribution parameters (`:12174–12222`) — heights `H_*`, means `mu_*`, and `cal_sig_sf`/`cal_Hmusig` calls for `sig_*`. Then the per-bin fragment number `n_f/n_s/n_d` (lognormal fragments + normal parents, `:12256–12308`) and fragment mass `m_f/m_s/m_d` (`:12313–12400`) are integrated over each liquid bin `[D_1,D_2]` using `getznorm2` (the normal CDF).

Assembly of total `dcon`/`dmass` per bin and mean-mass clamp (`:12404–12429`):

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

**The triangular pair-index and normalized storage into `bu_fd`/`bu_tmass`** (`:12432–12443`) — this is exactly the indexing `add_fragments_col_vec` and the kernel `ibreak` block read back:

```fortran
    i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2
    bu_tmass(i1d_pair)=m_coal

    mrat=m_coal/sum(dmass)
    Do ibin=1,g_1%N_BIN
       kk=(i1d_pair-1)*g_1%N_BIN+ibin
       bu_fd(2,kk)=mrat*dcon(ibin)
       bu_fd(1,kk)=mrat*dmass(ibin)
    end do
```

`bu_fd(2,·)` = fragment number distribution, `bu_fd(1,·)` = fragment mass distribution, mass-normalized so `sum(bu_fd(1,·))=m_coal`.

### 2.3 `add_fragments_col_vec` (`:15659–15774`) — full

```fortran
  subroutine add_fragments_col_vec(new_N_1,new_M_1, &
                    i,g_1,g_2,icond1,used_marker,E_coal,N_col, &
                    imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd)

    type (Group), intent(in)   :: g_1, g_2
    real(8), dimension(mxnbin,*),intent(inout)                :: new_N_1
    real(8), dimension(mxnbin,g_1%L,1+mxnmasscomp),intent(inout)  :: new_M_1
    integer,intent(in) :: i
    integer,dimension(mxnbin,*),intent(in) :: icond1
    integer,dimension(mxnbin,*),intent(in)      :: used_marker
    real(PS), dimension(mxnbin,mxnbin,*),intent(in)          :: E_coal
    real(PS),dimension(mxnbin,mxnbin,*),intent(in)   ::  N_col
    integer,intent(in) :: imin_bk,imax_bk,jmin_bk,jmax_bk
    real(PS),intent(in) :: bu_fd(2,*),bu_tmass(*)
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

Note the identical triangular index `i1d_pair=j-jmin_bk+1+(i-imin_bk)*(1+i-imin_bk)/2` and `kk=(i1d_pair-1)*g_1%N_BIN+k`; only pairs with `j<i` and both within `[imin_bk,imax_bk]×[jmin_bk,jmax_bk]` contribute. Fragments are added straight into `new_N_1`/`new_M_1` at the fragment bin `k` (no further remap).

### 2.4 `P_breakup` (`:19873–19904`) and `Q_breakup2` (`:19953–19997`) — full

```fortran
  function P_breakup(phase, a_star, max_dim) result(out)
    integer, intent(in)    :: phase           ! 1: water, 2: ice
    real(PS), intent(in)  :: a_star           ! radius of parent (mm)
    real(PS), intent(in)             :: max_dim ! max dimension (radius) rain/aggregate (cm)
    real(PS)             :: k, out
    if( phase == 1 ) then
       if( a_star >= max_dim ) then
          out = 1.0_PS
       else
          out = min( 2.94e-7_PS*exp(3.4_PS*a_star*10.0_PS), 1.0_PS )
       end if
    else if( phase == 2 ) then
       if( a_star >= max_dim ) then
          out = 1.0_PS
       else
          k = 15.03968607_PS/(max_dim*10.0_PS)
          out = min(2.94e-7_PS*exp(k*a_star*10.0_PS),1.0_PS)
       end if
    end if
  end function P_BREAKUP
```

```fortran
  function Q_breakup2(phase, a_star, a, m, AA, BB,switch) result(out)
    use mod_amps_utility, only: &
       fGM
    integer, intent(in)    :: phase           ! 1: water, 2: ice
    real(PS),intent(in)  :: a_star, a, m,AA,BB
    integer, intent(in)  :: switch
    real(PS)             :: out, x1, x2
    if( phase == 1 ) then
       if( switch == 1 ) then
          out = (AA*BB/3.0_PS/m)*(a/a_star)*exp(-BB*a/a_star)
       else if( switch == 2 ) then
          out = -AA*( exp(-BB*m/a_star) - exp(-BB*a/a_star))
       else if( switch == 3 ) then
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

(The older `Q_breakup` with hard-coded `AA=62.3, BB=7.0` is fully commented out at `:19906–19951`.)

---

## 3. Collision-efficiency lookup-table indexing — bilinear gather

The 7 tables are all indexed the same way inside `cal_collision_kernel_func`: a `col_lut_aux` descriptor `a<name>` carries grid metadata `%nr` (rows), `%nc` (cols), `%xs`/`%ys` (axis origins), `%dx`/`%dy` (axis steps); the table `<name>(i1,j1)` is a 2-D array. The runtime indices `i1` (Reynolds axis, `log10(Nre)`) and `j1` (radius-ratio or Reynolds-ratio axis) are computed with clamping, then a bilinear interpolation with weights `wx,wy` (not extrapolated).

**Warm rain-rain `drpdrp` gather** (`:16069–16104`), `x`=radius ratio `rrat_p` (folded to ≤1), `y`=`log10(Nre_L)`:

```fortran
          NreL_p=max(g_1%MS(i,n)%Nre,1.0e-10_RP)
          rrat=g_2%MS(j,n)%len/g_1%MS(i,n)%len
          if(rrat<0.0_RP) then
            em(i,j)=1
            rrat_p=0.0_RP
            NreL_p=max(g_1%MS(i,n)%Nre,1.0e-10_RP)
            rrat=g_2%MS(j,n)%len/g_1%MS(i,n)%len
            ...
          elseif(rrat>1.0_RP) then
            rrat_p=1.0/rrat
          else
            rrat_p=rrat
          endif

          j1=max(1,min(adrpdrp%nc-1 &
               ,int((rrat_p-adrpdrp%xs)/adrpdrp%dx)+1))
          x1=real(j1-1,PS_KIND)*adrpdrp%dx+adrpdrp%xs

          i1=max(1,min(adrpdrp%nr-1 &
               ,int((log10(NreL_p)-adrpdrp%ys)/adrpdrp%dy)+1))
          y1=real(i1-1,PS_KIND)*adrpdrp%dy+adrpdrp%ys

          wx=max(0.0_RP,min(1.0_RP,(rrat_p-x1)/adrpdrp%dx))
          wy=max(0.0_RP,min(1.0_RP,(log10(NreL_p)-y1)/adrpdrp%dy))

          E_c(i,j)=min(1.0_RP,max(ec_min, &
               (1.0_RP-wx)*(1.0_RP-wy)*drpdrp(i1,j1)+ &
               (1.0_RP-wx)*wy*drpdrp(i1,j1+1)+ &
               wx*(1.0_RP-wy)*drpdrp(i1+1,j1)+ &
               wx*wy*drpdrp(i1+1,j1+1)))

          if(E_c(i,j)>15.0_RP) then
            E_c(i,j)=15.0_PS
          endif
```

The other 6 tables (all in the token 2-1 ice-drop branch, `:16108–16244`) use exactly the same 4-line index+weight+gather pattern, differing only in axis variables:

- `hexdrp` (habit 1,5,6): `x=log10(Nre_rp)`, `y=log10(Nre_sp)` (`:16126–16142`)
- `bbcdrp` (habit 2,4): same axes (`:16145–16161`)
- `coldrp` (default habit): same axes (`:16164–16180`)
- `gp1drp` (den≤0.2): `x=rrat`, `y=log10(NreL_p)` where `rrat=g_2%len*0.5/max(semi_a,semi_c)` (`:16185–16204`)
- `gp4drp` (den≤0.6): `x=rrat`, `y=log10(NreL_p)` (`:16207–16223`)
- `gp8drp` (else): `x=rrat`, `y=log10(NreL_p)` (`:16226–16242`)

Canonical clamped-index formula (M2 target):
`j1 = max(1, min(nc-1, int((x - xs)/dx)+1))`, `i1 = max(1, min(nr-1, int((y - ys)/dy)+1))`, weights `wx = clamp01((x-x1)/dx)`, `wy = clamp01((y-y1)/dy)`, and `E = clamp(ec_min,1, (1-wx)(1-wy)T(i1,j1) + (1-wx)wy T(i1,j1+1) + wx(1-wy)T(i1+1,j1) + wx wy T(i1+1,j1+1))`.

---

## 4. Scatter helpers `add_simple_vec` / `add_samebin_vec`

### 4.1 `add_simple_vec` (`:15452–15577`) — bin-boundary search + scatter

Header + the target-bin search (`:15452–15499`):

```fortran
  subroutine add_simple_vec(nbin,L,icond,g,Np,Mp,ratio_Mp,new_N,new_M,new_Q,Qp)
    type (Group), intent(in)   :: g
    integer,intent(in) :: nbin,L
    integer,dimension(mxnbin+1,*),intent(in) :: icond
    real(PS),dimension(mxnbin+1,*),intent(inout) :: Np,Mp
    real(PS), dimension(mxnbin+1,L,mxnmasscomp)  :: ratio_Mp
    real(PS), dimension(mxnbin+1,L,mxnnonmc+2)   :: Qp
    real(8), dimension(mxnbin,*)                 :: new_N
    real(8), dimension(mxnbin,L,1+mxnmasscomp)   :: new_M
    real(8), dimension(mxnbin,L,mxnnonmc)        :: new_Q
    integer,dimension(nbin,L) :: jbin
    integer :: i,n,jj,j

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
```

The scatter with number-conservation clamp on `Np` (`:15505–15528`):

```fortran
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
    do jj=1,g%N_masscom
       do n = 1, L
       do i = 1, nbin
          if(icond(i,n)==1) then
             j=jbin(i,n)
             new_M(j,n,1+jj)=new_M(j,n,1+jj)+Mp(i,n)*ratio_Mp(i,n,jj)
          endif
       end do
       end do
    end do
```

The non-mass (`new_Q`) scatter follows for token 1/3 (`:15529–15544`) and token 2 (`:15545–15576`), the latter treating `iacr/iccr/idcr/iag/icg` axis-length components as `Np*Qp**3`.

### 4.2 `add_samebin_vec` (`:15579–15657`) — accumulate into same bin

```fortran
  subroutine add_samebin_vec(new_N,new_M,new_Q,nbin,L,icond,token,&
                    n_masscom,n_nonmass,Np,Mp,ratio_Mp,Qp)
    integer,intent(in) :: nbin,L,token,n_masscom,n_nonmass
    integer,dimension(mxnbin+1,*),intent(in) :: icond
    real(8),dimension(mxnbin+1,*),intent(in) :: Np,Mp
    real(PS), dimension(mxnbin+1,L,mxnmasscomp) :: ratio_Mp
    real(PS), dimension(mxnbin+1,L,mxnnonmc+2)  :: Qp
    real(8), dimension(mxnbin,*)                :: new_N
    real(8), dimension(mxnbin,L,1+mxnmasscomp)  :: new_M
    real(8), dimension(mxnbin,L,mxnnonmc)       :: new_Q
    integer :: i,n,jj

    do n = 1, L
    do i = 1, nbin
      if(icond(i,n)==1) then
        new_M(i,n,1)=new_M(i,n,1)+Mp(i,n)
        new_N(i,n)=new_N(i,n)+Np(i,n)
      endif
    enddo
    enddo
    do jj=1,N_masscom
       do n = 1, L
       do i = 1, nbin
          if(icond(i,n)==1) then
             new_M(i,n,1+jj)=new_M(i,n,1+jj)+Mp(i,n)*ratio_Mp(i,n,jj)
          endif
       end do
       end do
    end do
    if(token==2) then
      do jj=1,N_nonmass
         do n = 1, L
         do i = 1, nbin
            if(icond(i,n)==1) then
               if(jj.eq.ivcs.or.jj.eq.inex) then
                  new_Q(i,n,jj)=new_Q(i,n,jj)+Qp(i,n,jj)*Np(i,n)
               elseif(jj.eq.iacr.or.jj.eq.iccr.or.jj.eq.idcr.or.jj.eq.iag &
                    .or.jj.eq.icg) then
                  new_Q(i,n,jj)=new_Q(i,n,jj)+Np(i,n)*Qp(i,n,jj)**3
               else
                  LOG_ERROR("add_samebin_vec",*) "add_samebin: This parameter is not defined",jj
                  call PRC_abort
               endif
            end if
         end do
         end do
      end do
    else
      do jj=1,N_nonmass
         do n = 1, L
         do i = 1, nbin
            if(icond(i,n)==1) then
               new_Q(i,n,jj)=new_Q(i,n,jj)+Qp(i,n,jj)*Np(i,n)
            endif
         end do
         end do
      end do
    end if
  end subroutine add_samebin_vec
```

---

## 5. `cal_ratio_mass_col_vec` (`:17586–18036`) — mass-component ratios

Header + declarations (`:17586–17630`):

```fortran
  subroutine cal_ratio_mass_col_vec( g_1, g_2, ag, level, icond3, &
      n_all,Mc, dN_ice,ratio_M_2, dM_col, Mp, ratio_Mp)
    type (Group), intent(inout)    :: g_1, g_2
    type (AirGroup), intent(in) :: ag
    integer, intent(in)         :: level
    integer,intent(in),dimension(mxnbin+1,*) :: icond3
    real(PS), dimension(mxnbin+1,mxnmasscomp,*)            :: Mc
    integer,dimension(*),intent(in)         :: n_all
    real(PS), dimension(mxnbin+1,g_1%L,mxnmasscomp),intent(inout)   :: ratio_Mp
    real(PS), dimension(mxnbin,1+mxnmasscomp,*)            :: ratio_M_2
    real(PS), dimension(mxnbin+1,*)   ::  dN_ice
    real(PS), dimension(mxnbin+1,mxnbin,*)      :: dM_col
    real(PS),dimension(mxnbin+1,*)                    :: Mp
    real(PS),dimension(mxnbin+1,LMAX)       :: dM_1,dM_2
    integer                     :: i,j,n,n_all_max
    n_all_max=maxval(n_all(1:g_1%L))
```

**Warm rain-rain branch** (`g_1%token==1 .and. g_2%token==1`, `:17633–17687`). For `level>=4` it distributes collected total/soluble aerosol mass from `ratio_M_2` weighted by `dM_col`:

```fortran
    if( g_1%token == 1 .and. g_2%token == 1 ) then
      if(level>=4) then
        do n = 1, g_2%L
        do j = 1, n_all_max
           dM_1(j,n) = 0.0_PS   ! total aerosol
           dM_2(j,n) = 0.0_PS   ! soluble aerosol
        end do
        end do
        do i = 1, g_2%N_BIN
           do n = 1, g_2%L
           do j = 1, n_all_max
            if(icond3(j,n)==1.and.j<=n_all(n)) then
              if( dM_col(j,i,n) > 0.0_PS ) then
                dM_1(j,n)=dM_1(j,n)+ratio_M_2(i,rmat_m,n)*dM_col(j,i,n)
                dM_2(j,n)=dM_2(j,n)+ratio_M_2(i,rmas_m,n)*dM_col(j,i,n)
              end if
            endif
          enddo
          enddo
        end do
        do n = 1, g_2%L
        do j = 1, n_all_max
          if(icond3(j,n)==1.and.j<=n_all(n)) then
            ratio_Mp(j,n,rmat_m)=(Mc(j,rmat_m,n)+dM_1(j,n))/Mp(j,n)
            ratio_Mp(j,n,rmas_m)=(Mc(j,rmas_m,n)+dM_2(j,n))/Mp(j,n)
            ratio_Mp(j,n,rmai_m)=max(ratio_Mp(j,n,rmat_m)-ratio_Mp(j,n,rmas_m),0.0_PS)
          end if
        enddo
        enddo
      endif
```

(The `token 2-2` aggregation branch is `:17689–17874`; the `token 2-1/2-11` riming branch is `:17875–18025`; neither is on the warm path. For warm rain with `level<4` `ratio_Mp` stays zero — mass is pure water.)

---

## 6. `cal_needgive` (`class_Group.F90:9625–9817`) — inter-bin mass borrowing

Full verbatim. Group 1 is the higher group, group 2 the lower; negative-mass bins are repaired first from positive bins in the same group, then by borrowing from group 2, then from vapor.

```fortran
  subroutine cal_needgive(g_1, g_2, ngrid, m_t1, m_t2, need_m, give_m)
    ! group 1 is the higher group, and group 2 is lower one.
    ! (lower) - cloud drop - rain - solid_hydro (higher)
    type (group), intent(inout)         :: g_1, g_2
    integer, intent(in)       :: ngrid
    real(ps), intent(inout)      :: m_t1, m_t2
    real(ps)      :: nm_t2
    real(ps), intent(inout)   :: need_m, give_m
    real(ps), parameter           :: min_mt = 1.0e-20
    real(ps), parameter       :: ini_mcloud = (4.0_ps*pi/3.0_ps)*1.0e-9
    real(ps)   :: mod_dum, total_pos, total_neg
    real(ps),dimension(g_1%n_bin) :: old_mass,modrg
    integer    :: i,j,k
    integer,dimension(g_1%n_bin) :: iswitch

    ! calculate total positive and negative mass
    total_pos = 0.0_ps
    total_neg = 0.0_ps
    do i = 1, g_1%n_bin
      if( g_1%ms(i,ngrid)%mass(1) < 0.0_ps) then
        total_neg = total_neg + (-g_1%ms(i,ngrid)%mass(1))
      else
        total_pos = total_pos + g_1%ms(i,ngrid)%mass(1)
      end if
    end do

    if( total_neg == 0.0_ps .and. total_pos == 0.0_ps ) then
      g_1%mark_cm(ngrid) = 3
    else if( total_neg == 0.0_ps .and. total_pos > 0.0_ps ) then
      modrg=1.0
      call check_con(g_1,ngrid,modrg)
    else if( total_neg > 0.0_ps .and. total_pos > 0.0_ps ) then

      if( m_t1 >= min_mt ) then
        ! --- group can be fixed within the group ---
        iswitch=0
        do i=1,g_1%n_bin
          old_mass(i)=g_1%ms(i,ngrid)%mass(1)
        enddo
        do i=1,g_1%n_bin
          if( g_1%ms(i,ngrid)%mass(1) < 0.0 ) then
            g_1%ms(i,ngrid)%mark = 1
            mod_dum = (total_pos-(-g_1%ms(i,ngrid)%mass(1)))/total_pos
            total_pos = total_pos-(-g_1%ms(i,ngrid)%mass(1))
            do j = 1, g_1%n_bin
              if( i == j ) then
                iswitch(i)=1
                g_1%ms(j,ngrid)%mass(1)=0.0_PS
              else
                if( g_1%ms(j,ngrid)%mass(1) > 0.0_ps ) then
                  iswitch(i)=2
                  g_1%ms(j,ngrid)%mass(1)=g_1%ms(j,ngrid)%mass(1)*mod_dum
                  g_1%ms(j,ngrid)%mark = 1
                end if
              end if
            end do
          end if
        end do

        do i=1,g_1%n_bin
          modrg(i)=1.0_PS
          if(old_mass(i)>1.0e-30) then
            modrg(i)=g_1%ms(i,ngrid)%mass(1)/old_mass(i)
          endif
        enddo

        do i = 1, g_1%N_BIN
          if(iswitch(i)>0) then
             do k = 1, g_1%N_masscom
                g_1%ms(i,ngrid)%mass(1+k)=g_1%ms(i,ngrid)%mass(1+k)&
                     *modrg(i)
             enddo
          endif
        enddo

        ! recompute totals + concentration
        total_pos = 0.0_ps
        total_neg = 0.0_ps
        do i = 1, g_1%n_bin
          if( g_1%ms(i,ngrid)%mass(1) < 0.0_ps) then
            total_neg = total_neg + (-g_1%ms(i,ngrid)%mass(1))
          else
            total_pos = total_pos + g_1%ms(i,ngrid)%mass(1)
          end if
        end do
        call check_con(g_1,ngrid,modrg)

      else if( m_t1 < 0.0_ps ) then
        ! --- need to borrow mass from other groups ---
        iswitch=0
        do i = 1, g_1%n_bin
          modrg(i)=1.0
          if( g_1%ms(i,ngrid)%mass(1) < 0.0_ps) then
            iswitch(i)=1
            g_1%ms(i,ngrid)%mass(1) = 0.0_ps
            g_1%ms(i,ngrid)%con = 0.0_ps
            g_1%ms(i,ngrid)%mark = 2
            modrg(i)=0.0
          end if
        end do
        do i = 1, g_1%N_BIN
           if(iswitch(i)>0) then
              do k = 1, g_1%N_masscom
                 g_1%ms(i,ngrid)%mass(1+k)=0.0_PS
              enddo
          endif
        enddo
        call check_con(g_1,ngrid,modrg)

        if( m_t2 > min_mt ) then
          ! borrow from the second group
          if( total_neg <= m_t2 ) then
            nm_t2 = m_t2 - total_neg
          else if( total_neg > m_t2 ) then
            need_m = need_m + (total_neg - (m_t2-min_mt) )
            nm_t2 = min_mt
          end if
          mod_dum = nm_t2/m_t2
          m_t2 = nm_t2
          modrg=mod_dum
          do i = 1, g_2%N_BIN
          do k = 1, 1+g_2%N_masscom
            g_2%ms(i,ngrid)%mass(k)=g_2%ms(i,ngrid)%mass(k)*mod_dum
          end do
          end do
          call check_con(g_2,ngrid,modrg)
        else
          ! get mass from vapor
          need_m = need_m + total_neg
        end if
      else
        do i=1,g_1%n_bin
          g_1%ms(i,ngrid)%con=0.0_ps
          modrg(i)=0.0_PS
        end do
        do i = 1, g_1%N_BIN
        do k = 1, (1+g_1%N_masscom)
          g_1%ms(i,ngrid)%mass(k)=0.0_PS
        enddo
        enddo
        call check_con(g_1,ngrid,modrg)
        g_1%mark_cm(ngrid) = 3
        give_m = give_m + m_t1
      end if
    end if
  end subroutine cal_needgive
```

The `mark_cm(ngrid)=3` sentinel is exactly the flag that `coalescence` checks at `:1570–1573` to skip empty groups (`icycle_n(n)=1`).

---

## M2 porting notes (kernel structure recap)

- **Outer collection kernel** (`:1634`): triple loop `i(desc) × n × j`, guarded by `icycle_ijn==0`. Emit as SPLIT over `(i,j)` pairs per column `n`; the scatter target is a *scalar accumulation* into `used_N_2(j,n)`/`used_M_2(j,n)` (reduction over `i`) plus `N_col(i,j,n)` write. `N_col` clamped to `g_2%MS(j,n)%con`.
- **Collector scatter** (`collector_loop1`, `:1971`): for each collector bin `i`, the sub-bin set `n_all(n)` is variable (fortunate `n_f` + unfortunate 1 + high `n_h`), so the scatter is a gather over ragged sub-bins. Target original bin found by the linear-remap in `cal_transbin_vec` after `binb3d` boundaries `= g_1%binb(i)+dmass_min` / `g_1%binb(i+1)+dmass_max` and the `cal_lincubprms_vec`/`cal_linprms_vec_s` fit with `bbmf=0.2` widening on error. The simpler mean-mass boundary search (`add_simple_vec`, `:15475–15499`) is the fallback pattern: `jbin = the j with binb(j) < Mp/Np <= binb(j+1)`, with number clamped to `[Mp/(1.01·binb(j)), Mp/(0.99·binb(j+1))]`.
- **Breakup** is entirely additive at the end: `N_bk = N_col·max(0,1-E_coal)`, fragments read from the pre-tabulated `bu_fd(1:2, (i1d_pair-1)·N_BIN + k)` with `i1d_pair = j-jmin_bk+1 + (i-imin_bk)(1+i-imin_bk)/2` (triangular, `j<i` only), scaled by `mod_rat = (mean_mass_i + mean_mass_j)/bu_tmass(i1d_pair)`.