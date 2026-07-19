The extraction is complete. Below is the supplement doc filling G4's elisions and answering item 5.

---

# AMPS Collision–Coalescence — RAIN-RAIN — Gap-Fill Supplement to G4

Source: `/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_core.F90` (28357 lines); density round-trip in `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Group.F90`; caller in `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Cloud_Micro.F90`.

G4 (`docs/superpowers/facts/m2/coalescence.md`, 1347 lines) is correct and complete for §1.1–1.5, §2.1–2.4, §4, §5, §6. This supplement replaces every `...`/"map offset NNN"/"region" placeholder in G4 §1.6 and §3 with verbatim source, and adds the **missing item 5** (per-mass vs per-volume density convention), plus the exact rain-rain caller.

---

## 0. Rain-rain caller (`class_Cloud_Micro.F90:1015–1028`) — verbatim

```fortran
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
```

So for the warm case both group actuals are `CM%rain` (`g_1 === g_2`, `token==1`), `ibreak = CM%micexfg(18)` (breakup flag), and the returned diagnostics are `dM_auto_liq`/`dM_accr_liq`. Gated by `micexfg(2)==1` (collision-coalescence on) and `flagp_r>0` (rain present).

---

## 5. Per-mass vs per-volume — DENSITY CONVENTION (item 5, absent from G4)

**Answer: coalescence operates entirely on per-VOLUME quantities.** `gr%MS(i,n)%con` and `gr%MS(i,n)%mass(k)` are number- and mass-densities already multiplied by air density `ag%TV(n)%den` on entry to the microphysics step, and divided back out on exit. This is the OPPOSITE of the per-mass activation / vapor-deposition ports (which carry mixing ratios with no `/den`). M2b must therefore multiply the mixing-ratio state by `den` before the kernel and divide after, OR carry `den` through the kernel — coalescence itself contains no `den` factor.

**Entry conversion (state → MS, ×den), `class_Group.F90:756–758`:**

```fortran
          gr%MS(i,j)%mass(rmt)=XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den
          if(XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den>0.0) then
             gr%MS(i,j)%con=XR(rcon_q,i,NRCAT,j)*ag%TV(j)%den
```

Aerosol mass components (level≥4) likewise: `mass(rmat)=…XR(rmat_q,…)*ag%TV(j)%den…` (`:779–781`), `mass(rmas)=…*ag%TV(j)%den` (`:811`).

**Exit conversion (MS → state, ÷den), `class_Group.F90:3910–3921`:**

```fortran
          XR(rmt_q,i,nrcat,n)=gr%MS(i,n)%mass(rmt)/ag%tv(n)%den
          XR(rcon_q,i,nrcat,n)=gr%MS(i,n)%con/ag%tv(n)%den
             XR(rmat_q,i,nrcat,n)=gr%MS(i,n)%mass(rmat)/ag%tv(n)%den
             XR(rmas_q,i,nrcat,n)=gr%MS(i,n)%mass(rmas)/ag%tv(n)%den
```

(Ice group is identical: `class_Group.F90:888–890` in, `:3941–3960` out.) `XR(rcon_q,…)`/`XR(rmt_q,…)` are the mixing-ratio state (= `qrpv`-type). `%con = q_state * den`.

**The kernel confirms per-volume internally.** `KC(i,j)=E_c*(vtm_i−vtm_j)*A_c*con_j*dt` (`:16297–16298`) and `N_col(i,j,n)=con_i*KC(i,j,n)` (G4 §1.2, `:1638`). Expanding: `N_col ∝ con_i · con_j · kernel · dt` — a product of two per-volume number densities, i.e. the standard `n_i·n_j·K·dt` stochastic-collection rate. If M2b feeds mixing ratios (n/den) instead of number densities (n), each `con` is a factor `den` too small and `N_col` is `den²` too small; hence M2b must scale to per-volume for the kernel. No stray `/den` or `*den` appears anywhere inside `coalescence`, `cal_collision_kernel_func`, `cal_Coalescence_Efficiency`, or the scatter helpers — the density lives only in the `class_Group` wrapper.

---

## 3 (fill). `drpdrp` bilinear gather — VERBATIM including the `rrat<0` guard (`:16066–16107`)

G4 elided the `rrat<0` branch with `...`. Full loop:

```fortran
       do j = 1, g_2%N_BIN
       do i = 1, g_1%N_BIN
        if(icond1(i,j)==0) then
          NreL_p=max(g_1%MS(i,n)%Nre,1.0e-10_RP)
          rrat=g_2%MS(j,n)%len/g_1%MS(i,n)%len
          if(rrat<0.0_RP) then
            em(i,j)=1
            rrat_p=0.0_RP
            NreL_p=max(g_1%MS(i,n)%Nre,1.0e-10_RP)
            rrat=g_2%MS(j,n)%len/g_1%MS(i,n)%len
            if(debug) write(*,*) "something is not right getdrpdrp",i,j,rrat,NreL_p,&
                 g_2%MS(j,n)%len,g_1%MS(i,n)%len
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

!          call getdrpdrp(E_c(i,j),em(ij),real(g_1%MS(i,n)%Nre), real(g_2%MS(j,n)%len/g_1%MS(i,n)%len))
          if(E_c(i,j)>15.0_RP) then
            E_c(i,j)=15.0_PS
          endif
        endif
      enddo
      enddo
```

Note the branch opener at `:16061–16062` covers **both** `token 1-1` and `token 1-11`:
```fortran
    if( ( g_1%token == 1 .and. g_2%token == 1 ) .or. &
        ( g_1%token == 1 .and. g_2%token == 11 ) ) then
```
The `rrat<0` block is a degenerate guard (negative length ratio, "should not happen"): it sets `em=1`, forces `rrat_p=0`, then **recomputes `rrat` unchanged** and falls through to the same index math (so `x=rrat_p=0` → `j1=1`). The as_offset/tiled-table DSL concern: axis origins `%xs,%ys`, steps `%dx,%dy`, extents `%nc,%nr` come from the `col_lut_aux` descriptor at runtime; `i1,j1` are clamped to `[1, nr-1]×[1, nc-1]` and the gather touches the 2×2 stencil `drpdrp(i1:i1+1, j1:j1+1)`. `drpdrp` is `dimension(adrpdrp%nr,*)` (row-major first index = Reynolds axis `i1`, second = radius-ratio axis `j1`). The 6 ice tables (`:16108–16244`) reuse this identical 4-line index+weight+2×2-gather with different axis variables (G4 §3 bullet list is accurate).

Kernel assembly tail (`:16294–16310`) — verbatim, with commented Golovin/constant alternatives retained:

```fortran
    do j = 1, g_2%N_BIN
    do i = 1, g_1%N_BIN
      if(icond1(i,j)==0) then
        KC(i,j)=E_c(i,j)*( g_1%MS(i,n)%vtm - g_2%MS(j,n)%vtm )*A_c(i,j)* &
                g_2%MS(j,n)%con*g_1%dt
!!c        KC(i,j) = 1500.0_PS * (g_1%MS(i,n)%mean_mass + g_2%MS(j,n)%mean_mass) * g_2%MS(j,n)%con*g_1%dt ! Golovin kernel
!!c        E_coal(i,j) = 1.0D0 ! Golovin kernel
!!c        E_c(i,j) = 1.0_PS
!!c      KC(i,j) = 0.5_PS * g_2%MS(j,n)%con*g_1%dt ! constant kernel
      endif
    enddo
    enddo
```

---

## 1.6 (fill). `collector_loop1` scatter-to-target — VERBATIM (the hardest M2b kernel)

G4 §1.6 gave the categorization, the fortunate/unfortunate split, and named the "H"/"F"/"LM"/checker passes by "map offset". Here they are verbatim from source. Loop opens `:1971`:

```fortran
       collector_loop1: do i=g_1%N_bin,icolbin_min,-1
```

### 1.6a Mass-component seed `Mc` (`:2126–2154`, warm token-1 branch)

```fortran
         if(level>=4.and.g_1%token==1) then
            do n = 1, g_1%L
            do j = 1, n_all_max
             if(Mp(j,n)<1.0e-28_PS.and.j<=n_all(n)) then
               Mp(j,n)=0.0_PS
               Np(j,n)=0.0_PS
             endif
           enddo
           enddo
           do k = 1, g_1%N_masscom
              do n = 1, g_1%L
              do j = 1, n_all_max
                 if(j<=n_all(n)) then
                    if(Mp(j,n)<1.0e-28_PS) then
                       Mc(j,k,n)=0.0_PS
                    else
                       Mc(j,k,n) = Mp(j,n)*(g_1%MS(i,n)%mass(1+k)/g_1%MS(i,n)%mass(1))
                    endif
                 endif
              enddo
              enddo
           end do
```

Then `Mp_ini(j,n)=Mp(j,n)` snapshot (`:2200–2206`) for the autoconversion diagnostic.

### 1.6b High-rate collection pass "H" (`:2213–2260`) — VERBATIM

```fortran
         do jj = 1, n_h_max
            do n = 1, g_1%L
            do k = 1, n_all_max

             if(icond1(i,n)==1.and.jj<=n_h(n).and.k<=n_all(n)) then
               j=jbin_h(jj,n)
               if( col_ratio(j,n) >= 1.0_PS .and. col_ratio(j,n) <= 10.0_PS ) then
                 ndrop_B = aint(col_ratio(j,n))
               else
                 ndrop_B = col_ratio(j,n)
               end if
               dM_dum = ndrop_B * g_1%MS(i,n)%con * g_2%MS(j,n)%mean_mass * &
                        (Np(k,n)/max(left_N(i,n),1e-30_PS)) ! fraction of the subbin, total summation is one
               Mp(k,n) = Mp(k,n) + dM_dum
               dmass_max(k,n) = dmass_max(k,n) + ndrop_B * g_2%binb(j+1)
               dmass_min(k,n) = dmass_min(k,n) + ndrop_B * g_2%binb(j)

               dM_col(k,j,n) = dM_col(k,j,n) + dM_dum
               ! +++ add up what was used for coalescence process +++
               used_M_2(j,n) = used_M_2(j,n) + dM_dum
               if( used_M_2(j,n) >= g_2%MS(j,n)%mass(1) ) then
                 used_M_2(j,n) = g_2%MS(j,n)%mass(1)
                 used_marker(j,n) = 1
               end if

               N_drops(k,j,n) = N_drops(k,j,n) + ndrop_B
```

followed by the KiD auto/accr sub-classification (`:2243–2256`) keyed on `isplit_bin_liq`:

```fortran
               if(g_1%token==1.and.g_2%token==1) then
                 if(max(g_1%binb(i+1),g_2%binb(j+1))<=g_1%binb(isplit_bin_liq+1)) then
                   dMp_auto(k,n)=dMp_auto(k,n)+dM_dum
                   dmass_max_auto(k,n) = dmass_max_auto(k,n) + ndrop_B * g_2%binb(j+1)
                   dmass_min_auto(k,n) = dmass_min_auto(k,n) + ndrop_B * g_2%binb(j)
                 elseif(g_1%binb(i)>=g_1%binb(isplit_bin_liq+1).and.&
                        g_2%binb(j+1)<=g_1%binb(isplit_bin_liq+1)) then
                   dMp_accr(k,n)=dMp_accr(k,n)+dM_dum
                 endif
               endif
```

### 1.6c Fortunate pass "F" (`:2275–2322`) — VERBATIM

```fortran
         do n = 1, g_1%L
         do jj = 1, n_f_max
           if(icond1(i,n)==1.and.jj<=n_f(n)) then
             j = jbin_lm(jj,n)

             ndrop_F = 1.0_PS

             if( col_ratio(j,n) > 0.0_PS .and. col_ratio(j,n) < 1.0_PS ) then
               p_F = 0.0_PS
               dM_dum = N_col(i,j,n)*E_coal(i,j,n)*g_2%MS(j,n)%mean_mass
             else if( col_ratio(j,n) >= 1.0_PS .and. col_ratio(j,n) <= 10.0_PS ) then
               p_F = col_ratio(j,n)-aint(col_ratio(j,n))
               dM_dum = p_F * g_1%MS(i,n)%con * g_2%MS(j,n)%mean_mass
             else
               p_F = 0.0_PS
               dM_dum = 0.0_PS
             end if
             Mp(jj,n) = Mp(jj,n) + dM_dum
             dM_col(jj,j,n) = dM_col(jj,j,n) + dM_dum

             ! for fortunate growth group
             dmass_max(jj,n) = dmass_max(jj,n) + ndrop_F * g_2%binb(j+1)
             dmass_min(jj,n) = dmass_min(jj,n) + ndrop_F * g_2%binb(j)

             ! +++ add up what was used for coalescence process +++
             used_M_2(j,n) = used_M_2(j,n) + dM_dum
             if( used_M_2(j,n) >= g_2%MS(j,n)%mass(1) ) then
               used_M_2(j,n) = g_2%MS(j,n)%mass(1)
             end if

             N_drops(jj,j,n) = N_drops(jj,j,n) + ndrop_F
```
(same `dMp_auto/dMp_accr` block, `:2308–2317`, with `ndrop_F`.)

### 1.6d Left-over / low-medium pass "LM" (`:2335–2382`) — VERBATIM

```fortran
         do jj= (n_f_min+1), n_lm_max
            do n = 1, g_1%L
            do k = 1, n_all_max
             if(icond1(i,n)==1.and.k<=n_all(n).and.jj<=n_lm(n).and.jj>=n_f(n)+1) then
               j = jbin_lm(jj,n)
               if( col_ratio(j,n) > 0.0_PS .and. col_ratio(j,n) < 1.0_PS ) then
                 p_F = col_ratio(j,n)
               else if( col_ratio(j,n) >= 1.0_PS .and. col_ratio(j,n) <= 10.0_PS ) then
                 p_F = col_ratio(j,n)-aint(col_ratio(j,n))
               else
                 p_F = 0.0_PS
               end if
               ndrop_B = p_F

               dM_dum = p_F * g_1%MS(i,n)%con * g_2%MS(j,n)%mean_mass * &
                        (Np(k,n)/max(left_N(i,n),1e-30_PS))
               Mp(k,n) = Mp(k,n) + dM_dum
               dmass_max(k,n) = dmass_max(k,n) + ndrop_B * g_2%binb(j+1)
               dmass_min(k,n) = dmass_min(k,n) + ndrop_B * g_2%binb(j)

               dM_col(k,j,n) = dM_col(k,j,n) + dM_dum
               ! +++ add up what was used for coalescence process +++
               used_M_2(j,n) = used_M_2(j,n) + dM_dum
               if( used_M_2(j,n) >= g_2%MS(j,n)%mass(1) ) then
                 used_M_2(j,n) = g_2%MS(j,n)%mass(1)
               end if

               N_drops(k,j,n) = N_drops(k,j,n) + ndrop_B
```
(`dMp_auto/dMp_accr` block `:2367–2377`.)

### 1.6e Quality vars + no-collision rollback (`:2394–2457`)

Quality-var call (per-parent volume components) `:2394`:
```fortran
         call cal_Qp2_vec( g_1, g_2, icond1, i, n_all, &
                          N_drops, Qp, checker, dN_ice)
```
The `checker`-based rollback is **guarded `g_1%token==2 .and. g_2%token==1` (`:2428`) — NOT executed for warm rain** (G4 quoted this block but did not flag the guard). For rain-rain it is skipped entirely:
```fortran
         if( g_1%token == 2 .and. g_2%token == 1) then   ! 2014/10 tempei bug found
            do n = 1, g_2%L
            do k = 1, g_2%N_BIN
            do j = 1, n_all_max
             if(icond1(i,n)==1.and.j<=n_all(n)) then
               if( N_drops(j,k,n) > 0.0_PS .and. checker(j,k,n) == 0) then
                 dM_dum = N_drops(j,k,n)*&
                       g_1%MS(i,n)%con*g_2%MS(k,n)%mean_mass*(Np(j,n)/max(left_N(i,n),1e-30_PS))
                 Mp(j,n)=Mp(j,n)-dM_dum
                 dmass_max(j,n)=dmass_max(j,n)-N_drops(j,k,n)*g_2%binb(k+1)
                 dmass_min(j,n)=dmass_min(j,n)-N_drops(j,k,n)*g_2%binb(k)
                 used_N_2(k,n)=used_N_2(k,n)-N_drops(j,k,n)*g_1%MS(i,n)%con
                 used_M_2(k,n)=used_M_2(k,n)-dM_dum
               end if
             end if
           end do
           end do
           end do
         end if
```

### 1.6f `icond3` non-empty mask + `cal_ratio_mass_col_vec` (`:2468–2484`)

```fortran
         do n = 1, g_1%L
         do j = 1, n_all_max
           if(icond1(i,n)==1.and.j<=n_all(n)) then
             if( Np(j,n) > 1.0e-30_PS .and. Mp(j,n) > 1.0e-30_PS ) then
               icond3(j,n)=1
             endif
           endif
         enddo
         enddo
         ...
         ratio_Mp(:,:,:) = 0.0_PS
         call cal_ratio_mass_col_vec( g_1, g_2, ag, level, icond3, &
                 n_all, Mc, dN_ice, ratio_M_2, dM_col, Mp, ratio_Mp)
```

### 1.6g Shifted-boundary construction + linear-remap (`:2537–2608`) — VERBATIM

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

             a(1) = a2d(j,n,1)
             a(2) = a2d(j,n,2)
             a(3) = a2d(j,n,3)
             a(4) = a2d(j,n,4)
             call cal_linprms_vec_s(Npd(j,n),Mpd(j,n), &
                      binb3d(j,n,1),binb3d(j,n,2),a(:), &
                      error_number(j,n))
           endif
         enddo
         enddo
```

Second-chance remap when still erroring (`:2592–2608`) — resets `Npd = Mpd/(0.5·(binb(i+1)+binb(i)))` (mean-mass-preserving number) then retries `cal_linprms_vec_s`:

```fortran
         do n = 1, g_1%L
         do j = 1, n_all_max
           if(1<=error_number(j,n).and.error_number(j,n)<=4) then
             Npd(j,n) = Mpd(j,n)/(0.5_PS*(g_1%binb(i+1)+g_1%binb(i)))
             Np(j,n)=Npd(j,n)
             a(1) = a2d(j,n,1)
             a(2) = a2d(j,n,2)
             a(3) = a2d(j,n,3)
             a(4) = a2d(j,n,4)
             call cal_linprms_vec_s(Npd(j,n),Mpd(j,n), &
                      binb3d(j,n,1),binb3d(j,n,2),a(:), &
                      error_number(j,n))
           endif
         enddo
         enddo
```

If still failing, `icond3(j,n)=0` (`:2624`) → `error_number=10` sentinel (`:2730`) so `cal_transbin_vec` drops that sub-bin.

The warm branch has an **empty** `elseif( g_1%token == 1 .and. g_2%token == 1 )` at `:2717–2719` (the `token==2` apparent-density/`cal_xxx_p_v5_vec` block `:2649–2716` is ice-only, skipped for rain).

### 1.6h Scatter into original bins — main `cal_transbin_vec` (`:2746–2758`) — VERBATIM

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

### 1.6i Autoconversion-diagnostic second scatter (`:2770–2875`) — VERBATIM (warm-only)

Guarded `if(g_1%token==1.and.g_2%token==1)`. Rebuilds `binb3d` from the `_auto` shift accumulators and `Mp_ini+dMp_auto`, refits, and scatters to scratch `new_*_tmp`; then integrates `dM_auto`/`dM_accr`:

```fortran
            do n = 1, g_1%L
            do j = 1, n_all_max
             binb3d(j,n,1)=g_1%binb(i)
             binb3d(j,n,2)=g_1%binb(i+1)
             Npd(j,n)=0.0d+0
             Mpd(j,n)=0.0d+0
             if(icond3(j,n)==1.and.j<=n_all(n)) then
               if(Np(j,n)>1.0e-30_PS.and.dMp_auto(j,n)>1.0e-30_PS) then
                 binb3d(j,n,1)=g_1%binb(i)  +dmass_min_auto(j,n)
                 binb3d(j,n,2)=g_1%binb(i+1)+dmass_max_auto(j,n)
                 Npd(j,n)=Np(j,n)
                 Mpd(j,n)=Mp_ini(j,n)+dMp_auto(j,n)
               endif
             endif
           enddo
           enddo
           call cal_lincubprms_vec(mxnbin+1,n_all_max,g_1%L,Npd &
                               ,Mpd,binb3d &
                               ,a2d,error_number,"auto1")
           ...  ! same bbmf-widen retry, then:
           call cal_transbin_vec(g_1%token &
                           ,g_1%L,g_1%N_masscom &
                           ,g_1%N_bin,n_all_max &
                           ,g_1%binb &
                           ,error_number &
                           ,a2d,binb3d,mtend &
                           ,new_N_tmp,new_M_tmp,new_Q_tmp &
                           ,new_mtend_tmp &
                           ,ratio_Mp,den_ip_p,axr_p,spx_p &
                           ,habit_p,den_ic_p &
                           ,rag_p,rcg_p,n_exice_p &
                           ,actINF_p &
                           ,0)
           do j=isplit_bin_liq+1,g_1%n_bin
             do n=1,g_1%L
               dM_auto(n)=dM_auto(n)+new_M_tmp(j,n,1)
             enddo
           enddo
           do j=1,n_all_max
             do n=1,g_1%L
               if(icond3(j,n)==1.and.j<=n_all(n)) then
                 if(Np(j,n)>1.0e-30_PS.and.dMp_accr(j,n)>1.0e-30_PS) then
                   dM_accr(n)=dM_accr(n)+dMp_accr(j,n)
                 endif
               endif
             enddo
           enddo
         endif
```

### 1.6j Breakup fragments + loop close (`:2888–2895`) — VERBATIM

```fortran
         if(pro_type==2.and.ibreak==1) then
           call add_fragments_col_vec(new_N_1,new_M_1, &
                    i,g_1,g_2,icond1,used_marker,E_coal,N_col, &
                    imin_bk,imax_bk,jmin_bk,jmax_bk,bu_tmass,bu_fd)
         endif
       enddo collector_loop1
```

### 1.6k Post-loop left-over re-add + finalize (`:2899–2938`) — VERBATIM

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
       call check_csvolume_vec(g_1,new_N_1,new_M_1,new_Q_1)
       call assign_tendency_vec( g_1, pro_type,&
              new_N_1, new_M_1, new_Q_1, new_mtend, &
              g_2, used_N_2, used_M_2,ratio_M_2)
```

`assign_tendency_vec` (def `:20891`) converts the `new_*` per-volume accumulators back into per-bin tendencies on `g_1%MS`; no `den` scaling inside it — the `/den` happens only at `class_Group.F90:3910–3921` when MS is written back to the mixing-ratio state `XR`.

---

## M2b codegen notes (delta over G4's recap)

- **Bin-pair kernel** (G4 §1.2, `:1634`): SPLIT over `(i desc, j)` per column `n`; `N_col(i,j,n)=con_i·KC`, clamped to `con_j`; reduction accumulate into `used_N_2(j,n)`/`used_M_2(j,n)`. For warm rain only the `KC>0` branch fires (no `abs(KC)` riming).
- **Scatter (collector_loop1)** is 4 passes writing the SAME sub-bin accumulators `Mp(k,n)`, `dmass_min/max(k,n)`, `dM_col(k,j,n)`, `N_drops(k,j,n)`, `used_M_2(j,n)`: pass-H (`ndrop_B=aint` or full ratio), pass-F (`ndrop_F=1`, `dM_dum=N_col·E_coal·mean_mass`), pass-LM (fractional `p_F`). All three use `dM_dum ∝ con_i·mean_mass_j·(Np_k/left_N_i)`. Target original bin is found NOT by a direct boundary search here but by `cal_transbin_vec` after remapping the sub-bin PDF onto shifted boundaries `binb3d = [binb(i)+dmass_min, binb(i+1)+dmass_max]`. The `add_simple_vec` boundary-search (G4 §4.1) is the reference for the simpler mean-mass `jbin` search M2 may emit as a fallback.
- **Warm-path guards to hardcode**: skip the `checker` rollback (`token 2-1` only), skip `cal_xxx_p_v5_vec`/apparent-density (`token 2`), run the autoconv-diagnostic second `cal_transbin_vec` (`token 1-1`), and run `add_fragments_col_vec` only when `pro_type==2 .and. ibreak==1`.
- **Density**: multiply mixing ratios by `den` before the kernel, divide after (round-trip lives in `class_Group`, not in `coalescence`). Unlike the per-mass activation/vapor-dep ports, coalescence is per-volume end-to-end.