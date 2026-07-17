# SCALE↔AMPS SI/CGS Conversions, State Packing, and M1 State-Bundle Requirements

All line numbers are from the **post-instrumentation** working tree (branch `cloudlab`). Property-index values from `/Users/jcanton/projects/scale_amps/contrib/AMPS/par_amps.F90` (defaults; "redefined at the initialization based on the configuration"):

```fortran
! par_amps.F90:17-28 (ice qipv)
integer :: imt_q=1
integer :: imr_q=10,ima_q=15,imc_q=11,imw_q=12,imf_q=16
integer :: imat_q=13,imas_q=14
integer :: icon_q=2
integer :: ivcs_q=3,iacr_q=4,iccr_q=5,idcr_q=6
integer :: iag_q=7,icg_q=8,inex_q=9
! par_amps.F90:41-43 (liquid qrpv)      ! par_amps.F90:52-54 (aerosol qapv)
integer :: rmt_q=1                      integer :: amt_q=1
integer :: rmat_q=3,rmas_q=4            integer :: acon_q=2
integer :: rcon_q=2                     integer :: ams_q=3
! par_amps.F90:32-37 (Group%MS%mass and IS non-mass slots)
integer :: imt=1
integer :: imr=2,ima=3,imc=4,imw=8,imf=9, imat=5,imas=6,imai=7
integer :: ivcs=1,iacr=2,iccr=3,idcr=4,iag=5,icg=6,inex=7
integer :: rmt=1 ; rmat=2,rmas=3,rmai=4 ; amt=1 ; ams=2,ami=3
```

---

## 1. Pack block `Z_LOOP_01` — SCALE `QTRC` → AMPS `qrpv/qipv/qapv`

File: `/Users/jcanton/projects/scale_amps/scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90`, lines **1625–1828**.

### 1.1 factor_mxr1/2 and thermo state (1625–1672)

```fortran
       Z_LOOP_01: do k = KS, KE
          factor_mxr1 = (QDRY(k,i,j) + &
                         QTRC(k,i,j,I_QV)) ! moist air mixing ratio
          factor_mxr2 = (QDRY(k,i,j)*DENS(k,i,j) + &
                         QTRC(k,i,j,I_QV)*DENS(k,i,j))
          ! rho[kg m-3] = 0.001 rho[g cm-3]
          ...
          ! pressure
          ptotv(k) = PRES(k,i,j)
          ! temperature
          tv(k) = TEMP(k,i,j)
          ! potential temperature
          thv(k) = tv(k)*(PRE00/ptotv(k))**(Rdry/CPdry)
          ! exner function
          piv(k) = tv(k)/thv(k)*CPdry
          ! unknown
          pbv(k) = 0.0_RP
          ! density
          moist_denv(k) = DENS(k,i,j)*factor_mxr1
          ! vapor mixing ratio
          qvv(k) = QTRC(k,i,j,I_QV)/factor_mxr1
          qtotal(k) = qtotal(k) + qvv(k)*factor_mxr2
          Emoist(k,1) = - LHV0 * QTRC(k,i,j,I_QV) * DENS(k,i,j)
          ! virtual potential temperature
          thetav(k) = thv(k)*(1.0e+0+0.61*qvv(k))
          ! w at full grid
          wbv(k) = W(k,i,j)
          ! w at half grid
          momv(k) = MOMZ(k,i,j)
```

All AMPS mixing ratios are **per moist air** (SCALE tracers are per total mass incl. all hydrometeors), hence the `/factor_mxr1`. Density conversion kg m⁻³ ↔ g cm⁻³ (×0.001) is deferred: the interface passes SI `moist_denv` and only divides **non-mass (per-number/PPV-per-volume) properties** by `0.001_RP`; the g cm⁻³ conversion of density itself happens inside AMPS (`ini_AirGroup`, §5.2).

### 1.2 Liquid spectrum (1676–1702)

```fortran
          do icr = 1, ncr ! assume to be 1
             ipr_qpr = 0
             do ibr = 1, nbr
                qrpv(rmt_q,ibr,icr,k) = ( QTRC(k,i,j,I_QL+ibr-1) + QTRC(k,i,j,I_QPPVL+ipr_qpr) ) / factor_mxr1
                do ipr = rmat_q, rmas_q
                   qrpv(ipr,ibr,icr,k) = QTRC(k,i,j,I_QPPVL+ipr_qpr)/factor_mxr1
                   ipr_qpr = ipr_qpr + 1
                enddo
                ! number concentration, non-mass variable
                qrpv(rcon_q,ibr,icr,k) = QTRC(k,i,j,I_QPPVL+ipr_qpr)/factor_mxr1/0.001_RP
                ipr_qpr = ipr_qpr + 1
                do ipr = npr-1, npr
                   qrpv(ipr,ibr,icr,k) = 0.0
                enddo
                mmassrv(ibr,icr,k) = 0.0
                if(qrpv(rmt_q,ibr,icr,k)>1.0e-25) then
                   k1r = min(k1r,k-KS+1)
                   k2r = max(k2r,k-KS+1)
                endif
                qrov(ibr,icr,k) = qrpv(rmt_q,ibr,icr,k)
                ! total water content
                qtotal(k) = qtotal(k) + (qrpv(rmt_q,ibr,icr,k) - qrpv(rmat_q,ibr,icr,k))*factor_mxr2
             enddo
          enddo
```

**Aerosol-mass-into-drop-mass addition (liquid):** `rmt_q` (total drop mass) = `I_QL` (pure water) **+** `I_QPPVL+0` — and since `ipr_qpr=0` corresponds to the first PPVL tracer, which is `rmat_q` (total aerosol mass in drops), total drop mass = water + in-drop aerosol total mass. Conversely pure water = `rmt_q - rmat_q` (see 1699 and the unpack).
**PPVL tracer order per bin:** `[rmat_q, rmas_q]` mass PPVs, then number concentration (`rcon_q`).
**0.001 factor:** only `rcon_q` (number concentration) gets `/0.001_RP` (mass PPVs do not).

### 1.3 Ice spectrum (1704–1791)

```fortran
          do ici = 1, nci ! assumed to be 1
             ipi_qi = 0
             do ibi = 1, nbi
                ! total ice mass
                qipv(imt_q,ibi,ici,k) = ( QTRC(k,i,j,I_QI+ibi-1) &
                                        + QTRC(k,i,j,I_QW+ibi-1) &
                                        + QTRC(k,i,j,I_QPPVI+ipi_qi+4) ) / factor_mxr1
                ! rimed, aggregate, and crystal
                do ipi = imr_q, imc_q
                   qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1
                   ipi_qi = ipi_qi + 1
                enddo
                ! melt water
                qipv(imw_q,ibi,ici,k) = QTRC(k,i,j,I_QW+ibi-1)/factor_mxr1
                ! frozen water
                qipv(imf_q,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1
                ipi_qi = ipi_qi + 1

                ! aerosol mass
                do ipi = imat_q, imas_q
                   qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1
                   ipi_qi = ipi_qi + 1
                enddo

                ! number concentration and other non-mass variables
                if (l_gaxis_version == 1) then
                   ! -- ORIGINAL

                   do ipi = imas_q+1, npi-2
                      qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1/0.001_RP
                      ipi_qi = ipi_qi + 1
                   enddo

                else if (l_gaxis_version == 2) then
                   ! -- METHOD 1
                   do ipi = imas_q+1, npi-2
                      qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1/0.001_RP
                      ipi_qi = ipi_qi + 1
                   enddo

                   if (qipv(icon_q,ibi,ici,k) < 1.e-22_RP) then

                      qipv(iccr_q,ibi,ici,k) = 0.0_RP
                      qipv(idcr_q,ibi,ici,k) = 0.0_RP
                      qipv(icg_q,ibi,ici,k) = 0.0_RP

                   else
                      qipv(iccr_q,ibi,ici,k) = qipv(iccr_q,ibi,ici,k)*qipv(iacr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                      qipv(idcr_q,ibi,ici,k) = qipv(idcr_q,ibi,ici,k)*qipv(iacr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                      qipv(icg_q,ibi,ici,k) = qipv(icg_q,ibi,ici,k)*qipv(iag_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                   endif

                else if (l_gaxis_version == 3) then
                   ! -- METHOD 2
                   do ipi = imas_q+1, npi-2
                      qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1/0.001_RP
                      ipi_qi = ipi_qi + 1
                   enddo
                   if (qipv(icon_q,ibi,ici,k) < 1.e-22_RP) then
                      qipv(iag_q,ibi,ici,k) = 0.0_RP
                      qipv(icg_q,ibi,ici,k) = 0.0_RP
                   else
                      qipv(iag_q,ibi,ici,k) = qipv(iag_q,ibi,ici,k)*qipv(iacr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                      qipv(icg_q,ibi,ici,k) = qipv(icg_q,ibi,ici,k)*qipv(iccr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                   endif

                endif

             enddo

             do ibi = 1, nbi
                do ipi = npi-1, npi
                   qipv(ipi,ibi,ici,k) = 0.0_RP
                enddo
                mmassiv(ibi,ici,k)=0.0_RP
                if(qipv(imt_q,ibi,ici,k)>1.0e-25) then
                   k1i=min(k1i,k)
                   k2i=max(k2i,k)
                endif
                qiov(ibi,ici,k)=qipv(imt_q,ibi,ici,k)
                ! total water content
                qtotal(k) = qtotal(k) + (qipv(imt_q,ibi,ici,k) - qipv(imat_q,ibi,ici,k))*factor_mxr2
             enddo

          enddo
```

**Aerosol-mass-into-ice-mass addition:** `imt_q` = `I_QI` (pure ice) + `I_QW` (melt water tracer) + `I_QPPVI+ipi_qi+4` where at that point `ipi_qi=0`, so `+4` is the **5th PPVI tracer = imat_q** (total in-ice aerosol mass) — PPVI order per bin: `[rime, (agg,) crystal] (imr_q..imc_q)`, `frozen (imf_q)`, `[imat_q, imas_q]`, then non-mass `[icon_q, ivcs_q, iacr_q, iccr_q, idcr_q, iag_q, icg_q, inex_q]` (`imas_q+1 .. npi-2`).
**0.001 factor (ice):** applies to **all non-mass properties** `ipi = imas_q+1 .. npi-2`, i.e. number concentration `icon_q`, circumscribing volume `ivcs_q`, cubed axis lengths `iacr_q/iccr_q/idcr_q`, cubed CoG coords `iag_q/icg_q`, and `inex_q`.
Versions 2/3 additionally re-derive the pseudo-axes on pack: v2 stores in QTRC `c³/a³`, `d³/a³`, `cg³/ag³` ratios and multiplies back (`×iacr_q/icon_q`, `×iag_q/icon_q`); v3 stores `ag³/a³`, `cg³/c³` and multiplies back (`×iacr_q/icon_q`, `×iccr_q/icon_q`); both zero them when `icon_q < 1.e-22_RP`.

### 1.4 Aerosol spectrum (1793–1810)

```fortran
          !
          ! Initialize aerosol spectrum
          !
          ipa_qpa = 0
          do ica = 1, nca
             do iba = 1, nba
                !
                !do ipa=1,npa-2
                qapv(amt_q,iba,ica,k) = QTRC(k,i,j,I_QPPVA+ipa_qpa)/factor_mxr1
                qapv(acon_q,iba,ica,k) = QTRC(k,i,j,I_QPPVA+ipa_qpa+1)/factor_mxr1/0.001_RP
                qapv(ams_q,iba,ica,k) = QTRC(k,i,j,I_QPPVA+ipa_qpa+2)/factor_mxr1
                ipa_qpa = ipa_qpa + 3
                !enddo
                do ipa = npa-1, npa
                   qapv(ipa,iba,ica,k)=0.0_RP
                enddo
             enddo
          enddo
```

PPVA order per bin: `[amt (total mass), acon (number, ×1/0.001), ams (soluble mass)]`.

### 1.5 Emoist "before" state (1658, 1813–1826)

```fortran
          Emoist(k,1) = - LHV0 * QTRC(k,i,j,I_QV) * DENS(k,i,j)      ! line 1658
...
         if (.not. l_no_ice_heat) then
            do ibi = 1, nbi
               Emoist(k,1) = Emoist(k,1) + LHF0 * QTRC(k,i,j,I_QI+ibi-1) * DENS(k,i,j)
            end do
         else
            liquid_heat(k,1) = 0.0_RP
            do ibr = 1, nbr
               liquid_heat(k,1) = liquid_heat(k,1) + LHF0 * QTRC(k,i,j,I_QL+ibr-1) * DENS(k,i,j)
            end do
            ice_heat(k,1) = 0.0_RP
            do ibi = 1, nbi
               ice_heat(k,1) = ice_heat(k,1) + LHF0 * QTRC(k,i,j,I_QI+ibi-1) * DENS(k,i,j)
            end do
         endif
```

Underground extrapolation (1829–1834):

```fortran
       ! set underground, this is used for surface flux
       tv(KS-1) = tv(KS)
       pbv(KS-1) = 0.0_RP
       qvv(KS-1) = qvv(KS)
       ptotv(KS-1) = ptotv(KS)*exp( GRAV/Rdry/(tv(KS-1)*(1.0_RP + 0.61_RP*qvv(KS-1))) * ( CZ(KS,i,j) - CZ(KS-1,i,j) ) )
       moist_denv(KS-1) = ptotv(KS-1)/(Rdry*tv(KS-1)*(1.0_RP+0.61_RP*qvv(KS-1)))
```

---

## 2. Unpack / tendency block — AMPS state → `RHOQ_t`, `RHOE_t`, `CPtot_t/CVtot_t`

Same file, lines **2676–3010** (tendency), **2305–2352** (energy).

### 2.1 Conversion notebook (comment, 2676–2697)

```fortran
       ! tendency calculation notebook:
       ! X(n+1) = [ X(n) DEN(n) + RHOQ dt ] / DEN(n+1)
       ! RHOQ = [ X(n+1) DEN(n+1) - X(n) DEN(n) ] / dt
       ! X is SCALE variable, x is AMPS variable
       ! -----  mass  -----
       ! x is mixing ratio, its unit in AMPS is g cm-3/g cm-3
       ! X is mixing ratio, its unit in SCALE is kg m-3/kg m-3
       ! X is SCALE variable, x is AMPS variable
       ! no conversion needs to be done, as they are unitless, i.e. x = X
       ! thus tendency is RHOQ = [ x(n+1) DEN(n+1)  - X(n) DEN(n) ] / dt
       ! -----PPV mass-----
       ! x is mixing ratio, its unit in AMPS is g cm-3/g cm-3
       ! X is density, its unit in SCALE is kg m-3
       ! X is SCALE variable, x is AMPS variable
       ! conversion from x -> X must be performed by multiplying DEN (kg m-3)
       ! thus tendency is RHOQ = [ x(n+1) DEN(n+1) DEN(n+1) - X(n) DEN(n) ] / dt
       ! -----non-mass-----
       ! x is PPV per unit volume (cm-3) per density of moist air DEN (g cm-3)
       ! X is PPV per unit volume of cm-3, it is treated as a tracer
       ! conversion from x -> X must be performed by multiplying DEN (g cm-3)
       ! 1 kg m-3 = 0.001 g cm-3, thus X(n+1) = x(n+1)*DEN(n+1)*0.001
       ! thus tendency is RHOQ = [ x(n+1) DEN(n+1) 0.001 DEN(n+1)  - X(n) DEN(n) ] / dt
```

### 2.2 Density + vapor (2704–2713)

```fortran
       do k = KS, KE
          DENS_t(k,i,j) = DENS_t(k,i,j) + den_t(k)
          DENS_NEW(k) = DENS(k,i,j) + DENS_t(k,i,j)*dt
       enddo

       ! vapor tendency, mixing ratio -> density
       do k = KS, KE
          RHOQ_t(k,i,j,I_QV) = (qvv(k)*moist_denv(k) - &
               QTRC(k,i,j,I_QV)*DENS(k,i,j))/dt
       enddo
```

### 2.3 Liquid (2715–2737)

```fortran
       ! liquid tendency
       do k = KS, KE
          do icr = 1, ncr
             ipr_qpr = 0
             do ibr = 1, nbr
                RHOQ_t(k,i,j,I_QL+ibr-1) = ((qrpv(rmt_q,ibr,icr,k) - &
                     qrpv(rmat_q,ibr,icr,k))*moist_denv(k) -   &
                     QTRC(k,i,j,I_QL+ibr-1)*DENS(k,i,j))/dt

                do ipr = rmat_q, rmas_q
                   RHOQ_t(k,i,j,I_QPPVL+ipr_qpr) = &
                        (qrpv(ipr,ibr,icr,k)*moist_denv(k) - &
                        QTRC(k,i,j,I_QPPVL+ipr_qpr)*DENS(k,i,j))/dt
                   ipr_qpr = ipr_qpr + 1
                enddo
                ! number concentration, non-mass variable
                RHOQ_t(k,i,j,I_QPPVL+ipr_qpr) = &
                     (qrpv(rcon_q,ibr,icr,k)*moist_denv(k)*0.001_RP - &
                     QTRC(k,i,j,I_QPPVL+ipr_qpr)*DENS(k,i,j))/dt
                ipr_qpr = ipr_qpr + 1
             enddo
          enddo
       enddo
```

Note the exact inverse of the pack: `I_QL` gets back `rmt_q - rmat_q` (aerosol mass subtracted out of drop mass), and `rcon_q` re-multiplied by `0.001_RP`.

### 2.4 Ice — mass part + version-1 non-mass, verbatim (2739–2801)

```fortran
       ! ice tendency
       do k = KS, KE
          do ici = 1, nci
             ipi_qi = 0
             do ibi = 1, nbi
                ! total ice mass
                RHOQ_t(k,i,j,I_QI+ibi-1) = &
                     ((qipv(imt_q,ibi,ici,k) - qipv(imat_q,ibi,ici,k) - &
                     qipv(imw_q,ibi,ici,k))*moist_denv(k) - &
                     QTRC(k,i,j,I_QI+ibi-1)*DENS(k,i,j))/dt

                ! melt water mass
                RHOQ_t(k,i,j,I_QW+ibi-1) = &
                     (qipv(imw_q,ibi,ici,k)*moist_denv(k) - &
                     QTRC(k,i,j,I_QW+ibi-1)*DENS(k,i,j))/dt

                ! rimed ,aggregate, and crystal
                do ipi = imr_q, imc_q
                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(ipi,ibi,ici,k)*moist_denv(k) - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1
                enddo
                ! frozen water
                RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                     (qipv(imf_q,ibi,ici,k)*moist_denv(k) - &
                     QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                ipi_qi = ipi_qi + 1

                ! aerosol mass
                do ipi = imat_q, imas_q
                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(ipi,ibi,ici,k)*moist_denv(k) - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1
                enddo

                ! number concentration and other non-mass variables
                if (l_gaxis_version == 1) then
                   ! -- ORIGINAL
                   do ipi = imas_q+1, npi-2
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (qipv(ipi,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt

                      if ( ipi == iag_q .and. l_axis_limit .and. qipv(iacr_q,ibi,ici,k) > EPS ) then
                         if ( qipv(iag_q,ibi,ici,k) / qipv(iacr_q,ibi,ici,k) > 1.0_RP) then
                            RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                                 (qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP &
                                 - QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                         endif
                      endif

                      if ( ipi == icg_q .and. l_axis_limit .and. qipv(iccr_q,ibi,ici,k) > EPS ) then
                         if ( qipv(icg_q,ibi,ici,k) / qipv(iccr_q,ibi,ici,k) > 1.0_RP) then
                            RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                                 (qipv(iccr_q,ibi,ici,k)*moist_denv(k)*0.001_RP &
                                 - QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                         endif
                      endif

                      ipi_qi = ipi_qi + 1
                   enddo
```

`I_QI` gets back **pure ice** = `imt_q - imat_q - imw_q` (aerosol and melt water subtracted).

**Version 2 summary (2803–2879):** writes non-mass tracers individually in order `icon_q, ivcs_q, iacr_q`, then — if `iacr_q < 1.e-22` — zeros the next two tracers, else stores **ratios** `iccr_q*icon_q/iacr_q` and `idcr_q*icon_q/iacr_q` (i.e. QTRC holds c³/a³- and d³/a³-scaled quantities). Then `iag_q` (clipped to `iacr_q` when `l_axis_limit` and `(iag_q/iacr_q)^(1/3) > 1`; note: in the un-clipped else-branch `ipi_qi` increments but in the clipped branch it does not), then — if `iag_q < 1.e-22` — zero, else `iccr_q*icon_q/iag_q` (when cg-limit trips) or `icg_q*icon_q/iag_q`, finally `inex_q`. All with `*moist_denv(k)*0.001_RP`.
**Version 3 summary (2881–2949):** stores `icon_q, ivcs_q, iacr_q, iccr_q, idcr_q` directly, then ratio tracers `iag_q*icon_q/iacr_q` (or `icon_q` if the `l_axis_limit` `(iag_q/iacr_q)^(1/3)>1` cap trips; 0 if `iacr_q<1.e-22`) and `icg_q*icon_q/iccr_q` (analogous with `iccr_q`), then `inex_q`. All with `*moist_denv(k)*0.001_RP`.

### 2.5 Aerosol (2957–2974)

```fortran
       ! aerosol tendency
       do k = KS, KE
          ipa_qpa = 0
          do ica = 1, nca
             do iba = 1, nba
                RHOQ_t(k,i,j,I_QPPVA+ipa_qpa) = &
                      (qapv(amt_q,iba,ica,k)*moist_denv(k) - &
                      QTRC(k,i,j,I_QPPVA+ipa_qpa)*DENS(k,i,j))/dt
                RHOQ_t(k,i,j,I_QPPVA+ipa_qpa+1) = &
                      (qapv(acon_q,iba,ica,k)*moist_denv(k)*0.001_RP - &
                      QTRC(k,i,j,I_QPPVA+ipa_qpa+1)*DENS(k,i,j))/dt
                RHOQ_t(k,i,j,I_QPPVA+ipa_qpa+2) = &
                      (qapv(ams_q,iba,ica,k)*moist_denv(k) - &
                      QTRC(k,i,j,I_QPPVA+ipa_qpa+2)*DENS(k,i,j))/dt
                ipa_qpa = ipa_qpa + 3
             enddo
          enddo
       enddo
```

### 2.6 CPtot_t / CVtot_t (2976–3010)

```fortran
       ! specific heat tendency, loop over liquid and ice mass tendency, (mixing ratio)
       do k = KS, KE
          ! vapor difference
          dq = qvv(k) * moist_denv(k) / DENS_NEW(k) - QTRC(k,i,j,I_QV)
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CP_VAPOR * dq / dt
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CV_VAPOR * dq / dt
 
          ! liquid difference
          dq = 0.0_RP
          do ibr = 1, nbr
             ! liquid drop mass
             dq = dq + ( qrpv(rmt_q,ibr,1,k) - qrpv(rmat_q,ibr,1,k) ) * moist_denv(k) / DENS_NEW(k) &
                      - QTRC(k,i,j,I_QL+ibr-1)
          enddo
          do ibi = 1, nbi
             ! melt water mass
             dq = dq  + qipv(imw_q,ibi,1,k) * moist_denv(k) / DENS_NEW(k) - QTRC(k,i,j,I_QW+ibi-1)
          enddo
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CP_WATER * dq / dt
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CV_WATER * dq / dt
 
          ! ice difference
          dq = 0.0_RP
          do ibi = 1, nbi
             dq = dq + ( qipv(imt_q,ibi,1,k) - qipv(imw_q,ibi,1,k) - qipv(imat_q,ibi,1,k) ) * moist_denv(k) / DENS_NEW(k) &
                      - QTRC(k,i,j,I_QI+ibi-1)
          enddo
          if (.not. l_no_ice_heat) then
            CPtot_t(k,i,j) = CPtot_t(k,i,j) + CP_ICE * dq / dt
            CVtot_t(k,i,j) = CVtot_t(k,i,j) + CV_ICE * dq / dt
          else
            CPtot_t(k,i,j) = CPtot_t(k,i,j) + CP_WATER * dq / dt
            CVtot_t(k,i,j) = CVtot_t(k,i,j) + CV_WATER * dq / dt
          endif
       enddo
```

### 2.7 Emoist "after" and RHOE_t (2305–2352)

```fortran
          if (.not. l_no_ice_heat) then
            do k = KS, KE
               ! vapor difference
               Emoist(k,2) = - LHV0 * qvv(k) * moist_denv(k) 
               ! ice difference
               do ibi = 1, nbi
                  Emoist(k,2) = Emoist(k,2) &
                        + LHF0 * ( qipv(imt_q,ibi,1,k) - qipv(imw_q,ibi,1,k) - qipv(imat_q,ibi,1,k) ) * moist_denv(k)
               enddo
            enddo
          else
            do k = KS, KE
               ! vapor difference
               Emoist(k,2) = - LHV0 * qvv(k) * moist_denv(k) 
               ! (liquid_heat/ice_heat alternative is fully commented out, lines 2325-2340)
            enddo
          endif

          ! diabatic heating tendency, potential energy rho g h to be calculated in sedimentation later
          do k = KS, KE
             RHOE_t(k,i,j) = ( Emoist(k,2) - Emoist(k,1) ) / dt
          enddo
```

---

## 3. `moistthermo2_scale` (interface file 964–1170) and `moistthermo2` (`/Users/jcanton/projects/scale_amps/contrib/AMPS/mod_amps_utility.F90` 5354–5646)

### 3.1 `moistthermo2_scale` — thresholds (1029–1046)

```fortran
    ! mass mixing ratio limit for cloud droplets
    !  rad=0.5 micron and 1 1/m^3
    real(DP), parameter :: RCLMT=0.5d-15
    ! mass mixing ratio limit for rain
    !  rad=100 micron and 0.025 1/m^3
    real(DP), parameter :: RRLMT=1.0d-10
    ! mass mxing ratio limit for ice
    !  rad=0.5 micron and 1 1/m^3
    real(DP), parameter :: RILMT=1.0d-15
    ! for bin
    real(DP), parameter :: RRLMTB=1.0d-22
    real(DP), parameter :: RILMTB=1.0d-22

    cpr  = CPdry / Rdry
    aklv = LHV0 / CPdry
    akiv = LHS0 / CPdry
```

### 3.2 qc/qr/qi partition + cloudy-mask + qtp/thil, verbatim (1053–1117)

```fortran
      ibr_st=nbhzcl+1
      nwv=2
      do k = KS-1, KE

        qr(k)=0.0_RP
        qi(k)=0.0_RP
        qc(k)=0.0_RP
        micptr(k)=0
        cnr=0.0_RP
        qtotal = 0.0_RP
        do icr=1,ncr
          do ibr=1,ibr_st-1
            qtotal = qtotal + qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k)
            if(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
              qc(k)=qc(k)   &
                    +max(0.0_RP,qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k))
              micptr(k)=1
              cnr=cnr+qrp(rcon_q,ibr,icr,k)
            else
              do ipr=1,npr
                qrp(ipr,ibr,icr,k)=0.0_RP
              end do
            end if
          end do

          do ibr=ibr_st,nbr
            qtotal = qtotal + qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k)
            if(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
              qr(k)=qr(k) &
                   +max(0.0_RP,qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k))
              micptr(k)=1
              cnr=cnr+qrp(rcon_q,ibr,icr,k)
            else
              do ipr=1,npr
                qrp(ipr,ibr,icr,k)=0.0_RP
              end do
            endif
          end do
        enddo

        do ibi=1,nbi
          do ici=1,nci
            qtotal = qtotal + qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k)
            if(qip(imt_q,ibi,ici,k).ge.RILMTB) then
              qi(k)=qi(k) &
                   +max(0.0_RP,qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k)-qip(imw_q,ibi,ici,k))

              qr(k)=qr(k) &
                   +max(0.0_RP,qip(imw_q,ibi,ici,k))
              micptr(k)=1
            else
              do ipi=1,npi
                qip(ipi,ibi,ici,k)=0.0_RP
              end do
            endif
          enddo
        enddo

        ! calculate the total water
        qtotal = qtotal + qv(k)
        qtp(k) = qtotal

        ! calculate the theta il
        thp(k)=th(k)/(1.0_RP  +(aklv*(qr(k)+qc(k))   &
                     +akiv*qi(k))/max(t(k),253.0_RP))
```

**Cloudy-mask criteria (micptr):** any liquid bin with `qrp(rmt_q) ≥ RRLMTB=1e-22` or ice bin with `qip(imt_q) ≥ RILMTB=1e-22` sets `micptr(k)=1`; bins below the threshold are **zeroed in-place** (all npr/npi properties). Sub-threshold bin mass still contributes to `qtotal`/`qtp`.

### 3.3 Initial guess + supersaturation trigger (1119–1164)

```fortran
!------ This is initial guess ------
        th1=th(k)
        ov=qv(k)
        !deb_t = t(k)
        if(ivis<=1) then
          qv(k)=max(0.0_RP,qtp(k)-qr(k)-qi(k)-qc(k))
          !qv(k)=max(0.D0,qtotal-qr(k)-qi(k)-qc(k))

          til=thp(k)*pi(k)/CPdry
          t(k)=til*(1.0_RP+  &
               (aklv*(qr(k)+qc(k))+akiv*qi(k))/253.0_RP)
          if(t(k)>253.0_RP) then
            t(k)=0.5_RP*(til+sqrt(til**2+4.0_RP*til*  &
                 (aklv*(qr(k)+qc(k))+akiv*qi(k))))
            if(t(k)<253.0_RP-1.0e-3_RP) then
               LOG_ERROR("moistthermo2_scale",*) "mo2 t(k) is not right",t(k)
              call PRC_abort
            endif
          endif
          th1=thp(k)*(1.0_RP+(aklv*(qr(k)+qc(k)) &
              +akiv*qi(k))/max(t(k),253.0_RP))

          th(k)=th1
          theta=th1*(EPSvap+qv(k))/EPSvap/(1.0_RP+qv(k))

          es = get_sat_vapor_pres_lk(1, t(k), estbar, esitbar )/10.0_RP
          qs = es/(ptot(k) - es)*0.622_RP
          sw=ptot(k)*qv(k) / ( 0.622_RP * ( 1.0_RP+0.61_RP*qv(k) ) ) / es - 1.0_RP
          if(sw.gt.0.50_RP) then
             LOG_INFO("moistthermo2_scale",'("mo2 sup:",a,4I5,20ES15.6)')  &
                  from,isect,k,i,j, &
                  ptot(k),qv(k),qr(k)+qc(k),qi(k),t(k), &
                  thp(k),pi(k), &
                  sw
          endif

          if(ivis==1) then
            if(qv(k).gt.qs) micptr(k)=1
          end if
        end if

        if(ivis>=1.and.t(k).lt.273.16_RP) then
          es = get_sat_vapor_pres_lk(2, t(k), estbar, esitbar )/10.0_RP
          qvis = es/(ptot(k) - es)*0.622_RP
          if(qv(k).gt.qvis) micptr(k)=1
        end if
```

Note the `/10.0_RP` — `get_sat_vapor_pres_lk` returns CGS (dyn cm⁻² = 0.1 Pa), divided by 10 to compare against SI `ptot`.

### 3.4 AMPS-native `moistthermo2` — key differences (mod_amps_utility.F90 5418–5637)

Same bin partition/mask (`RRLMTB=1.0d-22`, `RILMTB=1.0d-22`, 5450–5493, incl. the same imw split: `qi += imt-imat-imw`, `qr += imw`) **without** the qtotal/qtp accumulation; instead it has an integrity rescale absent from the `_scale` version (5497–5524):

```fortran
            if(qr(k)+qi(k).gt.qtp(k)) then
               qtp(k)=max(qtp(k),0.0_RP)
               if( qc(k)+qr(k)+qi(k) .gt. 1.0e-20_RP ) then
                  fct=max(0.0_RP,qtp(k)/(qc(k)+qr(k)+qi(k)))
               else
                  qtp(k)=0.0_RP
                  fct=0.0_RP
                  nwv=0
               end if
               qc(k)=qc(k)*fct
               qr(k)=qr(k)*fct
               qi(k)=qi(k)*fct
               ...(qrp/qip all properties 1..np*-nwv scaled by fct)...
            endif
```

and diagnoses qv from a prognostic qtp (5528: `qv(k)=max(0.0_RP,qtp(k)-qr(k)-qi(k)-qc(k))`), pressure/T from Exner (`pi=(pi0(k)+pb(k))/cpd; ptot(k)=pi**cpr*1.e5_RP; t(k)=thp(k)*pi`, 5530–5532, guarded by `if(ivis.eq.0)`), uses table lookups `estbar/esitbar` directly (no /10; those tables are in the model's pressure unit), an `ivis.eq.2` saturation-adjustment iteration (5579–5607) with `if(qc(k).gt.0.0_RP) micptr(k)=1`, and the same ice-supersaturation trigger (5631–5637):

```fortran
            if(ivis>=1.and.t(k).lt.273.16_RP) then
               J=MAX(1,MIN(int(T(K))-163,110))
               wt=MAX(MIN(T(K)-(J+163),1.0_RP),0.0_RP)
               esi=esitbar(J)*(1.-wt)+esitbar(J+1)*wt
               qvis=0.622_RP*esi/max((ptot(k)-0.378_RP*esi),esi)
               if(qv(k).gt.qvis) micptr(k)=1
            end if
```

In `_scale`, density is **not** re-diagnosed (`!den(k)=ptot(k)/(r*theta(k)*pi) ! CHIARUI: DENSITY STAYS CONSTANT IN SCALE`, line 5557 of native version comment; `_scale` has no den update at all).

---

## 4. `ini_cloud_micro` + `return_output` — `/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Cloud_Micro.F90`

### 4.1 `ini_cloud_micro` (647–752) — verbatim core

```fortran
  subroutine ini_cloud_micro( &
       CM &
      ,XC &
      ,XR, NRTYPE,NRBIN, NRCAT  &
      ,XS, NSTYPE,NSBIN, NSCAT  &
      ,XA, NATYPE,NABIN, NACAT  &
      ,RV, DEN, PT, T, W  &
      ,L,qtp,ID,JD,KD)
```

Documented type layout (673–683):

```fortran
    ! ITYPE : 1. mixing ratio (g/g)
    !         2. concentration (#/g)
    !         3. volume of circumscribing sphere (cm^3/g)
    !         4. mixing ratio of a-axis production (g/g)
    !         5. mixing ratio of c-axis production (g/g)
    !         6. mixing ratio of d-axis production (g/g)
    !         7. mixing ratio of rossetta production (g/g)
    !         8. mixing ratio of riming production (g/g)
    !         9. volume of riming production (cm^3/g)
    !        10. mixing ratio of melt water (g/g)
```

Body (708–736, comments elided):

```fortran
    CM%L=L
    CM%rain%L=L
    CM%solid_hydro%L=L
    CM%aerosol(1:CM%ncat_a)%L=L

    CM%mes_rc=1

    ! initialize thermo_var object
    call ini_AirGroup(CM%air,L,RV,DEN,PT,T,W,CM%rdsd,CM%ihabit_gm_random)

    ! initialize the concentration and mass in each bin of all the groups
    call ini_group_all(CM%level_comp,CM%air,CM%rain,CM%solid_hydro,&
         CM%aerosol,CM%fcon_c,&
         L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT,&
         XC,XR,XS,XA,CM%flagp_c,CM%flagp_r,CM%flagp_s,CM%flagp_a,CM%coef_ap,&
         CM%eps_ap0,CM%den_apt0,CM%den_aps0,CM%den_api0,ID,JD,KD)
```

### 4.2 `return_output` (1427–1455) — verbatim

```fortran
  subroutine return_output (CM,L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       NATYPE,NABIN,NACAT,RV,XC,XR,XS,XA,qtp)
    use class_Group, only: &
       update_modelvars_all, &
       cal_needgive
    use class_Thermo_Var, only: &
       renew_rv_var
    type (Cloud_Micro), intent(in) :: CM
    integer,intent(in) :: L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT
    real(MP_KIND)    :: RV(*),XC(*),XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),&
         XA(NATYPE,NABIN,NACAT,*),qtp(*)
    integer                     :: i

    call update_modelvars_all(XC,XR,XS,XA,RV, &
             L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT,&
             CM%level_comp,CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
             CM%flagp_r,CM%flagp_s,CM%flagp_a)
  end subroutine return_output
```

### 4.3 The reverse packing — `update_modelvars_all` (`/Users/jcanton/projects/scale_amps/contrib/AMPS/class_Group.F90` 3866–4106), key verbatim

Vapor / rain (3896–3930):

```fortran
    ! +++ vapor +++
    do n=1,L
      rv(n)=ag%tv(n)%rv
    enddo
    !
    ! +++ for rain +++
    !
    if( flagp_r > 0 ) then
       do n = 1, L
       do i = 1, gr%N_BIN
          XR(rmt_q,i,nrcat,n)=gr%MS(i,n)%mass(rmt)/ag%tv(n)%den
          XR(rcon_q,i,nrcat,n)=gr%MS(i,n)%con/ag%tv(n)%den
       enddo
       enddo
       if(level>=4) then
          do n = 1, L
          do i = 1, gr%N_BIN
             XR(rmat_q,i,nrcat,n)=gr%MS(i,n)%mass(rmat)/ag%tv(n)%den
             XR(rmas_q,i,nrcat,n)=gr%MS(i,n)%mass(rmas)/ag%tv(n)%den
          enddo
          enddo
       endif
       ! for cloud droplets
       do n=1,L
          XC(n)=XR(rmt_q,1,nrcat,n)
       enddo
    endif
```

Ice masses (3939–3951) and non-mass cubed-length packing (3953–4035):

```fortran
       do n = 1, L
       do i = 1, gs%N_BIN
        XS(imt_q,i,nscat,n)=gs%MS(i,n)%mass(imt)/ag%tv(n)%den
        XS(icon_q,i,nscat,n)=gs%MS(i,n)%con/ag%tv(n)%den
        ! mass by ice crystals (g/g)
        XS(imc_q,i,nscat,n)=gs%MS(i,n)%mass(imc)/ag%tv(n)%den
        ! mass by riming process (g/g)
        XS(imr_q,i,nscat,n)=gs%MS(i,n)%mass(imr)/ag%tv(n)%den
        ! mass by agg process (g/g)
        XS(ima_q,i,nscat,n)=gs%MS(i,n)%mass(ima)/ag%tv(n)%den
      enddo
      enddo

      iq1=ivcs_q
      nonmass_loop: do k=1,gs%N_nonmass
        !      1. volume of circumscribing sphere (cm^3/g)
        !      2. con. weighted a-axis length^3 (cm^3/g)
        !      3. con. weighted c-axis length^3 (cm^3/g)
        !      4. con. weighted d-axis length^3 (cm^3/g)
        !      5. con. weighted r-axis length^3 (cm^3/g) (rosetta bulletes)
        !      6. con. weighted e-axis length^3 (cm^3/g) (polycrystals)

        if(k==ivcs) then
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%v_cs
        elseif(k==iacr) then
            aq(i,n) = gs%MS(i,n)%con*gs%MS(i,n)%a_len**3.0
        elseif(k==iccr) then
            aq(i,n) = gs%MS(i,n)%con*gs%MS(i,n)%c_len**3.0
        elseif(k==idcr) then
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%d**3.0
        elseif(k==iag) then
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%ag**3.0
        elseif(k==icg) then
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%cg**3.0
        elseif(k==inex) then
            aq(i,n) = gs%MS(i,n)%con*gs%IS(i,n)%n_exice
        endif
        ! (each branch is wrapped in "do n = 1, L / do i = 1, gs%N_BIN"; elided)
        do n = 1, L
        do i = 1, gs%N_BIN
          XS(iq1-1+k,i,nscat,n)=aq(i,n)/ag%tv(n)%den
        enddo
        enddo
      enddo nonmass_loop

      if(level>=6) then
          ! melt water mass  (g/g)
          XS(imw_q,i,nscat,n)=gs%MS(i,n)%mass(imw)/ag%tv(n)%den
          ! freezing nucleation mass  (g/g)
          XS(imf_q,i,nscat,n)=gs%MS(i,n)%mass(imf)/ag%tv(n)%den
      endif
      if(level>=4) then
          ! mass by total mass of aerosols (g/g)
          XS(imat_q,i,nscat,n)=gs%MS(i,n)%mass(imat)/ag%tv(n)%den
          ! mass by soluble mass of aerosols (g/g)
          XS(imas_q,i,nscat,n)=gs%MS(i,n)%mass(imas)/ag%tv(n)%den
      endif
```

Aerosol (4094–4102):

```fortran
        do n = 1, L
        do i = 1, ga(ic)%N_BIN
          XA(acon_q,i,ic,n)=ga(ic)%MS(i,n)%con/ag%tv(n)%den
          ! total mass of aerosols (g/g)
          XA(amt_q,i,ic,n)=ga(ic)%MS(i,n)%mass(amt)/ag%tv(n)%den
          ! mass by soluble mass of aerosols (g/g)
          XA(ams_q,i,ic,n)=ga(ic)%MS(i,n)%mass(ams)/ag%tv(n)%den
        enddo
        enddo
```

**Cubed-length convention:** the qipv/XS pseudo-axis tracers are *concentration-weighted cubes*: `XS(iacr_q) = con·a_len³/den`, `iccr_q = con·c_len³/den`, `idcr_q = con·d³/den`, `iag_q = con·ag³/den`, `icg_q = con·cg³/den`, `inex_q = con·n_exice/den` (no cube), `ivcs_q = con·V_cs/den`. Recovery is `len = (X*_q/icon_q)^(1/3)` (see §5.1).

---

## 5. `ini_group_all` (`class_Group.F90` 595–1682) and `ini_AirGroup` (`class_AirGroup.F90` 156–262)

### 5.1 `ini_group_all` — XR/XS/XA → Group essentials (verbatim, debug/comments elided)

Liquid mass+con (753–768):

```fortran
       do j = 1, gr%L
       do i = 1, gr%N_BIN
          gr%MS(i,j)%mass(rmt)=XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den
          if(XR(rmt_q,i,NRCAT,j)*ag%TV(j)%den>0.0) then
             gr%MS(i,j)%con=XR(rcon_q,i,NRCAT,j)*ag%TV(j)%den
          else
             gr%MS(i,j)%con=0.0_PS
             ...
          end if
       enddo
       enddo
```

Liquid aerosol masses, level≥4 (777–825, key lines):

```fortran
            if(XR(rmat_q,i,NRCAT,j)>0.0) then
              ! total aerosol mass
              gr%MS(i,j)%mass(rmat)=min( &
                     max(m_lmt_ap,min_fapt_r*gr%MS(i,j)%mass(rmt) &
                    ,XR(rmat_q,i,NRCAT,j)*ag%TV(j)%den),mx_aprat*gr%MS(i,j)%mass(rmt) )
            else
              gr%MS(i,j)%mass(rmat)=min(max(m_lmt_ap,min_fapt_r*gr%MS(i,j)%mass(rmt)),&
                     gr%MS(i,j)%mass(rmt))
            end if
...
            if(XR(rmas_q,i,NRCAT,j)>0.0) then
              ! soluble mass of aerosols
              gr%MS(i,j)%mass(rmas)=min(gr%MS(i,j)%mass(rmat),max(m_lmt_ap,&
                     max(&
                     min_faps_r*gr%MS(i,j)%mass(rmat)&
                    ,min(XR(rmas_q,i,NRCAT,j)*ag%TV(j)%den,gr%MS(i,j)%mass(rmat))) &
                       ))
            else
              gr%MS(i,j)%mass(rmas)=min(max(m_lmt_ap,min_faps_r*gr%MS(i,j)%mass(rmat)),&
                     gr%MS(i,j)%mass(rmat))
            end if
            ! insoluble mass is diagnosed in diag_pq
            gr%MS(i,j)%mass(rmai)=0.0_PS
          else   ! con <= 0
            gr%MS(i,j)%mass(rmat)=0.0_PS ; gr%MS(i,j)%mass(rmas)=0.0_PS ; gr%MS(i,j)%mass(rmai)=0.0_PS
```

Ice mass+con (885–896) — identical pattern (`mass(imt)=XS(imt_q,...)*den`, `con=XS(icon_q,...)*den` if mass>0 else 0). Ice aerosol masses `imat/imas` (909–960) mirror the rain formulas with `min_fapt_s/min_faps_s/mx_aprat`. Then mass components with consistency clamps (967–1115, key lines):

```fortran
          ! crystal mass component
          if( XS(imt_q,i,NSCAT,j) <= XS(imc_q,i,NSCAT,j)) then
            gs%MS(i,j)%mass(imc) = gs%MS(i,j)%mass(imt)
          else
            if(gs%MS(i,j)%mass(imt)>=gs%MS(i,j)%con*m_icmin) then
              gs%MS(i,j)%mass(imc)=max( &
                   min(XS(imc_q,i,NSCAT,j)*ag%TV(j)%den,gs%MS(i,j)%mass(imt)) &
                  ,gs%MS(i,j)%con*m_icmin)
            else
              gs%MS(i,j)%mass(imc) = min(XS(imc_q,i,NSCAT,j)*ag%TV(j)%den,gs%MS(i,j)%mass(imt))
            end if
          end if
...
          ! rime mass component
          if( gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc) < XS(imr_q,i,NSCAT,j)*ag%TV(j)%den ) then
            gs%MS(i,j)%mass(imr) = max(gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc),0.0_PS)
          else
            gs%MS(i,j)%mass(imr) = max(XS(imr_q,i,NSCAT,j)*ag%TV(j)%den,0.0_PS)
          end if
          ! aggregation mass component
          if( gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc)<XS(ima_q,i,NSCAT,j)*ag%TV(j)%den ) then
            gs%MS(i,j)%mass(ima) = max(gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc),0.0_PS)
          else
            gs%MS(i,j)%mass(ima) = max(XS(ima_q,i,NSCAT,j)*ag%TV(j)%den,0.0_PS)
          end if
...(level>=6: imw and imf analogous; imf clamped by mass(imc))...
          ! maintain the relative relation of predicted mass component
          mass_left(i,j)=max(0.0_PS,gs%MS(i,j)%mass(imt)-gs%MS(i,j)%mass(imc))
          mass_stot(i,j)=gs%MS(i,j)%mass(imr)+gs%MS(i,j)%mass(ima)+gs%MS(i,j)%mass(imw)
          if(mass_stot(i,j)>1.0e-25_PS) then
            rat=mass_left(i,j)/mass_stot(i,j)
            gs%MS(i,j)%mass(imr)=gs%MS(i,j)%mass(imr)*rat
            gs%MS(i,j)%mass(ima)=gs%MS(i,j)%mass(ima)*rat
            gs%MS(i,j)%mass(imw)=gs%MS(i,j)%mass(imw)*rat
          else
            gs%MS(i,j)%mass(imc)=gs%MS(i,j)%mass(imt)
            gs%MS(i,j)%mass(imr)=0.0_PS ; ...(ima)=0.0_PS ; ...(imw)=0.0_PS
          end if
```

Axis recovery — the cube-root convention (1120–1232, verbatim key lines):

```fortran
          if( XS(iacr_q,i,NSCAT,j) <= 0.0 .or. XS(iccr_q,i,NSCAT,j) <= 0.0) then
            if(gs%MS(i,j)%mass(imt)<=1.0e-30.or.gs%MS(i,j)%con<=1.0e-30) then
              gs%MS(i,j)%mass(imt)=0.0
              gs%MS(i,j)%con=0.0
            else
              call gen_length(gs%MS(i,j)%mass(imc)/gs%MS(i,j)%con,& !ag%TV(j)%T,
                      gs%MS(i,j)%a_len,gs%MS(i,j)%c_len)
              gs%IS(i,j)%d=0.0_PS ; gs%IS(i,j)%ag=0.0_PS ; gs%IS(i,j)%cg=0.0_PS ; gs%IS(i,j)%n_exice=0.0_PS
            endif
          else
            gs%MS(i,j)%a_len=min(max(1.0e-4_PS,&
                   (XS(iacr_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0)),5.0_PS)
            gs%MS(i,j)%c_len=min(max(1.0e-4_PS,&
                   (XS(iccr_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0)),5.0_PS)
          endif
...
            ! reality check of axis ratio  (phi_max=2.0d+1, phi_min=5.0d-3)
            if(gs%MS(i,j)%c_len>gs%MS(i,j)%a_len*phi_max) then
              gs%MS(i,j)%c_len=gs%MS(i,j)%a_len*phi_max
            elseif(gs%MS(i,j)%c_len<gs%MS(i,j)%a_len*phi_min) then
              gs%MS(i,j)%a_len=gs%MS(i,j)%c_len/phi_min
            end if
            ! d axis lengths of ice crystals should be bounded by a-axis length
            if( level >= 2 ) then
              if( XS(idcr_q,i,NSCAT,j) > 0.0.and.gs%MS(i,j)%a_len>=10.0e-4_PS ) then
                gs%IS(i,j)%d = max(min( &
                        (XS(idcr_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0) &
                        ,0.9_PS*gs%MS(i,j)%a_len),0.0_PS)
              else
                gs%IS(i,j)%d = 0.0_PS
              end if
            end if
            ! number of extra single ice crystals
            if( level==7 ) then
              gs%IS(i,j)%n_exice=max(0.0_PS,min(max_exice,XS(inex_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)))
            end if
            ! position of center of gravity: a-axis coordinate
            if( level==3.or.level==5.or. level==7 ) then
              if( XS(iag_q,i,NSCAT,j) > 0.0 ) then
                gs%IS(i,j)%ag=max(min( &
                     (XS(iag_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0/3.0) &
                     ,gs%MS(i,j)%a_len),0.0_PS)
              else
                gs%IS(i,j)%ag=0.0_PS
              end if
            end if
            ! position of center of gravity: c-axis coordinate
            if( level==3.or.level==5.or.level==7 ) then
              if( XS(icg_q,i,NSCAT,j) > 0.0_PS ) then
                gs%IS(i,j)%cg=max(min( &
                     (XS(icg_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j))**(1.0_PS/3.0_PS) &
                     ,gs%MS(i,j)%c_len),0.0_PS)
              else
                gs%IS(i,j)%cg=0.0_PS
              end if
            end if
```

plus a hexagonal-plate density reality check (1237–1271) that can rescale all five lengths from `a_min=(mass(imc)/con/(3*sq_three*phi*(1-psi)*den_i))**0.33333333`. V_cs (1280–1347):

```fortran
          if(XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j) <= 0.0_PS ) then
            ! ... zero the entire bin (con, all masses, V_cs, all lengths)
          elseif(XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)<V_csmin_hex) then
            gs%IS(i,j)%V_cs=V_csmin_hex
            gs%MS(i,j)%a_len = 1.0e-4_PS
            gs%MS(i,j)%c_len = 1.0e-4_PS
            ! ... reset bin to minimum crystal: mass(imc)=mass(imt)=m_icmin*con, aerosol masses re-derived
          elseif(XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)>V_csmax) then
            gs%IS(i,j)%V_cs=V_csmax
          else
            gs%IS(i,j)%V_cs=XS(ivcs_q,i,NSCAT,j)/XS(icon_q,i,NSCAT,j)
          end if
```

Aerosol groups — `flagp_a==1` (1382–1407) diagnostic-mass mode:

```fortran
          con1=XA(acon_q,i,k,j)*ag%TV(j)%den
          if(con1<0.0_PS) then
            ga(k)%MS(i,j)%con=0.0_PS ; ga(k)%MS(i,j)%mass(amt)=0.0_PS ; ga(k)%MS(i,j)%mass(ams)=0.0_PS
          else
            ga(k)%MS(i,j)%con=con1
            tmass1=con1*coef_ap(k)
            ga(k)%MS(i,j)%mass(amt)=tmass1
            ga(k)%MS(i,j)%mass(ams)=tmass1*ga(k)%MS(i,j)%eps_map
          end if
```

else-branch (predicted mass, 1409–1474):

```fortran
        m_lmt=coef4pi3*r3_lmt*den_apt0(k)
        ...
          ga(k)%MS(i,j)%con=max(n_lmt_ap,XA(acon_q,i,k,j)*ag%TV(j)%den)
          tmass1=XA(amt_q,i,k,j)*ag%TV(j)%den
          ga(k)%MS(i,j)%mass(amt)=tmass1
          if(level>=4) then
            tmass2=max(0.0_PS,min(XA(ams_q,i,k,j)*ag%TV(j)%den,tmass1))
            ga(k)%MS(i,j)%mass(ams)=tmass2
            ! k==1/3/4 (CCN): min(max(tmass1*sep_faps,tmass2),tmass1)
            ! k==2 (IN):      max(min(tmass1*0.99*sep_faps,tmass2),0.0_PS)
          end if
        ...
          tmass1=m_lmt*ga(k)%MS(i,j)%con
          if(ga(k)%MS(i,j)%mass(amt)<tmass1) then
            ga(k)%MS(i,j)%mass(amt)=tmass1
            ga(k)%MS(i,j)%mass(ams)=eps_ap0(k)*tmass1
          end if
```

Everything not set from XR/XS/XA is initialized by `initialize_ig`/`initialize_shape_ig` (1504–1680): all `dcondt/dmassdt/dvoldt/Ldmassdt = 0.0_PS`, all state fields (`mass`, `vol`, `con`, `mean_mass`, `len`, `den`, `a_len`, `c_len`, coefficients, `vtm`, ..., all `IS` shape fields) = `undef`, except `den_as=den_aps0(ica)`, `den_ai=den_api0(ica)`, `eps_map=eps_ap0(ica)`, `mark/inevp/inmlt=0`, `mark_cm/mark_er=0`.

### 5.2 `ini_AirGroup` — the SI→CGS conversion point (class_AirGroup.F90 191–207, verbatim)

```fortran
    do n=1,a%L
      a%TV(n)%rv = rv(n)
      a%TV(n)%den = den(n)*1.0e-3_PS      ! kg m-3 -> g cm-3
      a%TV(n)%P = PT(n)*1.0e+1_PS         ! Pa -> dyn cm-2
      a%TV(n)%T = T(n)
      a%TV(n)%W = W(n)*1.0e+2_PS          ! m s-1 -> cm s-1

      a%TV(n)%D_v = get_diffusivity(a%TV(n)%P,a%TV(n)%T)
      a%TV(n)%k_a = get_thermal_conductivity(a%TV(n)%T)
      a%TV(n)%d_vis = get_dynamic_viscosity( a%TV(n)%T)
      a%TV(n)%sig_wa=get_sfc_tension(a%TV(n)%T)
      a%TV(n)%e=a%TV(n)%P*a%TV(n)%rv/(Rdvchiarui+a%TV(n)%rv)

      a%TV(n)%den_a=M_a*(a%TV(n)%P-a%TV(n)%e)/(R_u*a%TV(n)%T)

      a%TV(n)%dmassdt=0.0_PS
```

(unit-conversion comments mine; remainder 211–260 derives `e_sat(1:2)` via `get_sat_vapor_pres_lk`, `s_v(1:2)`, `GTP(1:2)`, `rv_sat(1:2)=Rdvchiarui*e_sat/max(P-e_sat,e_sat)`, clamps `s_v` sign against `rv_sat` vs `rv`, zeroes `dmassdt_v(1:3)`, sets `T_n=T_m=T`, caps `e_sat(2)/rv_sat(2)` at their liquid values, copies `s_v_n/e_sat_n`, and sets `nuc_gmode` via `cal_growth_mode_inl_vec`.)

Density anywhere inside AMPS Groups is thus **g cm⁻³** (`ag%TV%den`), which is why `ini_group_all`/`update_modelvars_all` multiply/divide by `ag%TV(j)%den` without further factors, and why the SCALE interface applies the lone `0.001` only to non-mass tracers (which are per-cm³ in QTRC but per-(g cm⁻³ of moist air) in `q*pv`).

---

## 6. What M1's state bundles must hold

To reproduce a round trip (pack → AMPS → unpack) exactly, a per-column state bundle must capture, per (k,i,j):

**Pre-call SCALE state (inputs to both pack and tendency-differencing):** `DENS`, `QDRY` (⇒ `factor_mxr1/2`), `PRES`, `TEMP`, `W`, `MOMZ`, `CZ(KS-1:KS)`, `dt`, and the full `QTRC` slices: `I_QV`, `I_QL+0..nbr-1`, `I_QI+0..nbi-1`, `I_QW+0..nbi-1`, `I_QPPVL+0..3*nbr-1` (per bin: rmat, rmas, rcon), `I_QPPVI` block (per bin: rime..crystal, frozen, imat, imas, then non-mass icon/ivcs/iacr/iccr/idcr/iag/icg/inex), `I_QPPVA+0..3*nba*nca-1` (amt, acon, ams). Also entering tendencies `DENS_t` (accumulated) and namelist switches `l_gaxis_version`, `l_axis_limit`, `l_no_ice_heat`, `l_sediment`, `level`, `nbhzcl`.

**AMPS-side packed state (what `ini_cloud_micro` consumes, all CGS-facing, per moist air):** `qvv (rv)`, `moist_denv (den, SI — converted ×1e-3 inside)`, `ptotv (Pa — ×10 inside)`, `tv`, `wbv (×100 inside)`, plus `qrpv(npr,nbr,ncr)`, `qipv(npi,nbi,nci)`, `qapv(npa,nba,nca)` with the last two property slots zeroed, and the diagnostics `qtp`, `thp/thv/thetav/piv/pbv`, `micptrv`, `k1r/k2r/k1i/k2i`, `qrov/qiov`, `mmassrv/mmassiv`.

**Post-call state needed for the tendencies:** the updated `qvv`, `qrpv`, `qipv`, `qapv`, `moist_denv`, `den_t` (sedimentation), and `Emoist(:,1)`/`Emoist(:,2)` (and `liquid_heat/ice_heat` slot 1 when `l_no_ice_heat`) — outputs are `RHOQ_t` (all above tracer indices), `RHOE_t`, `CPtot_t`, `CVtot_t`, `DENS_t`/`DENS_NEW`, `QLIQ`, `QICE(:,1:6)`.

**Invariants a bundle must preserve** (used implicitly by the recipes): drop mass bookkeeping `I_QL = rmt_q − rmat_q`, ice bookkeeping `I_QI = imt_q − imat_q − imw_q`, `I_QW = imw_q`; the single `0.001` factor exactly on {`rcon_q`; ice `imas_q+1..npi-2`; `acon_q`}; and the concentration-weighted cubed-length convention `X = con·len³/den` with recovery `len = (X/Xcon)^(1/3)` for a/c/d/ag/cg axes (v_cs and n_exice weighted but not cubed).