All evidence gathered. Here is the report.

# Cloudlab temperature range & ice-entry report (branch `cloudlab_port`)

## (a) Domain temperature — initial state + drift bound

**Domain geometry:** `KMAX=60`, `DZ=20 m` → model top = 1200 m; cell centers span ~10–1190 m (`run.conf:20,34`).

**Initial T (from `env.txt`, column 3 = T in K), over 0–1200 m:**
- **Min T = 266.09 K (−7.06 °C) at z ≈ 600 m** (`env.txt`, the T minimum of the profile).
- **Max T = 272.16 K (−1.0 °C) at z ≈ 940 m.**
- Profile is non-monotonic: cools from surface (268.8 K at 200 m) to a 266.1 K minimum near 600 m, then warms to a ~272 K layer at 900–1000 m (cloud/inversion), then cools slightly to 271.4 K at the 1200 m top.
- **The entire domain is below 273.15 K** — a fully supercooled / mixed-phase boundary layer, roughly −7 °C to −1 °C initially.

**Plausible T drift over the 1 h spin-up + 20 min seeding (≈1.3 h, up to 2 h):**
- Radiation is ON (`ATMOS_PHY_RD_TYPE="MSTRNX"`, `run.conf:75`), stepped every 3 s (`TIME_DT_ATMOS_PHY_RD=3.0`, `run.conf:51`). Cloud-top longwave cooling is the dominant thermodynamic driver.
- Surface fluxes are small and constant (`ATMOS_PHY_SF_TYPE="CONST"`, `run.conf:73`): sensible SH = 7.98 W/m², latent LH = 2.86 W/m² (`run.conf:97–98,105–106`), surface temp fixed at 270.5 K (`run.conf:83`). These add only weak near-surface warming.
- A prescribed large-scale subsidence forcing (`largeScaledForce.txt`, applied in `mod_user.F90`) plus SMAGORINSKY turbulence redistribute heat and largely offset cooling in a well-mixed layer.
- **Estimate:** a cloud-top LW flux divergence of ~50–80 W/m² concentrated over a ~50–100 m top layer gives dT/dt ≈ ΔF/(ρ·cp·Δz) ≈ 1–3 K/hr, localized near cloud top. Over ~1.3 h that is ~1.5–4 K of local cooling, tempered by turbulent mixing and subsidence warming.
- **Bound:** domain minimum T realistically stays within a few K of the initial 266 K — i.e. it does **not** drop below ~262 K (−11 °C) nor rise near 0 °C. **The whole domain remains supercooled (below 273.15 K) for the entire run.** No configuration mechanism can push any level to warm (>0 °C) temperatures.

## (b) Seeding injection mechanics (`mod_user.F90`)

**Key finding: seeding injects INP *aerosol*, not ice.** Ice enters the simulation only indirectly — the injected ice-nucleating-particle aerosol is later converted to ice by the AMPS microphysics' heterogeneous nucleation, inside the MP scheme, on a subsequent step. No ice tracer/PPV is written at injection.

Injection loop (`mod_user.F90:490–541`), active only when `DO_CLOUD_SEEDING` and `time_now_second ∈ [TIME_INP_LOWER, TIME_INP_UPPER)`:

- Writes into **aerosol per-particle-variable (PPV) tracers**: index base `QS_MP + I_QPPVA` where `I_QPPVA` = "starting index of aerosol PPVs" (`scale_atmos_phy_mp_amps.F90:491`; `I_QPPVA` imported at `mod_user.F90:439`). It writes **no** liquid (`I_QPPVL`) or ice (`I_QPPVI`) PPV.
- Only aerosol category `ica==2` is seeded (`if ( ica /= 2 ) ... cycle`, line 509); loops all aerosol bins `iba=1..nba`.
- Three quantities per seeded bin (lines 513–518):
  - `QS_MP+I_QPPVA+ipa_qpa-1` ← aerosol **mass** = `RATE * coef_ap(ica) / dt * 1000`
  - `QS_MP+I_QPPVA+ipa_qpa`   ← aerosol **number concentration** = `RATE * DENS / dt`
  - `QS_MP+I_QPPVA+ipa_qpa+1` ← aerosol **soluble mass** = `RATE * coef_ap(ica) * eps_ap(ica) / dt * 1000`
  - `coef_ap`, `eps_ap` are fixed per-category aerosol coefficients from `com_amps` (line 441–442); `RATE = RELEASE_INP_CONC_TIME_RATE`.
- The tendency `RHOQ_t_SEED` is simply added to `RHOQ_t` (lines 529–539). **No call into any AMPS habit-selection or RNG routine occurs at injection** — the injected particles get no habit/shape/axis properties here; they are plain aerosol number/mass/soluble-mass increments. Habit, a/c/d-axis, ice mass etc. are assigned only later, when AMPS nucleates ice on these INP.

**Injection height:** single model level — the loop takes the first `k` with `DOMAIN_CZ(k) >= RELEASE_INP_Z_LOWER_LIMIT + CONST_EPS` and then `exit`s (lines 494–527). With `RELEASE_INP_Z_LOWER_LIMIT = 300.0` (`restart_run.conf`), that is the first cell center above 300 m, **z ≈ 310 m**. (The code comment at line 487 says "500 m per BAMS paper," but the active namelist value is 300 m.)

**Injection horizontal footprint:** one J-row only (`JS`, the first row; the CY < CDY guard, line 496) and X ∈ [949, 1051] m (`RELEASE_INP_X_LOWER/UPPER_LIMIT`, restart).

**Injection time window:** `restart_run.conf` sets HOUR 18→18, MIN 0→1, SEC 0→0 → `[64800 s, 64860 s)` = 18:00:00–18:01:00, i.e. **the first 60 s** of the restart run (which begins at 18:00:00). `RELEASE_INP_CONC_TIME_RATE = 10.0`.

## (c) Other temperature-relevant configuration facts

Spin-up = `run.conf`; seeding = `restart_run.conf`. Notable diffs / settings:

- **Two-phase design (`CLOUDLAB.md:9–17`):** spin-up 1 h with **ice physics off** (no INP); restart run continues from spin-up restart files with **ice physics on, 40 ice bins** (`num_h_bins=40,40` in restart vs `40,20` in spin-up).
- **Durations:** spin-up `TIME_DURATION=3600` starting 2023-01-24 17:00:00 (`run.conf:38,40`); restart `TIME_DURATION=1200` starting 18:00:00 (`restart_run.conf`). So seeding occupies first 60 s of a 20-min run.
- **`RESTART_SKIP_READING_AMPS_ICE=.true.`** in restart — ice state not carried from spin-up (ice starts clean).
- **`MP_do_precipitation=.false.`** in both (`run.conf:137`) — no precip sedimentation removal of hydrometeors from the microphysics driver (AMPS `l_sediment=.true.` still set).
- **Radiation:** MSTRNX with MIPAS climatology profile, fixed lat/lon 47.07 N / 7.87 E, non-fixed date 2023-01-24 (`run.conf:113–133`) — winter, low sun; longwave cloud-top cooling dominates over shortwave.
- **Surface:** CONST fluxes SH=7.98, LH=2.86 W/m², SFC_TEMP=270.5 K, albedo_SW=0.827 (snow) (`run.conf:82–107`).
- **Turbulence:** SMAGORINSKY, vertical only (`ATMOS_PHY_TB_SMG_horizontal=.false.`, `run.conf:110`).
- **`l_no_ice_heat=.false.`** in both configs — latent heating from ice processes IS active in the microphysics (relevant once ice forms after seeding). Recent branch commits (e.g. `85a10f8`, `dffd2e3`) were experiments toggling ice-heat / WBF handling, but the checked-in cloudlab config keeps it enabled.
- **Large-scale subsidence** forcing via `PARAM_USER` (`USER_file="largeScaledForce.txt"`, `SWITCH_RHOT=.true.`, `SWITCH_TEMP=.true.`) applied every step to maintain BL depth — a warming/thermodynamic term counteracting radiative cooling.

**Bottom line:** the cloudlab domain is supercooled top-to-bottom (initial ~266–272 K; realistically staying ~262–273 K through the run), and ice does not enter directly. Seeding deposits INP aerosol (number, mass, soluble mass) into aerosol PPV category 2 at one level ~310 m, over 60 s, with zero habit/RNG involvement at injection; ice appears only when AMPS subsequently nucleates on that INP.