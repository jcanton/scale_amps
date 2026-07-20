# Lead: liquid number-PPV QTRC slot off-by-one (init/diagnostic vs microphysics pack)

Investigation for congchia into why the warm-run reference dumps show the
droplet-number/aerosol liquid PPVs (`rcon_q`/`rmat_q`/`rmas_q`) as uniformly
trace (~1e-9, ~14 orders below physical) while the mass field (`rmt_q`)
validates exactly against the dump's own `qcvm`, even though the LOG shows
realistic clouds. All line numbers are in
`scalelib/src/atmosphere/physics/microphysics/scale_atmos_phy_mp_amps.F90`
unless noted.

## Ruled out (with evidence)

1. **NOT the dump instrumentation.** The gather copies all `ipr` incl. `rcon`
   (L2058-2064) and `AMPS_DUMP_micro` writes `qrpvm` verbatim (L5364), all
   AFTER `moistthermo2_scale` (L1937) and the pack. Mass `rmt_q` reads back
   exactly vs `qcvm`, so the reader/PPV-axis mapping is correct.
2. **NOT the pack/unpack round-trip.** Pack (L1692) divides `rcon` by
   `factor_mxr1/0.001`; the writeback tendency (L2736-2738) multiplies by
   `moist_denv*0.001`. Since `moist_denv = DENS*factor_mxr1` (L1658), a no-op
   microphysics gives `RHOQ_t = QTRC*(moist_denv/factor_mxr1 - DENS)/dt = 0`.
   The number round-trip is exactly conservative — it does not degrade number.

## Leading hypothesis: the number-PPV QTRC slot offset is inconsistent

Two conventions in the same file disagree by ONE slot on where the liquid
number concentration lives within each bin's PPV block (block stride =
`numberPPVL = num_h_moments(1)-1`, L493):

| Code path | number (rcon) QTRC index for bin `ibr` | total-ap (rmat) base | soluble-ap (rmas) base |
|---|---|---|---|
| **Microphysics pack** (`Z_LOOP_01`, L1684-1693) + **unpack/writeback** (L2723-2739) | `I_QPPVL + 2 + numberPPVL*(ibr-1)` | `I_QPPVL+0` | `I_QPPVL+1` |
| **Init** `qhyd2qtrc` (`lcon_index`, L4507) + **diagnostic/output** `qtrc2qhyd` (`liqConc_index = I_QPPVL+2-1`, L3979/L4751) | `I_QPPVL + 1 + numberPPVL*(ibr-1)` | `lamt_index = I_QPPVL-1` | `lams_index = I_QPPVL` |

The pack writes the per-bin QTRC block in the order `[rmat, rmas, rcon]`
(the `ipr_qpr` walk: `rmt` reuses slot 0, then `do ipr=rmat_q,rmas_q` fills
slots 0,1, then `rcon` at slot 2). The canonical qrpv property order that
fixes this is `rmt=1, rmat=2, rmas=3, rcon=4` (`mod_amps_utility.F90:1619-1636`)
— so `rcon` is genuinely the LAST of the three stored PPVs → offset +2.

The init/diagnostic base is shifted -1: `lamt_index = I_QPPVL-1` writes bin-1
total-aerosol mass to the slot BEFORE the liquid-PPV region (the last ice
PPV), and puts number at `I_QPPVL+1`. That `I_QPPVL-1` write is itself a red
flag that the init/diagnostic base is the mis-shifted one.

## Why this reconciles every symptom

- **LOG shows realistic clouds**: the output diagnostic `qtrc2qhyd` reads
  number from `I_QPPVL+1`, where the init actually put it — self-consistent,
  looks fine on paper.
- **Dump shows trace number**: the microphysics pack reads `rcon` from
  `I_QPPVL+2`, which the init left ~0 → `qrpv(rcon)` trace at micro entry →
  dumped trace.
- **Mass real**: bin water mass lives in `I_QL` (L1686), unaffected by the
  PPV offset — validates exactly vs `qcvm`.
- **`rmas > rmt` "everywhere" and `rmat_q` tracks `rcon_q`**: the off-by-one
  scrambles the aerosol/number slots between the two conventions (the pack's
  `rmas` slot = the init's number slot, etc.), exactly the correlated
  inert-residue signature the magnitude diagnostic found.

## Decisive check (congchia)

In one known-cloudy `(k,i,j)`, `LOG` both `QTRC(k,i,j,I_QPPVL+1)` and
`QTRC(k,i,j,I_QPPVL+2)` (for bin 1) next to `qcvm(k)`, right before the L2097
dump. If the real droplet number is in `+1` and `+2` is ~0, the pack/unpack
convention (`+2`) is the one out of step with init/diagnostic (`+1`), and the
two must be reconciled to a single per-bin sub-order. Then also confirm
`num_h_moments(1)` / `numberPPVL` so the stride is unambiguous, and check
whether `lamt_index = I_QPPVL-1` is corrupting the last ice PPV.

## Impact on the port

The icon4py port reads its reference from the DUMP (pack convention, `+2`),
so once the Fortran side is reconciled and re-dumped with real number, the
port's warm per-call validation can proceed. No port change is implied by
this lead — it is a scale_amps-side reconciliation.
