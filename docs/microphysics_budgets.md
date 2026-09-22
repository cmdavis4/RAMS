# Microphysics budget output

This note describes how RAMS outputs the microphysics budget processes on this
branch. It replaces an older convention, and the old variable names no longer
exist.

## Summary

- `IMBUDGET` (0 to 3) chooses which budget processes are computed. That part has
  not changed.
- Each process is written as an **instantaneous rate**, named `<NAME>R`. For
  example, `VAPLIQR` is vapor deposition onto liquid, in kg/kg/s.
- When `IMBUDGETACC = 1`, each process is also written as a **run total**, named
  `ACC<NAME>`. For example, `ACCVAPLIQ` is in kg/kg. The run totals behave like
  `ACCPR`: they start at zero and never reset.
- No budget value depends on `FRQSTATE`, `FRQLITE` or any other output frequency.

`IMBUDGETACC` belongs in the `$MODEL_OPTIONS` namelist next to `IMBUDGET`.
Like every RAMSIN key, it must be present:

```
   IMBUDGET    = 3,   ! Micro budgets: 0=Off,1=partial,2=all,3=all+dust/regen
   IMBUDGETACC = 0,   ! Micro budget run totals (ACC*): 0=off, 1=on
```

## Units

| Kind of process | Rate `<NAME>R` | Run total `ACC<NAME>` |
|---|---|---|
| Mass transfers (`VAPLIQ`, `RIMECLD`, `NUCCLDR`, ...) | kg/kg/s | kg/kg |
| Number transfers (`NUCCLDC`, `REGEN_AERO*_N`, ...) | #/kg/s | #/kg |
| Latent heating (`LATHEATVAP`, `LATHEATFRZ`) | K/s (of theta or temperature, per `LHRTHETA`) | K |

A rate is the amount the process added during the most recent model step,
divided by the step length. It therefore describes the state at the output
time. It is not a mean over the output interval.

## How it works

The bulk microphysics works on one column at a time. Each process adds its
increments to a per-column work vector, such as `xvapliq` (see `micphys.f90`).

1. `range_check()` (`mic_misc.f90`) sets every work vector to zero at the start
   of every step.
2. The microphysics routines add each increment to the work vector. These
   routines have not changed.
3. `copyback()` (`mic_driv.f90`) calls `store_budget_rate()`, which stores
   `x<name> / dtlt` in `micro%<name>r`.
4. When `IMBUDGETACC = 1`, `copyback()` also calls `add_to_budget_total()`, which
   adds `x<name>` to `micro%acc<name>`.

Two other groups of budgets use the same pattern at the grid level:

- **Aerosol regeneration** (`REGEN_AERO*_NR/_MR`) has two producers within one
  step, and the arrays accumulate across both:
  1. `vapflux()` (`mic_vap.f90`) returns aerosol as hydrometeors evaporate,
     whether the evaporation is partial or complete. It fills the per-column
     vectors `xregenaero1n/1m/2n/2m`, which `copyback()` adds to the arrays
     through `add_to_regen_budget()`. This is the larger of the two sources.
  2. `adj1()` (`mic_adj.f90`) restores what is left in a hydrometeor category
     driven below `rxmin` by some route other than evaporation, such as
     collection or melting. It writes the arrays directly. `negadj1` calls
     `adj1` twice per step, before and after the microphysics.

  `timestep()` (`rtimh.f90`) calls `regen_budget_step_start()` before the first
  `negadj1` to zero the arrays. It calls `regen_budget_step_finish()` after the
  second, which adds the step's total to the run totals, divides by `dtlt`, and
  sets the below-ground level from level 2. Both are in `mic_adj.f90`.

  Before September 2026 only the `adj1` producer was counted, so these
  variables reported a small residue rather than the regeneration source.
  Output written before that date understates regeneration.
- **Bin microphysics** (`LEVEL = 4`) adds increments directly to
  `micro%<name>r`. `zero_budgets()` zeroes those arrays each step, and
  `finish_bin_budgets()` converts them after `micro_proc()`. Both routines are
  in `micro_bin/microphysics.f90`.

The run-total arrays exist only when `IMBUDGETACC = 1`. Any code that reads one
must check `allocated()` first.

## Old names and new names

The old `*T` totals were reset after every state-file write, so their values
depended on `FRQSTATE`. The old `LATHEATVAP` and `LATHEATFRZ` were in K per
step. All of them are replaced as follows:

| Old output | Rate (always, when gated on) | Run total (IMBUDGETACC=1) | Gate |
|---|---|---|---|
| `LATHEATVAP` (K/step) and `LATHEATVAPT` | `LATHEATVAPR` | `ACCLATHEATVAP` | IMBUDGET>=1 |
| `LATHEATFRZ` (K/step) and `LATHEATFRZT` | `LATHEATFRZR` | `ACCLATHEATFRZ` | IMBUDGET>=1 |
| `NUCCLDRT` | `NUCCLDRR` | `ACCNUCCLDR` | IMBUDGET>=1 |
| `CLD2RAINT` | `CLD2RAINR` | `ACCCLD2RAIN` | IMBUDGET>=1 |
| `ICE2RAINT` | `ICE2RAINR` | `ACCICE2RAIN` | IMBUDGET>=1 |
| `NUCICERT` | `NUCICERR` | `ACCNUCICER` | IMBUDGET>=1 |
| `VAPLIQT` | `VAPLIQR` | `ACCVAPLIQ` | IMBUDGET>=1 |
| `VAPICET` | `VAPICER` | `ACCVAPICE` | IMBUDGET>=1 |
| `EVAPLIQT` | `EVAPLIQR` | `ACCEVAPLIQ` | IMBUDGET>=1 |
| `EVAPICET` | `EVAPICER` | `ACCEVAPICE` | IMBUDGET>=1 |
| `FREEZINGT` | `FREEZINGR` | `ACCFREEZING` | IMBUDGET>=1 |
| `MELTINGT` | `MELTINGR` | `ACCMELTING` | IMBUDGET>=1 |
| `MELTICET` | `MELTICER` | `ACCMELTICE` | IMBUDGET>=1 |
| `RIMECLDT` | `RIMECLDR` | `ACCRIMECLD` | IMBUDGET>=1 |
| `RAIN2ICET` | `RAIN2ICER` | `ACCRAIN2ICE` | IMBUDGET>=1 |
| `AGGREGATET` | `AGGREGATER` | `ACCAGGREGATE` | IMBUDGET>=1 |
| `INUCHOMRT` | `INUCHOMRR` | `ACCINUCHOMR` | IMBUDGET>=2 |
| `INUCCONTRT` | `INUCCONTRR` | `ACCINUCCONTR` | IMBUDGET>=2 |
| `INUCIFNRT` | `INUCIFNRR` | `ACCINUCIFNR` | IMBUDGET>=2 |
| `INUCHAZRT` | `INUCHAZRR` | `ACCINUCHAZR` | IMBUDGET>=2 |
| `VAPCLDT` | `VAPCLDR` | `ACCVAPCLD` | IMBUDGET>=2 |
| `VAPRAINT` | `VAPRAINR` | `ACCVAPRAIN` | IMBUDGET>=2 |
| `VAPPRIST` | `VAPPRISR` | `ACCVAPPRIS` | IMBUDGET>=2 |
| `VAPSNOWT` | `VAPSNOWR` | `ACCVAPSNOW` | IMBUDGET>=2 |
| `VAPAGGRT` | `VAPAGGRR` | `ACCVAPAGGR` | IMBUDGET>=2 |
| `VAPGRAUT` | `VAPGRAUR` | `ACCVAPGRAU` | IMBUDGET>=2 |
| `VAPHAILT` | `VAPHAILR` | `ACCVAPHAIL` | IMBUDGET>=2 |
| `VAPDRIZT` | `VAPDRIZR` | `ACCVAPDRIZ` | IMBUDGET>=2 |
| `EVAPCLDT` | `EVAPCLDR` | `ACCEVAPCLD` | IMBUDGET>=2 |
| `EVAPRAINT` | `EVAPRAINR` | `ACCEVAPRAIN` | IMBUDGET>=2 |
| `EVAPPRIST` | `EVAPPRISR` | `ACCEVAPPRIS` | IMBUDGET>=2 |
| `EVAPSNOWT` | `EVAPSNOWR` | `ACCEVAPSNOW` | IMBUDGET>=2 |
| `EVAPAGGRT` | `EVAPAGGRR` | `ACCEVAPAGGR` | IMBUDGET>=2 |
| `EVAPGRAUT` | `EVAPGRAUR` | `ACCEVAPGRAU` | IMBUDGET>=2 |
| `EVAPHAILT` | `EVAPHAILR` | `ACCEVAPHAIL` | IMBUDGET>=2 |
| `EVAPDRIZT` | `EVAPDRIZR` | `ACCEVAPDRIZ` | IMBUDGET>=2 |
| `MELTPRIST` | `MELTPRISR` | `ACCMELTPRIS` | IMBUDGET>=2 |
| `MELTSNOWT` | `MELTSNOWR` | `ACCMELTSNOW` | IMBUDGET>=2 |
| `MELTAGGRT` | `MELTAGGRR` | `ACCMELTAGGR` | IMBUDGET>=2 |
| `MELTGRAUT` | `MELTGRAUR` | `ACCMELTGRAU` | IMBUDGET>=2 |
| `MELTHAILT` | `MELTHAILR` | `ACCMELTHAIL` | IMBUDGET>=2 |
| `RIMECLDSNOWT` | `RIMECLDSNOWR` | `ACCRIMECLDSNOW` | IMBUDGET>=2 |
| `RIMECLDAGGRT` | `RIMECLDAGGRR` | `ACCRIMECLDAGGR` | IMBUDGET>=2 |
| `RIMECLDGRAUT` | `RIMECLDGRAUR` | `ACCRIMECLDGRAU` | IMBUDGET>=2 |
| `RIMECLDHAILT` | `RIMECLDHAILR` | `ACCRIMECLDHAIL` | IMBUDGET>=2 |
| `RAIN2PRT` | `RAIN2PRR` | `ACCRAIN2PR` | IMBUDGET>=2 |
| `RAIN2SNT` | `RAIN2SNR` | `ACCRAIN2SN` | IMBUDGET>=2 |
| `RAIN2AGT` | `RAIN2AGR` | `ACCRAIN2AG` | IMBUDGET>=2 |
| `RAIN2GRT` | `RAIN2GRR` | `ACCRAIN2GR` | IMBUDGET>=2 |
| `RAIN2HAT` | `RAIN2HAR` | `ACCRAIN2HA` | IMBUDGET>=2 |
| `AGGRSELFPRIST` | `AGGRSELFPRISR` | `ACCAGGRSELFPRIS` | IMBUDGET>=2 |
| `AGGRSELFSNOWT` | `AGGRSELFSNOWR` | `ACCAGGRSELFSNOW` | IMBUDGET>=2 |
| `AGGRPRISSNOWT` | `AGGRPRISSNOWR` | `ACCAGGRPRISSNOW` | IMBUDGET>=2 |
| `DUST1CLDRT` | `DUST1CLDRR` | `ACCDUST1CLDR` | IMBUDGET=3, IDUST>=1 |
| `DUST2CLDRT` | `DUST2CLDRR` | `ACCDUST2CLDR` | IMBUDGET=3, IDUST>=1 |
| `DUST1DRZRT` | `DUST1DRZRR` | `ACCDUST1DRZR` | IMBUDGET=3, IDUST>=1 |
| `DUST2DRZRT` | `DUST2DRZRR` | `ACCDUST2DRZR` | IMBUDGET=3, IDUST>=1 |
| `REGEN1CLDRT` | `REGEN1CLDRR` | `ACCREGEN1CLDR` | IMBUDGET=3, ICCNLEV>=2 |
| `REGEN2CLDRT` | `REGEN2CLDRR` | `ACCREGEN2CLDR` | IMBUDGET=3, ICCNLEV>=2 |
| `REGEN1DRZRT` | `REGEN1DRZRR` | `ACCREGEN1DRZR` | IMBUDGET=3, ICCNLEV>=2 |
| `REGEN2DRZRT` | `REGEN2DRZRR` | `ACCREGEN2DRZR` | IMBUDGET=3, ICCNLEV>=2 |
| `NUCCLDCT` | `NUCCLDCR` | `ACCNUCCLDC` | LEVEL=4, IMBUDGET>=1 or 2 |
| `NUCICECT` | `NUCICECR` | `ACCNUCICEC` | LEVEL=4, IMBUDGET>=1 or 2 |
| `INUCHOMCT` | `INUCHOMCR` | `ACCINUCHOMC` | LEVEL=4, IMBUDGET>=1 or 2 |
| `INUCIFNCT` | `INUCIFNCR` | `ACCINUCIFNC` | LEVEL=4, IMBUDGET>=1 or 2 |
| `REGEN_AERO1_DN` | `REGEN_AERO1_NR` | `ACCREGEN_AERO1_N` | IMBUDGET>=1, ICCNLEV>=2 |
| `REGEN_AERO1_DM` | `REGEN_AERO1_MR` | `ACCREGEN_AERO1_M` | IMBUDGET>=1, ICCNLEV>=2 |
| `REGEN_AERO2_DN` | `REGEN_AERO2_NR` | `ACCREGEN_AERO2_N` | IMBUDGET>=1, ICCNLEV>=2 |
| `REGEN_AERO2_DM` | `REGEN_AERO2_MR` | `ACCREGEN_AERO2_M` | IMBUDGET>=1, ICCNLEV>=2 |

The bin microphysics allocates only a subset of these processes. See the
`level == 4` block of `alloc_micro()` in `mem_micro.f90` for the exact set.

## Adding a new budget process

1. Add its base name to the declaration lists in `mem_micro.f90`: `<name>r` and
   `acc<name>`.
2. Allocate `<name>r` under the right `IMBUDGET` gate. Allocate `acc<name>` under
   the same gate combined with `imbudgetacc==1`.
3. Register both arrays with `vtables2`, and add both to `dealloc_micro()`.
4. Declare the work vector `x<name>` in `micphys.f90`, and zero it in
   `range_check()`.
5. Call `store_budget_rate()` and `add_to_budget_total()` for it in
   `copyback()`.
6. Add both arrays to `zero_micro_budget_point()` (`rbnd_nonscalar.f90`).
