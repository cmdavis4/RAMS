# Momentum budget output (IUVWTEND)

This note describes the wind-tendency output that `IUVWTEND` controls on this
branch. All values are rates in m/s².

## Flag values

| `IUVWTEND` | Output |
|---|---|
| 0 | No momentum budget output. |
| 1 | Only the instantaneous accelerations `UP_DUDT`, `VP_DVDT` and `WP_DWDT`. |
| 2 | The accelerations, plus every individual budget term (listed below). |

`IMBUDGET` no longer affects any momentum output. Before this change,
`IMBUDGET>=1` also turned on `WP_BUOY_THETA`, `WP_BUOY_COND` and `WP_ADVDIF`.

## The accelerations (IUVWTEND >= 1)

RAMS steps the winds with leapfrog. At each step the model evaluates the
tendency at time level n and applies it across two step lengths, from level
n−1 to level n+1. The accelerations are the centred difference over that
interval:

```
UP_DUDT = (u(n+1) - u(n-1)) / (2*dtlt)
```

`VP_DVDT` and `WP_DWDT` follow the same formula. u(n+1) is the `UC` that is
written with the same output time.

This is the model's own instantaneous time derivative, and it equals the sum
of the `IUVWTEND=2` terms. The result differs from the budget sum only at
points where boundary conditions reset the wind, such as the lowest and
highest levels and the lateral boundaries.

The accelerations are not `(UC - UP) / dtlt`. That quantity is the change
since the previous step, and you can compute it from the written `UC` and
`UP` fields. u(n−1) is overwritten during the step and never written, so
only the model can form the centred difference.

`save_uvw_before_acoustic()` in `raco.f90` stores u(n−1). It runs just before
the acoustic solver in `timestep()` (`rtimh.f90`). `compute_uvw_accel()` forms
the difference at the end of the step.

## The budget terms (IUVWTEND = 2)

| Component | Terms |
|---|---|
| u | `UP_PGFORCE`, `UP_ADVECTION`, `UP_DIFFUSION`, `UP_CORIOLIS`, `UP_RAYLEIGH` |
| v | `VP_PGFORCE`, `VP_ADVECTION`, `VP_DIFFUSION`, `VP_CORIOLIS`, `VP_RAYLEIGH` |
| w | `WP_PGFORCE`, `WP_ADVECTION`, `WP_DIFFUSION`, `WP_RAYLEIGH`, `WP_BUOY_THETA`, `WP_BUOY_COND` |

`WP_ADVDIF` is the whole w tendency before buoyancy is added. It is kept for
comparison with older runs and is not a separate term in the sum.

`WP_PGFORCE` is computed as a residual. The vertical pressure-gradient force
is applied partly implicitly inside the acoustic solver, so no single array
ever holds it. It is therefore the acoustic-solver increment in w minus every
other w term.
