The FLEXPARAMS ramsin argument allows for changing arbitrary values in the source code without needing to recompile RAMS, like the flex params in CM1. It is a 20-element real 1D array, so can only accept floating-point numbers. (It can easily be made longer than this if need be; this is the main reason I chose to implement this as an array rather than as 20 individual parameters like CM1 does.) It needs to be specified in the RAMSIN, but does not need to have all 20 values specified, it just needs at least one. If you're not using any flexparams then this value is meaningless. The values can be accessed within any subroutine by importing them via `use mem_flexparams, only: flexparams`

The purpose of this document is also to describe the particular usage of the FLEXPARAMS in the version of RAMS to which this document is attached. The use of each parameter 1-20, if used, should be described below.

1: The value of ccn1_release_start_z in mic_init.f90, which is the start of the height at which the CCN1 profile begins linearly decreasing to 0. Params 1-3 are all used in this file.
2: ccn1_release_end_z, the height at which the CCN1 profile decrease to its release value
3: ccn1_release_value, the concentration in #/mg of CCN1 above ccn1_release_end_z
4: The value of n_z_points_per_tracer in mic_init.f90, which is the number of z levels encompassed in each tracer species (used for determining vertical origins of air)
5: Height over which random initial theta perturbations (i.e. BUBBLE=3 [or BUBBLE=4 with my modifications]) linearly decrease to 0. Default RAMS value for this was 500 m. Specifically, sets random_perturbation_max_z in the `bubble` subroutine in ruser.f90.
6: n_bubble_shells in mic_init.f90: the number of equal-volume tracer shells filling the warm bubble (species 1 is the innermost). 0 (the default) keeps the original slab mode driven by param 4. Params 6-9 are all used by subroutine init_tracer, and require the cosine-squared bubble (IBUBBLE=2 or 4) because the shells are nested on the same normalized ellipsoidal radius the bubble is built on, read from the same IBD* namelist entries.
7: n_outer_shells, the number of equal-volume tracer shells outside the bubble, numbered outward after the bubble's own shells. These measure the near-field air a thermal entrains. 0 for none.
8: outer_shell_limit, the outermost normalized bubble radius those outer shells reach (e.g. 2.0 = twice the bubble radius). Must exceed 1 when param 7 is nonzero.
9: Set to 1 to make the species following the shells an environment tracer, filling everything beyond outer_shell_limit. That is what makes a boundary's purity measurable: how much of the air inside a fitted thermal boundary never belonged to the bubble.
10:
11:
12:
13:
14:
15:
16:
17:
18:
19:
20: 
