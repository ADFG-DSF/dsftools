# dsftools 0.1.2 (June 2024)

Added a couple of miscellaneous functions:

* `rp()` for calculating confidence or (relative) accuracy from a vector of simulated values and a true value

# dsftools 0.1.1 (Apr 2024)

Split out the simulation routine in `verify_ASL_table()` and made the original
function a wrapper.  Creation of two additional function:

* `simulate_ASL_table()` specifically for the simulation portion
* `rp_ASL_table()` for describing the relative precision of estimates via simulation

# dsftools 0.1.0 (Feb 2024)

Initial creation.  Functions currently include

* For ASL summaries: `ASL_table()`, `verify_ASL_table()`, and `ASL_boilerplate()`
* Miscellaneous functions: `logit()`, `expit()`, and `se()`
