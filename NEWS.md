# dsftools 0.1.4 (July 2024)

Added two functions:

* `binomial_detection()` and `multinomial_detection()`: in the context of a telemetry study, estimating the probability of detecting a single area or multiple areas used by some proportion of the marked population, given the sample size of instrumented fish, assuming random sampling

# dsftools 0.1.3 (July 2024)

Added some shortcut functions for vector subsetting:

* `%inside%` behaves similar to %in%, but with a numeric vector; `%inside()%`, `%inside[)%`, and `%inside(]%` allow specification of interval endpoints
* `%s_l%, `%s_leq%, `%s_g%, `%s_geq%, `%s_inside%, `%s_inside()%, `%s_inside[)%, and `%s_inside[)% to shortcut vector subsetting

# dsftools 0.1.2 (June 2024)

Added a couple of miscellaneous functions:

* `rp()` for calculating confidence or (relative) accuracy from a vector of simulated values and a true value
* `plotcor()` for plotting a correlation matrix

# dsftools 0.1.1 (Apr 2024)

Split out the simulation routine in `verify_ASL_table()` and made the original
function a wrapper.  Creation of two additional functions:

* `simulate_ASL_table()` specifically for the simulation portion
* `rp_ASL_table()` for describing the relative precision of estimates via simulation

# dsftools 0.1.0 (Feb 2024)

Initial creation.  Functions currently include

* For ASL summaries: `ASL_table()`, `verify_ASL_table()`, and `ASL_boilerplate()`
* Miscellaneous functions: `logit()`, `expit()`, and `se()`
