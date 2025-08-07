#' Construct an ASL Table from Vectors of Data
#' @description This function is intended to perform all the necessary calculations
#' to construct tables summarizing ASL (Age, Sex, and Length) for several data
#' scenarios:
#'
#' * Stratified vs. non-stratified
#' * Abundance known without error vs. estimated with error vs. unknown
#' * Data available for Age, Sex, and Length vs. some subset.
#'
#' The function will return a table with rows for each combination of Age, Sex, and
#' Length (as available), and will summarize the respective proportions, estimated
#' abundance, and associated Lengths (as available).
#'
#' Generally, `NULL` values in a given data vector indicate that the vector will
#' not be used.
#'
#' Methods are also described in the document found here:
#'
#' `https://github.com/ADFG-DSF/dsftools/blob/main/addl_documentation/asl_equations.Rmd`
#'
#' @param age Vector of ages.  If the default (`NULL`) is accepted, age will be
#' omitted from the resulting table.
#' @param sex Vector of sex.  If the default (`NULL`) is accepted, sex will be
#' omitted from the resulting table.
#' @param length Vector of length, treated as NUMERIC.  If the default (`NULL`) is accepted, length
#' will be omitted from the resulting table.
#' @param length_cat Vector of length, treated as CATEGORICAL.  If the default (`NULL`) is accepted, length
#' will be omitted from the resulting table.
#' @param stratum Optional vector of stratum, if stratified estimators will be
#' used.  This must be formatted as positive whole numbers (1, 2, 3, etc).
#' Defaults to `NULL`.
#' @param Nhat Optional vector of abundance for each stratum, if stratified
#' estimators will be used.  The length of this vector must correspond to the
#' maximum value in the `stratum=` argument.  Defaults to `NULL`.
#' @param se_Nhat Optional vector of the standard error of abundance for each
#' stratum, if stratified estimators will be used and by-stratum abundance is
#' considered to have measurement error.    The length of this vector must
#' correspond to the maximum value in the `stratum=` argument.  Defaults to `NULL`.
#' @param stratum_weights Optional vector of weights for each stratum, if
#' relative weights are known but abundance is not.  Using this argument rather
#' than `Nhat` will also omit the Finite Population Correction factor (FPC) from
#' calculations, but weights are treated as constant (without error), likely
#' underestimating variance.  Defaults to `NULL`.
#' @param verbose Whether to print messages corresponding to the method used.
#' Defaults to `FALSE`.
#' @param FPC Whether to incorporate the Finite Population Correction factor (FPC) in
#' variance calculations.  Allowed values are
#' * `"ifknown"` (the default), which will use the FPC only if abundance is known without error
#' * `"always"`, which will use the FPC wherever possible (i.e. if there is an estimate of abundance)
#' * `"never"`, which will always ignore the FPC.
#' used if the abundance is considered to be known without error.
#' @return A data.frame with rows corresponding to categories of age, sex, and/or length,
#' depending on data inputs.
#'
#' ## Stratified - If abundance is known without error
#'
#' ### If there are proportions
#'
#' The proportion of each age, sex, and/or length category *z* will be estimated for each sampling stratum *t* as follows:
#'
#'  \deqn{\hat{p}_{tz}=\frac{n_{tz}}{n_t}}
#'
#' in which \eqn{n_{tz}} equals the number of fish sampled during sampling stratum \eqn{t} classified as age, sex, and/or length category \eqn{z}, and \eqn{n_t} equals the number of fish sampled for age, sex, and/or length determination within sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\hat{p}_{tz}} will be estimated as the following (Cochran 1977), in which \eqn{N_t} represents the total abundance of fish in sampling stratum *t*:
#'
#' \deqn{\hat{var}[\hat{p}_{tz}]=\frac{\hat{p}_{tz}(1-\hat{p}_{tz})}{n_t-1}\left(\frac{N_t-n_t}{N_t-1}\right)}
#'
#' if the finite population correction factor (FPC) is used; otherwise, as the following:
#'
#' \deqn{\hat{var}[\hat{p}_{tz}]=\frac{\hat{p}_{tz}(1-\hat{p}_{tz})}{n_t-1}}
#'
#' The total abundance by age, sex, and/or length category in each sampling stratum will be estimated as follows:
#'
#' \deqn{\hat{N}_{tz}=N_t\hat{p}_{tz}}
#'
#' with variance estimated as
#'
#' \deqn{\hat{var}[\hat{N}_{tz}]=N_t^2\hat{var}[\hat{p}_{tz}]}
#'
#' The total abundance by age, sex, and/or length category and its variance will then be estimated by summation as follows:
#'
#' \deqn{\hat{N}_z=\sum_{t=1}^{L}\hat{N}_{tz}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{N}_{z}]=\sum_{t=1}^{L}\hat{var}[\hat{N}_{tz}]}
#'
#' where \eqn{L} equals the number of sampling strata.
#'
#' Finally, the overall proportion by age, sex, and/or length category and its variance will be estimated as follows:
#'
#' \deqn{\hat{p}_z=\frac{\hat{N}_z}{N}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{p}_z]=\frac{\hat{var}[\hat{N}_z]}{N^2}}
#'
#' where \eqn{N} is the total abundance across all sampling periods.
#'
#' The mean length by age, sex, and/or length for each sampling stratum will be estimated as follows:
#'
#' \deqn{\bar{x}_{tz}=\frac{\sum_{i=1}^{n_{tz}}x_{tzi}}{n_{tz}}}
#'
#' where \eqn{x_{tzi}} is the length of the *i*th fish sampled of age, sex, and/or length \eqn{z} during sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\bar{x}_{tz}} will be estimated as
#'
#' \deqn{\hat{var}[\bar{x}_{tz}]=\frac{\sum_{i=1}^{n_{tz}}(x_{tzi}-\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}\left(\frac{
#'     \hat{N}_{tz}-n_{tz}}{\hat{N}_{tz}-1}\right)}
#'
#' if the finite population correction factor (FPC) will be used; otherwise, as the following:
#'
#' \deqn{\hat{var}[\bar{x}_{tz}]=\frac{\sum_{i=1}^{n_{tz}}(x_{tzi}-\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}}
#'
#' The mean length by age, sex, and/or length category will then be estimated as follows:
#'
#' \deqn{\bar{x}_z=\sum_{t=1}^{L}\frac{\hat{N}_{tz}}{\hat{N}_z}\bar{x}_{tz}}
#'
#' with its variance approximated using a Taylor's series expansion (Mood et al. 1974):
#'
#' \deqn{\hat{var}[\bar{x}_z]\approx\sum_{t=1}^{L}\frac{\hat{N}_{tz}^2}{\hat{N}_z^2}\hat{var}[\bar{x}_{tz}]+\sum_{t=1}^{L}\frac{\left(\bar{x}_{tz}\hat{N}_z-\left(\sum_{u=1}^{L}\bar{x}_{uz}\hat{N}_{uz}\right)\right)^2}{\hat{N}_z^4}\hat{var}[\hat{N}_{tz}]}
#'
#' ### If there are no proportions to estimate
#'
#' The mean length for each sampling stratum will be estimated as follows, where \eqn{x_{ti}} is the length of the *i*th fish sampled within sampling stratum \eqn{t}, and \eqn{n_t} is the number of fish in stratum *t* sampled for length:
#'
#' \deqn{\bar{x}_{t}=\frac{\sum_{i=1}^{n_{t}}x_{ti}}{n_{t}}}
#'
#' The sampling variance of \eqn{\bar{x}_{t}} will be estimated as
#'
#' \deqn{\hat{var}[\bar{x}_{t}]=\frac{\sum_{i=1}^{n_{t}}(x_{ti}-\bar{x}_{t})^2}{n_{t}(n_{t}-1)}\left(\frac{
#' N_{t}-n_{t}}{N_{t}-1}\right)}
#'
#' if abundance per stratum \eqn{N_t} is known and if the finite population correction factor is used, otherwise as:
#'
#' \deqn{\hat{var}[\bar{x}_{t}]=\frac{\sum_{i=1}^{n_{t}}(x_{ti}-\bar{x}_{t})^2}{n_{t}(n_{t}-1)}}
#'
#' Stratified estimates of mean length will be calculated as follows, in which \eqn{N_t} represents the abundance associated with sampling stratum \eqn{t}, \eqn{N} represents the total abundance, and *L* represents the number of sampling strata:
#'
#' \deqn{\bar{x}=\frac{1}{N}\sum_{t=1}^L N_t\bar{x}_t}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}]=\sum_{t=1}^L\left(\frac{N_t}{N}\right)^2 \hat{var}[\bar{x}_t]}
#'
#' ## Stratified - If abundance is estimated with error
#'
#' ### If there are proportions
#'
#' The proportion of each age, sex, and/or length category *z* will be estimated for each sampling stratum *t* as follows:
#'
#' \deqn{\hat{p}_{tz}=\frac{n_{tz}}{n_t}}
#'
#' in which \eqn{n_{tz}} equals the number of fish sampled during sampling stratum \eqn{t} classified as age, sex, and/or length category \eqn{z}, and \eqn{n_t} equals the number of fish sampled for age, sex, and/or length determination within sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\hat{p}_{tz}} will be estimated as the following (Cochran 1977) in which \eqn{\hat{N}_t} is the estimated abundance of fish in sampling stratum \eqn{t}:
#'
#' \deqn{\hat{var}[\hat{p}_{tz}]=\frac{\hat{p}_{tz}(1-\hat{p}_{tz})}{n_t-1}\left(\frac{\hat{N}_t-n_t}{\hat{N}_t-1}\right)}
#'
#' if the finite population correction factor (FPC) is used; otherwise, as the following:
#'
#' \deqn{\hat{var}[\hat{p}_{tz}]=\frac{\hat{p}_{tz}(1-\hat{p}_{tz})}{n_t-1}}
#'
#' The total abundance by age, sex, and/or length category in each sampling stratum will be estimated as follows:
#'
#' \deqn{\hat{N}_{tz}=\hat{N}_t\hat{p}_{tz}}
#'
#' with variance estimated as (Goodman 1960):
#'
#' \deqn{\hat{var}[\hat{N}_{tz}]=\hat{N}_t^2\hat{var}[\hat{p}_{tz}] + \hat{p}_{tz}^2\hat{var}[\hat{N}_t]-\hat{var}[\hat{p}_{tz}]\hat{var}[\hat{p}_{tz}]}
#'
#' The total abundance by age, sex, and/or length category \eqn{z} and its variance will then be estimated by summation as follows:
#'
#' \deqn{\hat{N}_z=\sum_{t=1}^{L}\hat{N}_{tz}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{N}_{z}]=\sum_{t=1}^{L}\hat{var}[\hat{N}_{tz}]}
#'
#' where \eqn{L} equals the number of sampling strata.
#'
#' Finally, the overall proportion by age, sex, and/or length category and its variance will be estimated as follows:
#'
#' \deqn{\hat{p}_z=\frac{\hat{N}_z}{\sum_{t=1}^{L}\hat{N}_t}}
#'
#' with variance estimated by the delta method (Casella & Berger 2002) as:
#'
#' \deqn{\hat{var}[\hat{p}_z] \approx \left(\frac{\hat{N}_z}{\sum_{t=1}^{L}\hat{N}_t}\right)^2\left(\frac{\hat{var}[\hat{N}_z]}{\hat{N}_z^2} + \frac{\sum_{t=1}^{L}\hat{var}[\hat{N}_t]}{(\sum_{t=1}^{L}\hat{N}_t)^2} - 2\frac{\hat{cov}[\hat{N}_z,\sum_{t=1}^{L}\hat{N}_t]}{\hat{N}_z\sum_{t=1}^{L}\hat{N}_t}\right)}
#'
#' in which
#'
#' \deqn{\hat{cov}[\hat{N}_z,\sum_{t=1}^{L}\hat{N}_t]=\sum_{t=1}^{L}\hat{p}_{tz}\hat{var}[\hat{N_t}]}
#'
#' The mean length by age, sex, and/or length for each sampling stratum will be estimated as follows:
#'
#' \deqn{\bar{x}_{tz}=\frac{\sum_{i=1}^{n_{tz}}x_{tzi}}{n_{tz}}}
#'
#' where \eqn{x_{tzi}} is the length of the *i*th fish sampled of age, sex, and/or length \eqn{z} during sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\bar{x}_{tz}} will be estimated as
#'
#' \deqn{\hat{var}[\bar{x}_{tz}]=\frac{\sum_{i=1}^{n_{tz}}(x_{tzi}-\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}\left(\frac{\hat{N}_{tz}-n_{tz}}{\hat{N}_{tz}-1}\right)}
#'
#' if the finite population correction factor (FPC) is used; otherwise, as
#'
#' \deqn{\hat{var}[\bar{x}_{tz}]=\frac{\sum_{i=1}^{n_{tz}}(x_{tzi}-\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}}
#'
#' The mean length by age, sex, and/or length category will then be estimated as follows:
#'
#' \deqn{\bar{x}_z=\sum_{t=1}^{L}\frac{\hat{N}_{tz}}{\hat{N}_z}\bar{x}_{tz}}
#'
#' with its variance approximated using a Taylor's series expansion (Mood et al. 1974):
#'
#' \deqn{\hat{var}[\bar{x}_z]\approx\sum_{t=1}^{L}\frac{\hat{N}_{tz}^2}{\hat{N}_z^2}\hat{var}[\bar{x}_{tz}]+\sum_{t=1}^{L}\frac{\left(\bar{x}_{tz}\hat{N}_z-\left(\sum_{u=1}^{L}\bar{x}_{uz}\hat{N}_{uz}\right)\right)^2}{\hat{N}_z^4}\hat{var}[\hat{N}_{tz}]}
#'
#' ### If there are no proportions to estimate
#'
#' The mean length for each sampling stratum will be estimated as follows, where \eqn{x_{ti}} is the length of the *i*th fish sampled within sampling stratum \eqn{t}, and \eqn{n_t} is the number of fish in stratum *t* sampled for length:
#'
#' \deqn{\bar{x}_{t}=\frac{\sum_{i=1}^{n_{t}}x_{ti}}{n_{t}}}
#'
#' The sampling variance of \eqn{\bar{x}_{t}} will be estimated as
#'
#' \deqn{\hat{var}[\bar{x}_{t}]=\frac{\sum_{i=1}^{n_{t}}(x_{ti}-\bar{x}_{t})^2}{n_{t}(n_{t}-1)}\left(\frac{
#' \hat{N}_{t}-n_{t}}{\hat{N}_{t}-1}\right)}
#'
#' if the finite population correction factor (FPC) is used, otherwise as:
#'
#' \deqn{\hat{var}[\bar{x}_{t}]=\frac{\sum_{i=1}^{n_{t}}(x_{ti}-\bar{x}_{t})^2}{n_{t}(n_{t}-1)}}
#'
#' Stratified estimates of mean length will be calculated as follows, in which \eqn{\hat{N}_t} and \eqn{\bar{x}_t} represent the estimated abundance and mean length associated with stratum *t*, respectively:
#'
#' \deqn{\bar{x}=\frac{\sum_{t=1}^L \hat{N}_t\bar{x}_t}{\sum_{t=1}^L \hat{N}_t}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}] \approx \left(\frac{\sum_{t=1}^L \hat{N}_t\bar{x}_t}{\sum_{t=1}^L \hat{N}_t}\right)^2\left(\frac{\hat{var}[\sum_{t=1}^L \hat{N}_t\bar{x}_t]}{\left(\sum_{t=1}^L \hat{N}_t\bar{x}_t\right)^2}+\frac{\sum_{t=1}^L \hat{var}[\hat{N}_t]}{\left(\sum_{t=1}^L \hat{N}_t\right)^2}-2\frac{\hat{cov}[\sum_{t=1}^L \hat{N}_t,\sum_{t=1}^L N_t\bar{x}_t]}{\left(\sum_{t=1}^L \hat{N}_t\right)\left(\sum_{t=1}^L \hat{N}_t\bar{x}_t\right)}\right)}
#'
#' in which
#'
#' \deqn{\hat{cov}[\sum_{t=1}^L \hat{N}_t,\sum_{t=1}^L \hat{N}_t\bar{x}_t]=\sum_{t=1}^L \bar{x}_t\hat{var}[\hat{N}_t]}
#'
#' and
#'
#' \deqn{\hat{var}[\sum_{t=1}^L \hat{N}_t\bar{x}_t]=\sum_{t=1}^L\hat{N}_t^2\hat{var}[\bar{x}_t] + \bar{x}_t^2\hat{var}[\hat{N}_t]-\hat{var}[\hat{N}_t]\hat{var}[\bar{x}_t]}
#'
#' by means of the delta method (Casella & Berger 2002) and Goodman (1960), respectively.
#'
#' ## Stratified - If abundance is unknown and sample weights are used
#'
#' ### If there are proportions
#'
#' The proportion of each age, sex, and/or length category *z* will be estimated for each sampling stratum *t* as follows:
#'
#' \deqn{\hat{p}_{tz}=\frac{n_{tz}}{n_t}}
#'
#' in which \eqn{n_{tz}} equals the number of fish sampled during sampling stratum \eqn{t} classified as age, sex, and/or length category \eqn{z}, and \eqn{n_t} equals the number of fish sampled for age, sex, and/or length determination within sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\hat{p}_{tz}} will be estimated as the following (Cochran 1977):
#'
#' \deqn{\hat{var}[\hat{p}_{tz}]=\frac{\hat{p}_{tz}(1-\hat{p}_{tz})}{n_t-1}}
#'
#' The overall proportion by age, sex, and/or length category and its variance will be estimated as follows, in which \eqn{w_t} represents the sampling weight associated with stratum *t* and *L* equals the number of strata.  It is worth noting that weights \eqn{w_t} are treated as constant (i.e. known without error), therefore all variance estimates must be interpreted as minima without further assumptions.
#'
#' \deqn{\hat{p}_z=\frac{\sum_{t=1}^Lw_t\hat{p}_{tz}}{\sum_{t=1}^Lw_t}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{p}_z]=\frac{\sum_{t=1}^Lw_t^2\hat{var}[\hat{p}_{tz}]}{\left(\sum_{t=1}^Lw_t\right)^2}}
#'
#' The mean length by age, sex, and/or length for each sampling stratum will be estimated as follows:
#'
#' \deqn{\bar{x}_{tz}=\frac{\sum_{i=1}^{n_{tz}}x_{tzi}}{n_{tz}}}
#'
#' where \eqn{x_{tzi}} is the length of the *i*th fish sampled of age, sex, and/or length \eqn{z} during sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\bar{x}_{tz}} will be estimated as:
#'
#' \deqn{\hat{var}[\bar{x}_{tz}]=\frac{\sum_{i=1}^{n_{tz}}(x_{tzi}-\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}}
#'
#' The mean length by age, sex, and/or length category will then be estimated as follows:
#'
#' \deqn{\bar{x}_z=\frac{\sum_{t=1}^{L}w_t\hat{p}_{tz}\bar{x}_{tz}}{\sum_{t=1}^{L}w_t\hat{p}_{tz}}}
#'
#' with its variance approximated using a Taylor's series expansion (Mood et al. 1974):
#'
#' \deqn{\hat{var}[\bar{x}_z]\approx\frac{\sum_{t=1}^{L}w_t\hat{p}_{tz}^2\hat{var}[\bar{x}_{tz}]}{\left(\sum_{t=1}^Lw_t\hat{p}_{tz}\right)^2}+\frac{\sum_{t=1}^{L}\left(\bar{x}_{tz}\sum_{u=1}^Lw_u\hat{p}_{uz}-\left(\sum_{u=1}^{L}\bar{x}_{uz}w_u\hat{p}_{uz}\right)\right)^2w_t^2\hat{var}[\hat{p}_{tz}]}{\left(\sum_{t=1}^Lw_t\hat{p}_{tz}\right)^4}}
#'
#' ### If there are no proportions to estimate
#'
#' The mean length for each sampling stratum will be estimated as follows, where \eqn{x_{ti}} is the length of the *i*th fish sampled within sampling stratum \eqn{t}, and \eqn{n_t} is the number of fish in stratum *t* sampled for length:
#'
#' \deqn{\bar{x}_{t}=\frac{\sum_{i=1}^{n_{t}}x_{ti}}{n_{t}}}
#'
#' The sampling variance of \eqn{\bar{x}_{t}} will be estimated as:
#'
#' \deqn{\hat{var}[\bar{x}_{t}]=\frac{\sum_{i=1}^{n_{t}}(x_{ti}-\bar{x}_{t})^2}{n_{t}(n_{t}-1)}}
#'
#' Stratified estimates of mean length will be calculated as follows, in which \eqn{w_t} and \eqn{\bar{x}_t} represent the sampling weight and average length associated with sampling stratum \eqn{t}, respectively.  It is worth noting that weights \eqn{w_t} are treated as constant (i.e. known without error), therefore all variance estimates must be interpreted as minima without further assumptions.
#'
#' \deqn{\bar{x}=\frac{\sum_{t=1}^L w_t\bar{x}_t}{\sum_{t=1}^L w_t}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}]=\frac{\sum_{t=1}^L w_t^2\hat{var}[\bar{x}_t]}{\left(\sum_{t=1}^L w_t\right)^2}}
#'
#' ## Pooled (not stratified) - If abundance is known without error
#'
#' ### If there are proportions
#'
#' #' Proportions of each age, sex, and/or length category \eqn{z} will be estimated as follows (Cochran 1977):
#'
#' \deqn{\hat{p}_z=\frac{n_z}{n}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{p}_z]=\frac{\hat{p}(1-\hat{p})}{n-1}\left(\frac{N-n}{N-1}\right)}
#'
#' if total abundance \eqn{N} is known and the finite population correction factor (FPC) is used; otherwise as
#'
#' \deqn{\hat{var}[\hat{p}_z]=\frac{\hat{p}_z(1-\hat{p}_z)}{n-1}}
#'
#' in which \eqn{n_z} denotes the number of fish sampled in age, sex, and/or length category \eqn{z}, \eqn{n} denotes the total number of fish sampled, and *N* denotes the total abundance.
#'
#' Total abundance for age, sex, and/or length category \eqn{z} will be estimated as
#'
#' \deqn{\hat{N}_z=N\hat{p}_z} and
#'
#' \deqn{\hat{var}[\hat{N}_z]=N^2\hat{var}[\hat{p}_z]}
#'
#' The mean length associated with age, sex, and/or length category \eqn{z} will be estimated as the following, in which \eqn{x_{zi}} represents the length of the *i*th fish in age, sex, and/or length category *z* and \eqn{n_z} represents the number of fish in age, sex, and/or length category *z* with an associated length measurement:
#'
#' \deqn{\bar{x}_z=\frac{\sum_{i=1}^{n_z}x_{zi}}{n_z}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}_z]=\frac{\sum_{i=1}^{n_z}(x_{zi}-\bar{x_z})^2}{n_z(n_z-1)}\left(\frac{\hat{N}_z-n_z}{\hat{N}_z-1}\right)}
#'
#' if the finite population correction factor (FPC) is used, otherwise as:
#'
#' \deqn{\hat{var}[\bar{x}_z]=\frac{\sum_{i=1}^{n_z}(x_{zi}-\bar{x_z})^2}{n_z(n_z-1)}}
#'
#' ### If there are no proportions to estimate
#'
#' The mean length of all fish will be estimated as the following, in which \eqn{x_{i}} represents the length of the *i*th fish, *n* represents the number of fish with an associated length measurement, and *N* represents the total abundance:
#'
#' \deqn{\bar{x}=\frac{\sum_{i=1}^{n}x_{i}}{n}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}]=\frac{\sum_{i=1}^{n}(x_{i}-\bar{x})^2}{n(n-1)}\left(\frac{N-n}{N-1}\right)}
#'
#' if the finite population correction factor (FPC) is used, otherwise as:
#'
#' \deqn{\hat{var}[\bar{x}]=\frac{\sum_{i=1}^{n}(x_{i}-\bar{x})^2}{n(n-1)}}
#'
#' ## Pooled (not stratified) - If abundance is estimated with error
#'
#' ### If there are proportions
#'
#' Proportions of each age, sex, and/or length category \eqn{z} will be estimated as follows (Cochran 1977):
#'
#' \deqn{\hat{p}_z=\frac{n_z}{n}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{p}]=\frac{\hat{p}(1-\hat{p})}{n-1}\left(\frac{\hat{N}-n}{\hat{N}-1}\right)}
#'
#' if the finite population correction factor (FPC) is used, otherwise as
#'
#' \deqn{\hat{var}[\hat{p}]=\frac{\hat{p}(1-\hat{p})}{n-1}}
#'
#' in which \eqn{n_z} denotes the number of fish sampled in age, sex, and/or length category \eqn{z}, \eqn{n} denotes the total number of fish sampled, and \eqn{\hat{N}} denotes the estimated abundance.
#'
#' Total abundance for age, sex, and/or length category \eqn{z} will be estimated as follows (Goodman 1960):
#'
#' \deqn{\hat{N}_z=\hat{N}\hat{p}_z}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{N}_z]=\hat{N}^2\hat{var}[\hat{p}_z] + \hat{p}_z^2\hat{var}[\hat{N}]-\hat{var}[\hat{N}]\hat{var}[\hat{p}_z]}
#'
#' The mean length associated with age, sex, and/or length category \eqn{z} will be estimated as the following, in which \eqn{x_{zi}} represents the length of fish *i* within age, sex, and/or length category category *z*, and \eqn{n_z} denotes the number of fish within age, sex, and/or length category *z* with an associated length measurement:
#'
#' \deqn{\bar{x}_z=\frac{\sum_{i=1}^{n_z}x_{zi}}{n_z}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}_z]=\frac{\sum_{i=1}^{n_z}(x_{zi}-\bar{x_z})^2}{n_z(n_z-1)}\left(\frac{\hat{N}_z-n_z}{\hat{N}_z-1}\right)}
#'
#' if the finite population correction factor (FPC) is used; otherwise as
#'
#' \deqn{\hat{var}[\bar{x}_z]=\frac{\sum_{i=1}^{n_z}(x_{zi}-\bar{x_z})^2}{n_z(n_z-1)}}
#'
#' ### If there are no proportions to estimate
#'
#' The mean length of all fish will be estimated as the following, in which \eqn{x_{i}} represents the length of the *i*th fish, *n* represents the number of fish with an associated length measurement, and \eqn{\hat{N}} represents the estimated total abundance:
#'
#' \deqn{\bar{x}=\frac{\sum_{i=1}^{n}x_{i}}{n}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}]=\frac{\sum_{i=1}^{n}(x_{i}-\bar{x})^2}{n(n-1)}\left(\frac{\hat{N}-n}{\hat{N}-1}\right)}
#'
#' if the finite population correction factor (FPC) is used, otherwise as:
#'
#' \deqn{\hat{var}[\bar{x}]=\frac{\sum_{i=1}^{n}(x_{i}-\bar{x})^2}{n(n-1)}}
#'
#' ## Pooled (not stratified) - If abundance is unknown
#'
#' ### If there are proportions
#'
#' Proportions of each age, sex, and/or length category \eqn{z} will be estimated as follows (Cochran 1977):
#'
#' \deqn{\hat{p}_z=\frac{n_z}{n}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{p}_z]=\frac{\hat{p}_z(1-\hat{p}_z)}{n-1}}
#'
#' in which \eqn{n_z} denotes the number of fish sampled in age, sex, and/or length category \eqn{z}, and \eqn{n} denotes the total number of fish sampled.
#'
#' The mean length associated with age, sex, and/or length category \eqn{z} will be estimated as the following, in which \eqn{x_{zi}} denotes the length measurement associated with the *i*th fish in age, sex, and/or length category *z*:
#'
#' \deqn{\bar{x}_z=\frac{\sum_{i=1}^{n_z}x_{zi}}{n_z}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}_z]=\frac{\sum_{i=1}^{n_z}(x_{zi}-\bar{x_z})^2}{n_z(n_z-1)}}
#'
#' ### If there are no proportions to estimate
#'
#' The mean length will be estimated as the following, in which \eqn{x_{i}} denotes the length measurement associated with the *i*th fish and *n* denotes the number of fish with an associated length measurement:
#'
#' \deqn{\bar{x}=\frac{\sum_{i=1}^{n}x_{i}}{n}}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}]=\frac{\sum_{i=1}^{n}(x_{i}-\bar{x})^2}{n(n-1)}}
#'
#' @references Casella, George and Roger L. Berger. 2002. *Statistical Inference*. Australia ; Pacific Grove, CA, Thomson Learning
#'
#' Cochran, W. G.  1977.  *Sampling techniques*. 3rd edition.  John Wiley and Sons, New York.
#'
#' Goodman, L.A., 1960. On the exact variance of products. *Journal of the American statistical association, 55*(292), pp.708-713.
#'
#' Mood, A. M., F. A. Graybill, and D. C. Boes.  1974.  *Introduction to the theory of statistics*. 3rd edition.  McGraw-Hill Book Co., New York.
#' @author Matt Tyers
#' @seealso [verify_ASL_table]
#' @examples
#' ### a look at input
#' head(sim_data$data)
#' sim_data$abundance
#'
#' ### tables for possible scenarios
#'
#' # Stratified with error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           stratum=sim_data$data$stratum,
#'           Nhat=sim_data$abundance$Nhat,
#'           se_Nhat=sim_data$abundance$se_Nhat)
#'
#' # Stratified without error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           stratum=sim_data$data$stratum,
#'           Nhat=sim_data$abundance$Nhat)
#'
#' # Pooled (not stratified) with error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           Nhat=100*1000,
#'           se_Nhat=10*1000)
#'
#' # Pooled (not stratified) without error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           Nhat=100*1000)
#' @importFrom grDevices adjustcolor grey.colors
#' @importFrom graphics boxplot legend mosaicplot par points
#' @importFrom stats rnorm sd var
#' @export
ASL_table <- function(age=NULL,
                       sex=NULL,
                       length=NULL,
                       length_cat=NULL,
                       stratum=NULL,
                       Nhat=NULL,
                       se_Nhat=NULL,
                       stratum_weights = NULL,
                       verbose=FALSE,
                       FPC=c("ifknown", "always", "never")) {

  # -------- globally dealing with inputs ----------
  # combining age/sex categories as available
  notnulls <- sum(c(!is.null(sex), !is.null(age), !is.null(length_cat)))
  # if(!is.null(sex) & !is.null(age) & !is.null(length_cat)) {
  if(notnulls==3) {
    cats <- paste(sex, age, length_cat)
  } #else {
  if(notnulls==2) {
    if(!is.null(sex) & !is.null(age)) cats <- paste(sex, age)
    if(!is.null(sex) & !is.null(length_cat)) cats <- paste(sex, length_cat)
    if(!is.null(age) & !is.null(length_cat)) cats <- paste(age, length_cat)
  }
  if(notnulls==1) {
    if(!is.null(sex)) cats <- sex
    if(!is.null(age)) cats <- age
    if(!is.null(length_cat)) cats <- length_cat
  }
  if(notnulls==0) {
    cats <- rep("Total", length(length))
  }

  if(!is.null(stratum)) {
    # taking out data rows where strata value is NA
    length <- length[!is.na(stratum)]
    cats <- cats[!is.na(stratum)]
    stratum <- stratum[!is.na(stratum)]
  }

  ### need to constrain FPC to c("ifknown","always","never")
  FPC <- match.arg(FPC)

  # this will become output dataframe
  out <- data.frame(n=tapply(X=seq_along(cats), INDEX=cats, FUN=base::length))

  # --------- Stratified ----------
  if(!is.null(stratum)) {
    if(verbose) cat("\n", "Stratified","\n")
    # checking that inputs make sense

    # if(!inherits(stratum, "integer")) {
    #   stop("Stratum values must be positive whole numbers (1, 2, 3, etc.)")  ############### this fails
    # }


    if(length(stratum_weights) != max(stratum) & length(Nhat) != max(stratum)) {
      stop("Length of stratum weights does not equal max stratum value")
    }

    # -- Stratified without error in Nhats --
    if(is.null(se_Nhat)) {
      if(verbose) cat("\n", "without error in Nhat", "\n")
      if(is.null(Nhat) & !is.null(stratum_weights)) {
        Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
      } else {
        Nt <- Nhat # consistency with report eqns
      }

      # pulling these equations from the Jim Creek report
      if(length(unique(cats)) > 1) {

        # --- if there are proportions to estimate ---
        ntz <- table(stratum, cats, useNA="no")
        nt <- rowSums(ntz)
        if((is.null(Nhat) & !is.null(stratum_weights)) | FPC=="never") {
          FPC_vec <- 1  # this will ignore FPC in further calculations
          if(verbose) cat("\n","Finite Population Correction factor NOT used for proportions","\n")
        } else {
          FPC_vec <- (Nt-nt)/(Nt-1)
          if(verbose) cat("\n","Finite Population Correction factor USED for proportions","\n")
        }
        ptz <- ntz/nt
        vptz <- ptz*(1-ptz)/(nt-1)*FPC_vec
        Ntz <- Nt*ptz
        vNtz <- (Nt^2)*vptz
        Nz <- colSums(Ntz)
        vNz <- colSums(vNtz)
        pz <- Nz/sum(Nt)
        vpz <- vNz/(sum(Nt)^2)
        out$phat <- pz
        out$se_phat <- sqrt(vpz)
        if(!is.null(Nhat)) {
          out$Nhat <- Nz
          out$se_Nhat <- sqrt(vNz)
        }

        # summarizing length
        if(!is.null(length)) {
          xbartz <- tapply(length, list(stratum, cats), mean, na.rm=TRUE)
          vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          if((!is.null(Nhat) & is.null(stratum_weights)) & FPC=="always") {
            FPC_vec2 <- (Ntz-vxbartz_denom)/(Ntz-1)
            if(verbose) cat("\n","Finite Population Correction factor USED for means","\n")
          } else {
            FPC_vec2 <- 1#rep(1, length(nt))
            if(verbose) cat("\n","Finite Population Correction factor NOT used for means","\n")
          }
          vxbartz <- tapply(length, list(stratum, cats), var, na.rm=TRUE)/vxbartz_denom * FPC_vec2

          Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=T)

          # Taylor series expansion (Mood et al. 1974)
          prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
            colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)

          # ## estimate xbartz using my ugly derivation
          # vsumNtzx <- colSums((Ntz^2)*vxbartz + (xbartz^2)*vNtz - vNtz*vxbartz)
          # #vNz is calculated above
          # covthing <- colSums(xbartz*vNtz)
          # vxbarz <- ((colSums(Ntz*xbartz)/Nz)^2)*(vsumNtzx/(colSums(Ntz*xbartz)^2) + vNz/(Nz^2) - 2*covthing/(colSums(Ntz*xbartz)*Nz))

          out$n_length <- tapply(!is.na(length), cats, sum)
          out$mn_length <- xbarz
          out$se_length <- sqrt(vxbarz)
          out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          out$max_length <- tapply(length, cats, max, na.rm=TRUE)
        }
      } else {
        # --- if there are NO proportions to estimate ---
        nt <- table(stratum, is.na(length))[,1]    ### now excludes NA values in length
        if((is.null(Nhat) | !is.null(stratum_weights)) | FPC=="never") {
          if(verbose) cat("\n","Finite Population Correction factor NOT used for means", "\n")
          FPC_vec <- 1  # this will ignore FPC in further calculations
        } else {
          if(verbose) cat("\n","Finite Population Correction factor USED for means", "\n")
          FPC_vec <- (Nt-nt)/(Nt-1)
        }
        if(!is.null(stratum_weights)) {
          Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
        } else {
          Nt <- Nhat # consistency with report eqns
        }
        xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt*FPC_vec
        out$n_length <- sum(!is.na(length))
        out$mn_length <- sum(Nt*xbart/sum(Nt))
        out$se_length <- sqrt(sum(((Nt/sum(Nt))^2)*vxbart))
        out$min_length <- min(length, na.rm=TRUE)
        out$max_length <- max(length, na.rm=TRUE)
      }
    } else {
      # -- Stratified with error in Nhats --
      # estimating proportions..
      # estimating abundance
      # summarizing length
      if(verbose) cat("\n", "WITH error in Nhat", "\n")
      Nt <- Nhat # consistency with report eqns
      # modifying these equations from the Jim Creek report

      if(length(unique(cats)) > 1) {
        # --- if there are proportions to estimate ---
        ntz <- table(stratum, cats, useNA="no")
        nt <- rowSums(ntz)  # this will be robust to the presence of NA
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          stop("Nhat must be supplied when se_Nhat is used.")
        } else {
          Nt <- Nhat # consistency with report eqns
          # FPC_vec <- (Nt-nt)/(Nt-1)
        }
        if(is.null(stratum_weights) & FPC=="always") {
          FPC_vec <- (Nt-nt)/(Nt-1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for proportions", "\n")
        } else {
          FPC_vec <- 1#rep(1, length(Nt))
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for proportions", "\n")
        }
        ptz <- ntz/nt
        vptz <- ptz*(1-ptz)/(nt-1)*FPC_vec
        Ntz <- Nt*ptz
        vNtz <- ((Nt^2)*vptz) + ((ptz^2)*(se_Nhat^2)) - (vptz*(se_Nhat^2)) ### Goodman
        Nz <- colSums(Ntz)
        vNz <- colSums(vNtz)
        pz <- Nz/sum(Nt)

        # delta method
        covNzSum <- colSums((se_Nhat^2)*ptz, na.rm=TRUE)

        vpz <- ((Nz/sum(Nz))^2)*(vNz/(Nz^2) + sum(se_Nhat^2)/((sum(Nz))^2) - 2*covNzSum/(Nz*sum(Nz)))
        out$phat <- pz
        out$se_phat <- sqrt(vpz)

        out$Nhat <- Nz
        out$se_Nhat <- sqrt(vNz)

        # summarizing length
        if(!is.null(length)) {
          xbartz <- tapply(length, list(stratum, cats), mean, na.rm=TRUE)
          vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          if(is.null(stratum_weights) & FPC=="always") {
            FPC_vec2 <- (Ntz-vxbartz_denom)/(Ntz-1)
            if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
          } else {
            FPC_vec2 <- 1
            if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
          }
          vxbartz <- tapply(length, list(stratum, cats), var, na.rm=TRUE)/vxbartz_denom * FPC_vec2

          Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=TRUE)
          # Taylor series expansion (Mood et al. 1974)
          prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
            colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)
          out$n_length <- tapply(!is.na(length), cats, sum)
          out$mn_length <- xbarz
          out$se_length <- sqrt(vxbarz)
          out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          out$max_length <- tapply(length, cats, max, na.rm=TRUE)
        }
      } else {
        # --- if there are NO proportions to estimate ---
        nt <- table(stratum, is.na(length))[,1] # table(stratum)   ### now excludes NA values in length
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          stop("Nhat must be supplied when se_Nhat is used.")
        }
        if(is.null(stratum_weights) & FPC=="always") {
          FPC_vec <- (Nt-nt)/(Nt-1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
        } else {
          FPC_vec <- 1#rep(1, length(Nt))
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
        }
        Nt <- Nhat # consistency with report eqns
        xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt*FPC_vec
        out$n_length <- sum(nt)
        out$mn_length <- sum(Nt*xbart/sum(Nt))
        # NoverSumN <- Nt/sum(Nt)
        # # from delta method  -- validated by simulation   # currently no fpc anywhere
        # vNoverSumN <- (NoverSumN^2)*(((se_Nhat^2)/(Nt^2)) +
        #                                (sum(se_Nhat^2)/(sum(Nt)^2)) -
        #                                (2*(se_Nhat^2)/(Nt*sum(Nt))))
        # # from Goodman
        # out$se_length <- sqrt(sum(((NoverSumN^2)*vxbart) +   # vxbart validated by simulation
        #                             ((xbart^2)*vNoverSumN) -
        #                             (vxbart*vNoverSumN))) #### there is huge covariance in the sum!!
        # from covariance identities
        covSumNtSum <- sum(xbart*(se_Nhat^2))
        # from Goodman
        vSumNtXbart <- sum(((Nt^2)*vxbart) +
                             (xbart^2)*(se_Nhat^2) -
                             ((se_Nhat^2)*vxbart))
        # from delta method
        out$se_length <- sqrt(((sum(Nhat*xbart)/sum(Nhat))^2)*((vSumNtXbart/((sum(Nhat*xbart))^2)) +
                                                                 (sum(se_Nhat^2)/((sum(Nhat))^2)) -
                                                                 ((2*covSumNtSum)/(sum(Nhat)*sum(Nhat*xbart)))))
        out$min_length <- min(length, na.rm=TRUE)
        out$max_length <- max(length, na.rm=TRUE)
      }
    }
  }
  # --------- NOT Stratified---------
  if(is.null(stratum)) {
    if(verbose) cat("\n", "not stratified", "\n")

    # checking that inputs make sense
    if(length(Nhat) > 1 | length(se_Nhat) > 1 | length(stratum_weights) > 1) {
      stop("Stratum totals or weights are given, but no strata")
    }

    # estimating proportions..
    if(length(unique(cats)) > 1) {

      if(!is.null(Nhat) & (FPC=="always" | (FPC!="never" & is.null(se_Nhat)))) {
        FPC_prop <- (Nhat-sum(out$n))/(Nhat-1)
        if(verbose) cat("\n", "Finite Population Correction factor USED for proportions", "\n")
      } else {
        FPC_prop <- 1
        if(verbose) cat("\n", "Finite Population Correction factor NOT used for proportions", "\n")
      }

      out_ntot <- sum(out$n, na.rm=TRUE)      ########### this should have a na.rm=TRUE
      out$phat <- out$n/out_ntot
      out$se_phat <- sqrt(out$phat*(1-out$phat)/(out_ntot-1)*FPC_prop)

      # estimating abundance..
      if(!is.null(Nhat)) {

        out$Nhat <- Nhat * out$phat
        if(is.null(se_Nhat)) {
          if(verbose) cat("\n", "without error in Nhat", "\n")
          out$se_Nhat <- Nhat * out$se_phat   # verify this??
        } else {
          # Goodman's formula
          if(verbose) cat("\n", "with error in Nhat", "\n")
          out$se_Nhat <- sqrt(((Nhat^2)*(out$se_phat^2)) +
                                ((out$phat^2)*(se_Nhat^2)) -
                                ((out$se_phat^2)*(se_Nhat^2)))
        }
      }
    }

    # and summarizing lengths..
    if(!is.null(length)) {

      # populating FPC
      if(length(unique(cats)) > 1) {  # if there are proportions
        if(!is.null(Nhat) & FPC=="always") {
          FPC_mn <- (out$Nhat - out$n)/(out$Nhat - 1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
        } else {
          FPC_mn <- 1
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
        }
      } else {  # no proportions
        if(!is.null(Nhat) & ((FPC!="never" & is.null(se_Nhat)) | (FPC=="always"))) {
          FPC_mn <- (Nhat-sum(out$n))/(Nhat-1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
        } else {
          FPC_mn <- 1
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
        }
      }

      n_notNA <- tapply(!is.na(length), cats, sum)
      out$n_length <- n_notNA
      out$mn_length <- tapply(length, cats, mean, na.rm=TRUE)
      out_sd_length <- tapply(length, cats, sd, na.rm=TRUE) * sqrt(FPC_mn)
      out$se_length <- out_sd_length/sqrt(n_notNA)
      out$min_length <- tapply(length, cats, min, na.rm=TRUE)
      out$max_length <- tapply(length, cats, max, na.rm=TRUE)
    }
  }
  return(out)
}

# # simulating some fake data..
# nn <- 100
# sex <- sample(c("F","M"), nn, replace=T)
# age <- sample(c(1.1,1.2,1.3,2.1,2.2,2.3), nn, replace=T)
# length <- rnorm(nn, mean=700, sd=100)
#
# length[1] <- NA
#
# stratum <- sample(1:5, nn, replace=T)
# Nhat <- c(1000, 2000, 3000, 4000, 5000)
# ASL_table(
#   age=age,
#   # sex=sex,
#   length=length,
#   # Nhat=10000,
#   # se_Nhat=1000,
#   Nhat=Nhat,
#   se_Nhat=500,
#   # se_Nhat=Nhat/10,
#   # stratum_weights=Nhat/10,
#   stratum=stratum
# )




# verify_ASL_table <- function(case=NULL,   # should this default to NULL?
#                              nstrata,
#                              nage,
#                              nt,   # c(100, 100, 100, 100), # sample size for each stratum
#                              Nt,   # c(10000, 20000, 30000, 40000),  # abundance for each stratum
#                              stratum_weights = NULL,
#                              se_Nt = NULL,   # 0.2*Nt, # (possible) SE for abundance by stratum,
#                              mn_length = NULL,   # c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
#                              sd_length = NULL,   # c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
#                              ptz = NULL,   # matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                                           # c(1,2,5,5,2),
#                                           # c(2,5,3,2,1),
#                                           # c(5,4,3,1,1)),
#                                           # byrow=TRUE,
#                                           # nrow=4,
#                                           # ncol=5),
#                              nNA=10,
#                              nsim=1000,   # 1000,    # number of simulated replicates
#                              plot_pop=TRUE,   # whether to make summary plots of pop & one sample
#                              verbose=TRUE, # whether to output cases in sim function
#                              print_table=FALSE) {
#
#   ### pre-populate cases
#   if(!is.null(case)) {
#     if(case=="stratified_witherror_lengthage") {
#       nstrata <- 4
#       nage <- 5
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
#       se_Nt <- 0.2*c(10000, 20000, 30000, 40000) # (possible) SE for abundance by stratum,
#       mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
#       sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
#       ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                      c(1,2,5,5,2),
#                      c(2,5,3,2,1),
#                      c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 5,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
# se_Nt = 0.2*c(10000, 20000, 30000, 40000), # (possible) SE for abundance by stratum
# mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
# sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
# ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                c(1,2,5,5,2),
#                c(2,5,3,2,1),
#                c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#
# ")
#       }
#     }
#     if(case=="stratified_witherror_age") {
#       nstrata <- 4
#       nage <- 5
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
#       se_Nt <- 0.2*c(10000, 20000, 30000, 40000) # (possible) SE for abundance by stratum,
#       mn_length <- NULL # mean length FOR EACH AGE
#       sd_length <- NULL # sd length FOR EACH AGE
#       ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                       c(1,2,5,5,2),
#                       c(2,5,3,2,1),
#                       c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#       # mn_length=NULL, sd_length=NULL
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 5,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
# se_Nt = 0.2*c(10000, 20000, 30000, 40000), # (possible) SE for abundance by stratum
# mn_length = NULL, # mean length FOR EACH AGE
# sd_length = NULL, # sd length FOR EACH AGE
# ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                 c(1,2,5,5,2),
#                 c(2,5,3,2,1),
#                 c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
# ")
#       }
#     }
#     if(case=="stratified_witherror_length") {
#       nstrata <- 4
#       nage <- 1
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
#       se_Nt <- 0.2*c(10000, 20000, 30000, 40000) # (possible) SE for abundance by stratum,
#       mn_length <- 300 # mean length FOR EACH AGE
#       sd_length <- 70 # sd length FOR EACH AGE
#       ptz <- NULL
#       # mn_length=300, sd_length=70, ptz=NULL
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 1,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
# se_Nt = 0.2*c(10000, 20000, 30000, 40000), # (possible) SE for abundance by stratum
# mn_length = 300, # mean length FOR EACH AGE
# sd_length = 70, # sd length FOR EACH AGE
# ptz = NULL
# ")
#       }
#     }
#     if(case=="stratified_lengthage") {
#       nstrata <- 4
#       nage <- 5
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
#       sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
#       ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                       c(1,2,5,5,2),
#                       c(2,5,3,2,1),
#                       c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 5,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
# sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
# ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                 c(1,2,5,5,2),
#                 c(2,5,3,2,1),
#                 c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
# ")
#       }
#       # se_Nt = NULL
#     }
#     if(case=="stratified_age") {
#       nstrata <- 4
#       nage <- 5
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- NULL # mean length FOR EACH AGE
#       sd_length <- NULL # sd length FOR EACH AGE
#       ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                       c(1,2,5,5,2),
#                       c(2,5,3,2,1),
#                       c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#       # se_Nt = NULL, mn_length=NULL, sd_length=NULL
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 5,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = NULL, # mean length FOR EACH AGE
# sd_length = NULL, # sd length FOR EACH AGE
# ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                 c(1,2,5,5,2),
#                 c(2,5,3,2,1),
#                 c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
# ")
#       }
#     }
#     if(case=="stratified_length") {
#       nstrata <- 4
#       nage <- 1
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- 300 # mean length FOR EACH AGE
#       sd_length <- 70 # sd length FOR EACH AGE
#       ptz <- NULL
#       # se_Nt = NULL, mn_length=300, sd_length=70, ptz=NULL
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 1,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = 300, # mean length FOR EACH AGE
# sd_length = 70, # sd length FOR EACH AGE
# ptz = NULL
# ")
#       }
#     }
#     if(case=="pooled_witherror_lengthage") {
#       nstrata <- 1
#       nage <- 5
#       nt <- 400 # sample size for each stratum
#       Nt <- 10000  # abundance for each stratum
#       se_Nt <- 0.2*10000 # (possible) SE for abundance by stratum,
#       mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
#       sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
#       ptz <- 1:5
#       # ptz=1:5, Nt=10000, nt=100
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 5,
# nt = 400, # sample size for each stratum
# Nt = 10000,  # abundance for each stratum
# se_Nt = 0.2*10000, # (possible) SE for abundance by stratum
# mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
# sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
# ptz = 1:5
# ")
#       }
#     }
#     if(case=="pooled_witherror_age") {
#       nstrata <- 1
#       nage <- 5
#       nt <- 400 # sample size for each stratum
#       Nt <- 10000  # abundance for each stratum
#       se_Nt <- 0.2*10000 # (possible) SE for abundance by stratum,
#       mn_length <- NULL # mean length FOR EACH AGE
#       sd_length <- NULL # sd length FOR EACH AGE
#       ptz <- 1:5
#       # mn_length=NULL, sd_length=NULL, ptz=1:5, Nt=10000, nt=100
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 5,
# nt = 400, # sample size for each stratum
# Nt = 10000,  # abundance for each stratum
# se_Nt = 0.2*10000, # (possible) SE for abundance by stratum
# mn_length = NULL, # mean length FOR EACH AGE
# sd_length = NULL, # sd length FOR EACH AGE
# ptz = 1:5
# ")
#       }
#     }
#     if(case=="pooled_witherror_length") {
#       nstrata <- 1
#       nage <- 1
#       nt <- 400 # sample size for each stratum
#       Nt <- 10000  # abundance for each stratum
#       se_Nt <- 0.2*10000 # (possible) SE for abundance by stratum,
#       mn_length <- 300 # mean length FOR EACH AGE
#       sd_length <- 70 # sd length FOR EACH AGE
#       ptz <- NULL
#       # mn_length=300, sd_length=70, ptz=NULL, Nt=10000, nt=100
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 1,
# nt = 400, # sample size for each stratum
# Nt = 10000,  # abundance for each stratum
# se_Nt = 0.2*10000, # (possible) SE for abundance by stratum
# mn_length = 300, # mean length FOR EACH AGE
# sd_length = 70, # sd length FOR EACH AGE
# ptz = NULL
# ")
#       }
#     }
#     if(case=="pooled_lengthage") {
#       nstrata <- 1
#       nage <- 5
#       nt <- 400 # sample size for each stratum
#       Nt <- 10000  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
#       sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
#       ptz <- 1:5
#       # se_Nt = NULL, ptz=1:5, Nt=10000, nt=100
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 5,
# nt = 400, # sample size for each stratum
# Nt = 10000,  # abundance for each stratum
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
# sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
# ptz = 1:5
# ")
#       }
#     }
#     if(case=="pooled_age") {
#       nstrata <- 1
#       nage <- 5
#       nt <- 400 # sample size for each stratum
#       Nt <- 10000  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- NULL # mean length FOR EACH AGE
#       sd_length <- NULL # sd length FOR EACH AGE
#       ptz <- 1:5
#       # se_Nt = NULL, mn_length=NULL, sd_length=NULL, ptz=1:5, Nt=10000, nt=100
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 5,
# nt = 400, # sample size for each stratum
# Nt = 10000,  # abundance for each stratum
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = NULL, # mean length FOR EACH AGE
# sd_length = NULL, # sd length FOR EACH AGE
# ptz = 1:5
# ")
#       }
#     }
#     if(case=="pooled_length") {
#       nstrata <- 1
#       nage <- 1
#       nt <- 400 # sample size for each stratum
#       Nt <- 10000 # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- 300 # mean length FOR EACH AGE
#       sd_length <- 70 # sd length FOR EACH AGE
#       ptz <- NULL
#       # se_Nt = NULL, mn_length=300, sd_length=70, ptz=NULL, Nt=10000, nt=100
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 1,
# nt = 400, # sample size for each stratum
# Nt = 10000, # abundance for each stratum
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = 300, # mean length FOR EACH AGE
# sd_length = 70, # sd length FOR EACH AGE
# ptz = NULL
# ")
#       }
#     }
#     if(case=="stratified_Nunknown_lengthage") {
#       nstrata <- 4
#       nage <- 5
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- NULL
#       stratum_weights <- 1:4  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
#       sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
#       ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                       c(1,2,5,5,2),
#                       c(2,5,3,2,1),
#                       c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 5,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = NULL,  # abundance for each stratum
# stratum_weights = 1:4, # stratum weights instead of abundance
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
# sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
# ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                c(1,2,5,5,2),
#                c(2,5,3,2,1),
#                c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#
# ")
#       }
#     }
#     if(case=="stratified_Nunknown_length") {
#       nstrata <- 4
#       nage <- 1
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- NULL
#       stratum_weights <- 1:4  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- 300 # mean length FOR EACH AGE
#       sd_length <- 70 # sd length FOR EACH AGE
#       ptz <- NULL  # matrix of probabilities of each age BY stratum
#
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 1,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = NULL,  # abundance for each stratum
# stratum_weights = 1:4, # stratum weights instead of abundance
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = 300, # mean length FOR EACH AGE
# sd_length = 70, # sd length FOR EACH AGE
# ptz = NULL  # matrix of probabilities of each age BY stratum
#
# ")
#       }
#     }
#     if(case=="stratified_Nunknown_age") {
#       nstrata <- 4
#       nage <- 5
#       nt <- c(100, 100, 100, 100) # sample size for each stratum
#       Nt <- NULL
#       stratum_weights <- 1:4  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- NULL # mean length FOR EACH AGE
#       sd_length <- NULL # sd length FOR EACH AGE
#       ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                       c(1,2,5,5,2),
#                       c(2,5,3,2,1),
#                       c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#       if(verbose) {
#         cat("
# nstrata = 4,
# nage = 5,
# nt = c(100, 100, 100, 100), # sample size for each stratum
# Nt = NULL,  # abundance for each stratum
# stratum_weights = 1:4, # stratum weights instead of abundance
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = NULL, # mean length FOR EACH AGE
# sd_length = NULL, # sd length FOR EACH AGE
# ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
#                c(1,2,5,5,2),
#                c(2,5,3,2,1),
#                c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
#
# ")
#       }
#     }
#     if(case=="pooled_Nunknown_lengthage") {
#       nstrata <- 1
#       nage <- 5
#       nt <- 400 # sample size for each stratum
#       Nt <- NULL
#       stratum_weights <- 1  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
#       sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
#       ptz <- 1:5
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 5,
# nt = 400, # sample size for each stratum
# Nt = NULL,  # abundance for each stratum
# stratum_weights = 1, # stratum weights instead of abundance
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
# sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
# ptz = 1:5
#
# ")
#       }
#     }
#     if(case=="pooled_Nunknown_length") {
#       nstrata <- 1
#       nage <- 1
#       nt <- 400 # sample size for each stratum
#       Nt <- NULL
#       stratum_weights <- 1  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- 300 # mean length FOR EACH AGE
#       sd_length <- 70 # sd length FOR EACH AGE
#       ptz <- NULL  # matrix of probabilities of each age BY stratum
#
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 1,
# nt = 400, # sample size for each stratum
# Nt = NULL,  # abundance for each stratum
# stratum_weights = 1, # stratum weights instead of abundance
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = 300, # mean length FOR EACH AGE
# sd_length = 70, # sd length FOR EACH AGE
# ptz = NULL  # matrix of probabilities of each age BY stratum
#
# ")
#       }
#     }
#     if(case=="pooled_Nunknown_age") {
#       nstrata <- 1
#       nage <- 5
#       nt <- 400 # sample size for each stratum
#       Nt <- NULL
#       stratum_weights <- 1  # abundance for each stratum
#       se_Nt <- NULL # (possible) SE for abundance by stratum,
#       mn_length <- NULL # mean length FOR EACH AGE
#       sd_length <- NULL # sd length FOR EACH AGE
#       ptz <- 1:5
#       if(verbose) {
#         cat("
# nstrata = 1,
# nage = 5,
# nt = 400, # sample size for each stratum
# Nt = NULL,  # abundance for each stratum
# stratum_weights = 1, # stratum weights instead of abundance
# se_Nt = NULL, # (possible) SE for abundance by stratum
# mn_length = NULL, # mean length FOR EACH AGE
# sd_length = NULL, # sd length FOR EACH AGE
# ptz = 1:5
#
# ")
#       }
#     }
#   }
#
#   # pre-fixing ptz
#   if(is.null(dim(ptz)) & !is.null(ptz)) {
#     ptz <- t(as.matrix(ptz))
#   }
#
#   # check inputs: length(nt) vs length(Nt) vs. length(se_Nt) vs nrow(ptz)
#   if(length(nt)!=nstrata) stop("length(nt) should equal nstrata.")
#   if(length(Nt)!=nstrata & length(stratum_weights)!=nstrata) stop("length(Nt) or length(stratum_weights) should equal nstrata.")
#   if(!is.null(se_Nt) & (length(se_Nt)!=nstrata)) stop("length(se_Nt) should equal nstrata.")
#   if(!is.null(ptz)) {
#     if(nrow(ptz)!=nstrata) stop("nrow(ptz) should equal nstrata.")
#   }
#
#   # check inputs: length(mn_length) vs length(sd_length) vs ncol(ptz)
#   if(!is.null(mn_length) & (length(mn_length)!=nage)) stop("length(mn_length) should equal nage.")
#   if(!is.null(sd_length) & (length(sd_length)!=nage)) stop("length(sd_length) should equal nage.")
#   if(length(sd_length) != length(mn_length)) stop("mn_length and sd_length should have corresponding length.")
#   if(!is.null(ptz)) {
#     if(nage>1 & ncol(ptz)!=nage) stop("ncol(ptz) or length(ptz) must equal nage.")
#   }
#
#   #### no age means mn_length and sd_length have length 1, and ptz will be NULL
#   #### no length means mn_length and sd_length will be NULL
#   #### no strata means nt and Nt will be length 1 and ptz will be a vector
#
#   # nstrata <- length(nt)
#
#
#   ## constructing vectors at the POPULATION level
#   if(!is.null(Nt)) {
#   # stratum
#   t <- rep(1,Nt[1])
#   if(length(Nt) > 1) {
#     for(i in 2:length(Nt)) t <- c(t, rep(i, Nt[i]))
#   }
#
#   # age
#   # nage <- ifelse(!is.null(mn_length), length(mn_length), ncol(ptz))
#   if(!is.null(ptz)) {
#     age_sim <- NA*t
#     for(i_t in 1:nstrata) {
#       age_sim[t==i_t] <- sample(x=1:nage, size=sum(t==i_t), prob=ptz[i_t,],replace=TRUE)
#     }
#     # age_sim[t==1] <- sample(x=1:5, size=sum(t==1), prob=c(1,2,3,4,5), replace=T)
#   } else {
#     age_sim <- rep(1, length(t))
#   }
#
#   # length
#   if(!is.null(mn_length)) {
#     length_sim <- NA*t
#     for(i_z in 1:nage) {
#       length_sim[age_sim==i_z] <- rnorm(sum(age_sim==i_z), mn_length[i_z], sd_length[i_z])
#     }
#   } else {
#     length_sim <- NULL
#   }
#
#
#
#   # # storing graphics state to save
#   # parmfrow <- par("mfrow")
#   # on.exit(par(mfrow=parmfrow))  # making sure to re-set graphics state
#
#   # # plotting simulated population
#   # mosaicplot(table(t,age_sim), xlab="Stratum", ylab="Age")
#   # boxplot(length_sim~age_sim, xlab="Age", ylab="Length")
#   # boxplot(length_sim~t, xlab="Stratum", ylab="Length")
#
#
#
#   # simulate once to set up output array (yes this is a hack)
#   thesample <- sample(seq_along(t)[t==1], size=nt[1])
#   if(length(nt) > 1) {
#     for(i_t in 2:nstrata) {
#       thesample <- c(thesample, sample(seq_along(t)[t==i_t], size=nt[i_t]))
#     }
#   }
#   # table(t[thesample])  # making sure it worked
#   } else {
#     # t is stratum at pop level
#     # age from stratum
#     age_sim <- NULL
#     t <- NULL
#     for(i_nt in 1:length(nt)) {
#       if(!is.null(ptz)) {
#         age_sim <- c(age_sim, sample(1:nage, size=nt[i_nt], prob=ptz[i_nt, ], replace=TRUE))
#       } else {
#         age_sim <- c(age_sim, rep(1, nt[i_nt]))
#       }
#       t <- c(t, rep(i_nt, nt[i_nt]))
#     }
#
#     # length from age
#     if(!is.null(mn_length)) {
#       length_sim <- rnorm(length(t), mean=mn_length[age_sim], sd=sd_length[age_sim])
#     } else {
#       length_sim <- NULL
#     }
#
#     # dummy vector to maintain compatibility with existing code
#     thesample <- seq_along(t)
#   }
#
#   # plotting population & one sample (optionally)
#   if(plot_pop) {
#     if(prod(dim(table(t,age_sim))) > 1) {  #
#       if(!is.null(Nt)) mosaicplot(table(t,age_sim), xlab="Stratum", ylab="Age", main="Population", col=grey.colors(nage, rev=TRUE))
#       mosaicplot(table(t[thesample],age_sim[thesample]), xlab="Stratum", ylab="Age", main="Sample", col=grey.colors(nage, rev=TRUE))
#     }
#     if(!is.null(length_sim)) {
#       if(!is.null(Nt)) boxplot(length_sim~age_sim, xlab="Age", ylab="Length", main="Population")
#       boxplot(length_sim[thesample]~age_sim[thesample], xlab="Age", ylab="Length", main="Sample")
#       if(!is.null(Nt)) boxplot(length_sim~t, xlab="Stratum", ylab="Length", main="Population")
#       boxplot(length_sim[thesample]~t[thesample], xlab="Stratum", ylab="Length", main="Sample")
#     }
#   }
#
#   if(is.null(se_Nt)) {
#     Nhat_sim <- Nt
#   } else {
#     Nhat_sim <- rnorm(length(Nt), mean=Nt, sd=se_Nt) # could move this into function args
#   }
#
#   ### LEAVE IT OPEN FOR age_sim=NULL and length_sim=NULL
#   if(nage == 1) age_sim <- NULL
#   # if(nstrata == 1) t <- NULL
#
#   thestratum <- as.integer(t[thesample])
#   if(all(thestratum==1)) thestratum <- NULL
#
#   theage <- age_sim[thesample]
#   thelength <- length_sim[thesample]
#
#
#   ## imputing some NA values to make sure the function is robust to NA!!
#   if(!is.null(theage)) theage[sample(seq_along(theage), nNA)] <- NA
#   if(!is.null(thelength)) thelength[sample(seq_along(thelength), nNA)] <- NA
#   if(!is.null(thestratum)) thestratum[sample(seq_along(thestratum), nNA)] <- NA
#
#
#   thetable <- ASL_table(age=theage,
#                         length=thelength,
#                         stratum=thestratum,
#                         Nhat=Nhat_sim,
#                         stratum_weights = stratum_weights,
#                         se_Nhat=se_Nt,
#                         verbose=verbose)
#
#   # thetable <- ASL_table(age=age_sim[thesample],
#   #                       length=length_sim[thesample],
#   #                       stratum=thestratum,
#   #                       Nhat=Nhat_sim,
#   #                       se_Nhat=se_Nt,
#   #                       verbose=verbose)  # find a way to add stratum_weights???
#   if(print_table) print(thetable)
#
#   # initiate all stuff to store
#   # results <- array(dim=c(nrow(thetable), ncol(thetable), nsim))
#   results <- array(dim=c(nsim, nrow(thetable), ncol(thetable)))
#
#   ############### then loop nsim times!!! ###############
#   for(i_sim in 1:nsim) {
#     if(!is.null(Nt)) {
#     # first take a sample of fish
#     thesample <- sample(seq_along(t)[t==1], size=nt[1])
#     if(length(nt) > 1) {
#       for(i_t in 2:nstrata) {
#         thesample <- c(thesample, sample(seq_along(t)[t==i_t], size=nt[i_t]))
#       }
#     }
#     # table(t[thesample])  # making sure it worked
#     } else {
#       # t is stratum at pop level
#       # age from stratum
#       age_sim <- NULL
#       t <- NULL
#       for(i_nt in 1:length(nt)) {
#         if(!is.null(ptz)) {
#           age_sim <- c(age_sim, sample(1:nage, size=nt[i_nt], prob=ptz[i_nt, ], replace=TRUE))
#         } else {
#           age_sim <- c(age_sim, rep(1, nt[i_nt]))
#         }
#         t <- c(t, rep(i_nt, nt[i_nt]))
#       }
#
#       # length from age
#       if(!is.null(mn_length)) {
#         length_sim <- rnorm(length(t), mean=mn_length[age_sim], sd=sd_length[age_sim])
#       } else {
#         length_sim <- NULL
#       }
#
#       # dummy vector to maintain compatibility with existing code
#       thesample <- seq_along(t)
#       if(nage==1) age_sim <- NULL
#     }
#
#     if(is.null(se_Nt)) {
#       Nhat_sim <- Nt
#     } else {
#       Nhat_sim <- rnorm(length(Nt), mean=Nt, sd=se_Nt) # could move this into function args
#     }
#
#     thestratum <- as.integer(t[thesample])    #### this is new
#     if(all(thestratum==1)) thestratum <- NULL    #### this is new
#
#     theage <- age_sim[thesample]
#     thelength <- length_sim[thesample]
#
#
#     ## imputing some NA values to make sure the function is robust to NA!!
#     if(!is.null(theage)) theage[sample(seq_along(theage), nNA)] <- NA
#     if(!is.null(thelength)) thelength[sample(seq_along(thelength), nNA)] <- NA
#     if(!is.null(thestratum)) thestratum[sample(seq_along(thestratum), nNA)] <- NA
#
#
#     thetable <- ASL_table(age=theage,
#                           length=thelength,
#                           stratum=thestratum,
#                           Nhat=Nhat_sim,
#                           stratum_weights=stratum_weights,
#                           se_Nhat=se_Nt)  # find a way to add stratum_weights???
#     # results[,,i_sim] <- as.matrix(thetable)
#     results[i_sim,,] <- as.matrix(thetable)
#   }
#   # dimnames(results)[1:2] <- dimnames(thetable)
#   dimnames(results)[2:3] <- dimnames(thetable)
#
#   ### plotting simulation results, overlayed with true values
#   if("phat" %in% dimnames(results)[[3]]) {
#     phat_sim <- results[,,dimnames(results)[[3]]=="phat"]
#     if(!is.null(Nt)) {
#       phat_true <- as.numeric(table(age_sim)/sum(Nt))
#     } else {
#       phat_true <- colSums(ptz*stratum_weights)/sum(ptz*stratum_weights)
#     }
#
#     boxplot(phat_sim, ylim=range(phat_true, phat_sim, na.rm=TRUE),
#             main="p_hat", xlab="Age",
#             border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
#     points(phat_true, col=2, pch=16)
#     legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
#            fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
#
#     se_phat_sim <- results[,,dimnames(results)[[3]]=="se_phat"]
#     se_phat_true <- apply(phat_sim, 2, sd, na.rm=TRUE)
#     boxplot(se_phat_sim, ylim=range(se_phat_true, se_phat_sim, na.rm=TRUE),
#             main="se(p_hat)", xlab="Age",
#             border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
#     points(se_phat_true, col=2, pch=16)
#     legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
#            fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
#   }
#
#   if("Nhat" %in% dimnames(results)[[3]]) {
#     Nhat_sim <- results[,,dimnames(results)[[3]]=="Nhat"]
#     Nhat_true <- as.numeric(table(age_sim))
#     boxplot(Nhat_sim, ylim=range(Nhat_true, Nhat_sim, na.rm=TRUE),
#             main="N_hat", xlab="Age",
#             border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
#     points(Nhat_true, col=2, pch=16)
#     legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
#            fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
#
#     se_Nhat_sim <- results[,,dimnames(results)[[3]]=="se_Nhat"]
#     se_Nhat_true <- apply(Nhat_sim, 2, sd, na.rm=TRUE)
#     boxplot(se_Nhat_sim, ylim=range(se_Nhat_true, se_Nhat_sim, na.rm=TRUE),
#             main="se(N_hat)", xlab="Age",
#             border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
#     points(se_Nhat_true, col=2, pch=16)
#     legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
#            fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
#   }
#
#   if("mn_length" %in% dimnames(results)[[3]]) {
#     mn_length_sim <- as.matrix(results[,,dimnames(results)[[3]]=="mn_length"])
#     mn_length_true <- mn_length
#     boxplot(mn_length_sim, ylim=range(mn_length_true, mn_length_sim, na.rm=TRUE),
#             main="mn_length", xlab="Age",
#             border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
#     points(mn_length_true, col=2, pch=16)
#     legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
#            fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
#
#     se_mn_length_sim <- results[,,dimnames(results)[[3]]=="se_length"]
#     se_mn_length_true <- apply(mn_length_sim, 2, sd, na.rm=TRUE)
#     boxplot(se_mn_length_sim, ylim=range(se_mn_length_true, se_mn_length_sim, na.rm=TRUE),
#             main="se(mn_length)", xlab="Age",
#             border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
#     points(se_mn_length_true, col=2, pch=16)
#     legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
#            fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
#   }
# }







#' Example data: sim_data
#'
#' A simulated dataset intended to illustrate \link{ASL_table}.  This dataset
#' is formatted as a list with two components:
#' * `sim_data$data` is a data frame with 400 observations of `$age`, `$length`,
#' and `$stratum`.
#' * `sim_data$abundance` is a data frame with 4 rows of `$Nhat` and `$se_Nhat`,
#' corresponding to each stratum.
"sim_data"






#' ASL Boilerplate
#' @description Automagically creates boilerplate text to accompany an ASL
#' summary analysis, or methods for an Operational Plan.
#'
#' This can be done two ways:
#' * (1) Directly from arguments `stratified=`, `abundance=`, and `data=` (see argument description below)
#' * (2) From possible inputs to \link{ASL_table}.  This method is provided for convenience.
#'
#' Text and equations are taken directly from the document found here:
#'
#' `https://github.com/ADFG-DSF/dsftools/blob/main/addl_documentation/asl_equations.Rmd`
#' @param stratified Whether stratified estimators will be used (`TRUE`) or not (`FALSE`)
#' @param abundance Possible values of `c("known", "estimated", "unknown")`,
#' defining how abundance is treated.
#' @param data A vector containing some subset of `c("age", "sex", "length", "length_cat")`,
#' depending on what data are present
#' @param species Optional character to use for species.  Defaults to `"fish"`.
#' @param tense Which verb tense to use.  If the default `NA` is accepted, this
#' will be imputed as `"future"` if text is constructed from method 1 above,
#' or `"past"` if text is constructed from data (method 2 above).  Possible
#' values are `c(NA, "past", "present", "future")`
#' @param age Argument to `ASL_table()`.  Defaults to `NULL`.
#' @param sex Argument to `ASL_table()`.  Defaults to `NULL`.
#' @param length Argument to `ASL_table()`.  Defaults to `NULL`.
#' @param length_cat Argument to `ASL_table()`.  Defaults to `NULL`.
#' @param stratum Argument to `ASL_table()`.  Defaults to `NULL`.
#' @param Nhat Argument to `ASL_table()`.  Defaults to `NULL`.
#' @param se_Nhat Argument to `ASL_table()`.  Defaults to `NULL`.
#' @param stratum_weights Argument to `ASL_table()`.  Not used.
#' @param verbose Argument to `ASL_table()`.  Not used.
#' @param FPC Argument to `ASL_table()`.  Possible arguments are `c("ifknown",
#' "always", "never")`.  Defaults to `"ifknown"`.
#' @return NULL
#' @seealso \link{ASL_table}
#' @author Matt Tyers
#' @examples
#' \dontrun{
#'
#' ## An example Rmarkdown code chunk might look like the following, if
#' ## text is to be created from direct input
#'
#' ```{r, results='asis', echo=FALSE}
#' library(dsftools)
#' ASL_boilerplate(data=c("age", "length"), stratified=TRUE, abundance="known")
#' ```
#'
#' ## Alternately, an Rmarkdown code chunk might look like the following, if
#' ## text is to be created from data
#'
#' ```{r, results='asis', echo=FALSE}
#' library(dsftools)
#' ASL_boilerplate(age = sim_data$data$age,
#'                 length = sim_data$data$length,
#'                 stratum = sim_data$data$stratum,
#'                 Nhat = sim_data$abundance$Nhat)
#' ```
#' }
#' @importFrom utils citation
#' @importFrom utils capture.output
#' @export
ASL_boilerplate <- function(stratified=NULL,   # logical TRUE or FALSE
                             abundance=NULL,   # c("known", "estimated", "unknown")
                             data=NULL,   # c("age","sex","length")
                             species="fish",
                             tense=c(NA,"past","present","future"),

                             age=NULL,
                             sex=NULL,
                             length=NULL,
                             length_cat=NULL,
                             stratum=NULL,
                             Nhat=NULL,
                             se_Nhat=NULL,
                             stratum_weights = NULL,
                             verbose=FALSE,
                             FPC=c("ifknown", "always", "never")) {

  if(!is.null(stratified) | !is.null(abundance) | !is.null(data)) {
    if(is.null(stratified) | is.null(abundance) | is.null(data)) {
      stop("need inputs to stratified=, abundance=, and data=")
    }
    fromData <- FALSE
  } else {
    fromData <- TRUE
    if(is.null(age) & is.null(sex) & is.null(length) & is.null(length_cat)) {
      stop("need inputs to at least one of age=, sex=, length=, or length_cat=")
    }
  }


  if(fromData) {
    stratified <- !is.null(stratum)
    abundance <- ifelse(is.null(Nhat), "unknown",
                        ifelse(!is.null(se_Nhat), "estimated", "known"))
  }
  FPC <- match.arg(FPC)



  ### constructing some hidden switches
  if(fromData) {
    doAge <- !is.null(age)
    doSex <- !is.null(sex)
    doLength <- !is.null(length)
    doLengthcat <- !is.null(length_cat)
  } else {
    doAge <- "age" %in% data
    doSex <- "sex" %in% data
    doLength <- "length" %in% data
    doLengthcat <- "length_cat" %in% data
  }

  tense <- match.arg(tense)
  if(is.na(tense)) {
    if(!fromData) {
      tense <- "future"
    } else {
      tense <- "past"
    }
  }

  if(tense == "past") {
    verb1 <- "was"
    verb2 <- "was then"
    verb3 <- "were"
    verb4 <- "were then"
  }
  if(tense == "present") {
    verb1 <- "is"
    verb2 <- "is then"
    verb3 <- "are"
    verb4 <- "are then"
  }
  if(tense == "future") {
    verb1 <- "will be"
    verb2 <- "will then be"
    verb3 <- "will be"
    verb4 <- "will then be"
  }

  doProp <- (doAge | doSex| doLengthcat)  # if there are proportions

  if(doProp) {   # if there are proportions
    howmany <- doAge + doSex + doLengthcat
    if(howmany==3) {
      agesex <- "age, sex, and length"
    }
    if(howmany==2) {
      if(doAge & doSex) agesex <- "age and sex"
      if(doAge & doLengthcat) agesex <- "age and length"
      if(doSex & doLengthcat) agesex <- "sex and length"
    }
    if(howmany==1) {
      if(doAge) agesex <- "age"
      if(doSex) agesex <- "sex"
      if(doLengthcat) agesex <- "length"
    }
  }

  if(abundance == "known") {
    doFPCprop <- (FPC != "never")
    doFPCmean <- (FPC == "always")
    if(!doAge & !doSex) {
      doFPCmean <- (FPC != "never")
    }
  }
  if(abundance == "estimated") {
    doFPCprop <- (FPC == "always")
    doFPCmean <- (FPC == "always")
  }
  if(abundance == "unknown") {
    doFPCprop <- FALSE
    doFPCmean <- FALSE
  }

  # initializing reference switches, these will be turned to TRUE if needed
  doCasella <- FALSE
  doCochran <- doProp
  doGoodman <- (abundance == "estimated")
  doMood <- FALSE


  ### ---------- filling in boilerplate text!!! ---------- ###

  if(stratified) {
    if(abundance == "known") {
      if(doProp) {   # if there are proportions
        cat("The proportion of each ",
            agesex,
            " category *z* ", verb1, " estimated for each sampling stratum *t* as follows:

$$\\hat{p}_{tz}=\\frac{n_{tz}}{n_t}$$

in which $n_{tz}$ equals the number of ",
            species,
            " sampled during sampling stratum $t$ classified as ",
            agesex,
            " category $z$, and $n_t$ equals the number of ",
            species,
            " sampled for ",
            agesex,
            " determination within sampling stratum $t$.

The sampling variance of $\\hat{p}_{tz}$ ", verb1, " estimated as the following (Cochran 1977)")
        if(doFPCprop) {
          cat(", in which $N_t$ represents the total abundance of ",
              species,
              " in sampling stratum *t*:

$$\\hat{var}[\\hat{p}_{tz}]=\\frac{\\hat{p}_{tz}(1-\\hat{p}_{tz})}{n_t-1}\\left(\\frac{N_t-n_t}{N_t-1}\\right)$$
")
        } else {
          cat(":

$$\\hat{var}[\\hat{p}_{tz}]=\\frac{\\hat{p}_{tz}(1-\\hat{p}_{tz})}{n_t-1}$$
")
        }
        cat("
The total abundance by ",
            agesex,
            " category in each sampling stratum ", verb1, " estimated as follows")
        if(!doFPCprop) cat(", in which $N_t$ represents the total abundance of ",species," in sampling stratum *t*")
        cat(":

$$\\hat{N}_{tz}=N_t\\hat{p}_{tz}$$

with variance estimated as

$$\\hat{var}[\\hat{N}_{tz}]=N_t^2\\hat{var}[\\hat{p}_{tz}]$$

The total abundance by ",
            agesex,
            " category and its variance ", verb4, " estimated by summation as follows:

$$\\hat{N}_z=\\sum_{t=1}^{L}\\hat{N}_{tz}$$

and

$$\\hat{var}[\\hat{N}_{z}]=\\sum_{t=1}^{L}\\hat{var}[\\hat{N}_{tz}]$$

where $L$ equals the number of sampling strata.

Finally, the overall proportion by ",
            agesex,
            " category and its variance ", verb3, " estimated as follows:

$$\\hat{p}_z=\\frac{\\hat{N}_z}{N}$$

and

$$\\hat{var}[\\hat{p}_z]=\\frac{\\hat{var}[\\hat{N}_z]}{N^2}$$

where $N$ is the total abundance across all sampling periods.
")
        if(doLength) {
          doMood <- TRUE
          cat("
The mean length by ",
              agesex,
              " for each sampling stratum ", verb1, " estimated as follows:

$$\\bar{x}_{tz}=\\frac{\\sum_{i=1}^{n_{tz}}x_{tzi}}{n_{tz}}$$

where $x_{tzi}$ is the length of the *i*th ",
              species,
              " sampled of ",
              agesex,
              " category $z$ during sampling stratum $t$.

The sampling variance of $\\bar{x}_{tz}$ ", verb1, " estimated as
")
          if(doFPCmean) {
            cat("
$$\\hat{var}[\\bar{x}_{tz}]=\\frac{\\sum_{i=1}^{n_{tz}}(x_{tzi}-\\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}\\left(\\frac{\\hat{N}_{tz}-n_{tz}}{\\hat{N}_{tz}-1}\\right)$$
")
          } else {
            cat("
$$\\hat{var}[\\bar{x}_{tz}]=\\frac{\\sum_{i=1}^{n_{tz}}(x_{tzi}-\\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}$$
")
          }
          cat("
The mean length by ",
              agesex,
              " category ", verb2, " estimated as follows:

$$\\bar{x}_z=\\sum_{t=1}^{L}\\frac{\\hat{N}_{tz}}{\\hat{N}_z}\\bar{x}_{tz}$$

with its variance approximated using a Taylor's series expansion (Mood et al. 1974):

$$\\hat{var}[\\bar{x}_z]\\approx\\sum_{t=1}^{L}\\frac{\\hat{N}_{tz}^2}{\\hat{N}_z^2}\\hat{var}[\\bar{x}_{tz}]+\\sum_{t=1}^{L}\\frac{\\left(\\bar{x}_{tz}\\hat{N}_z-\\left(\\sum_{u=1}^{L}\\bar{x}_{uz}\\hat{N}_{uz}\\right)\\right)^2}{\\hat{N}_z^4}\\hat{var}[\\hat{N}_{tz}]$$")

        }
      } else { # if there are no proportions
        cat("The mean length for each sampling stratum ", verb1, " estimated as follows, where $x_{ti}$ is the length of the *i*th ",
            species,
            " sampled within sampling stratum $t$, and $n_t$ is the number of ",
            species,
            " in stratum *t* sampled for length:

$$\\bar{x}_{t}=\\frac{\\sum_{i=1}^{n_{t}}x_{ti}}{n_{t}}$$

The sampling variance of $\\bar{x}_{t}$ ", verb1, " estimated as")
        if(doFPCmean) {
          cat(" the following, in which $N_t$ represents the abundance associated with sampling stratum $t$:

$$\\hat{var}[\\bar{x}_{t}]=\\frac{\\sum_{i=1}^{n_{t}}(x_{ti}-\\bar{x}_{t})^2}{n_{t}(n_{t}-1)}\\left(\\frac{N_{t}-n_{t}}{N_{t}-1}\\right)$$
")
        } else {
          cat("

$$\\hat{var}[\\bar{x}_{t}]=\\frac{\\sum_{i=1}^{n_{t}}(x_{ti}-\\bar{x}_{t})^2}{n_{t}(n_{t}-1)}$$
")
        }
        cat("
Stratified estimates of mean length ", verb3, " calculated as follows, in which ")
        if(!doFPCmean) cat("$N_t$ represents the abundance associated with sampling stratum $t$, ")
        cat("$N$ represents the total abundance, and *L* represents the number of sampling strata:

$$\\bar{x}=\\frac{1}{N}\\sum_{t=1}^L N_t\\bar{x}_t$$

and

$$\\hat{var}[\\bar{x}]=\\sum_{t=1}^L\\left(\\frac{N_t}{N}\\right)^2 \\hat{var}[\\bar{x}_t]$$")
      }
    }

    if(abundance == "estimated") {
      if(doProp) {
        cat("The proportion of each ", agesex, " category *z* ", verb1, " estimated for each sampling stratum *t* as follows:

$$\\hat{p}_{tz}=\\frac{n_{tz}}{n_t}$$

in which $n_{tz}$ equals the number of ", species, " sampled during sampling stratum $t$ classified as ", agesex, " category $z$, and $n_t$ equals the number of ", species, " sampled for ", agesex, " determination within sampling stratum $t$.

The sampling variance of $\\hat{p}_{tz}$ ", verb1, " estimated as the following (Cochran 1977)")
        if(doFPCprop) {
          cat(" in which $\\hat{N}_t$ is the estimated abundance of ", species, " in sampling stratum $t$:

$$\\hat{var}[\\hat{p}_{tz}]=\\frac{\\hat{p}_{tz}(1-\\hat{p}_{tz})}{n_t-1}\\left(\\frac{\\hat{N}_t-n_t}{\\hat{N}_t-1}\\right)$$
")
        } else {
          cat(":

$$\\hat{var}[\\hat{p}_{tz}]=\\frac{\\hat{p}_{tz}(1-\\hat{p}_{tz})}{n_t-1}$$")
        }
        doCasella <- TRUE
        cat("

The total abundance by ", agesex, " category in each sampling stratum ", verb1, " estimated as follows")
        if(!doFPCprop) cat(", in which $\\hat{N}_t$ represents the estimated abundance of ",species," in sampling stratum *t*")
        cat(":

$$\\hat{N}_{tz}=\\hat{N}_t\\hat{p}_{tz}$$

with variance estimated as (Goodman 1960):

$$\\hat{var}[\\hat{N}_{tz}]=\\hat{N}_t^2\\hat{var}[\\hat{p}_{tz}] + \\hat{p}_{tz}^2\\hat{var}[\\hat{N}_t]-\\hat{var}[\\hat{p}_{tz}]\\hat{var}[\\hat{p}_{tz}]$$

The total abundance by ", agesex, " category $z$ and its variance ", verb4, " estimated by summation as follows:

$$\\hat{N}_z=\\sum_{t=1}^{L}\\hat{N}_{tz}$$

and

$$\\hat{var}[\\hat{N}_{z}]=\\sum_{t=1}^{L}\\hat{var}[\\hat{N}_{tz}]$$

where $L$ equals the number of sampling strata.

Finally, the overall proportion by ", agesex, " category and its variance ", verb3, " estimated as follows:

$$\\hat{p}_z=\\frac{\\hat{N}_z}{\\sum_{t=1}^{L}\\hat{N}_t}$$

with variance approximated by the delta method (Casella & Berger 2002) as:

$$\\hat{var}[\\hat{p}_z] \\approx \\left(\\frac{\\hat{N}_z}{\\sum_{t=1}^{L}\\hat{N}_t}\\right)^2\\left(\\frac{\\hat{var}[\\hat{N}_z]}{\\hat{N}_z^2} + \\frac{\\sum_{t=1}^{L}\\hat{var}[\\hat{N}_t]}{(\\sum_{t=1}^{L}\\hat{N}_t)^2} - 2\\frac{\\hat{cov}[\\hat{N}_z,\\sum_{t=1}^{L}\\hat{N}_t]}{\\hat{N}_z\\sum_{t=1}^{L}\\hat{N}_t}\\right)$$

in which

$$\\hat{cov}[\\hat{N}_z,\\sum_{t=1}^{L}\\hat{N}_t]=\\sum_{t=1}^{L}\\hat{p}_{tz}\\hat{var}[\\hat{N_t}]$$

")
        if(doLength) {
          doMood <- TRUE
          cat("
The mean length by ", agesex, " for each sampling stratum ", verb1, " estimated as follows:

$$\\bar{x}_{tz}=\\frac{\\sum_{i=1}^{n_{tz}}x_{tzi}}{n_{tz}}$$

where $x_{tzi}$ is the length of the *i*th ", species, " sampled of ", agesex, " $z$ during sampling stratum $t$.

The sampling variance of $\\bar{x}_{tz}$ ", verb1, " estimated as

")
          if(doFPCmean) {
            cat("$$\\hat{var}[\\bar{x}_{tz}]=\\frac{\\sum_{i=1}^{n_{tz}}(x_{tzi}-\\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}\\left(\\frac{\\hat{N}_{tz}-n_{tz}}{\\hat{N}_{tz}-1}\\right)$$

")
          } else {
            cat("$$\\hat{var}[\\bar{x}_{tz}]=\\frac{\\sum_{i=1}^{n_{tz}}(x_{tzi}-\\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}$$

")
          }
          cat("
The mean length by ", agesex, " category ", verb2, " estimated as follows:

$$\\bar{x}_z=\\sum_{t=1}^{L}\\frac{\\hat{N}_{tz}}{\\hat{N}_z}\\bar{x}_{tz}$$

with its variance approximated using a Taylor's series expansion (Mood et al. 1974):

$$\\hat{var}[\\bar{x}_z]\\approx\\sum_{t=1}^{L}\\frac{\\hat{N}_{tz}^2}{\\hat{N}_z^2}\\hat{var}[\\bar{x}_{tz}]+\\sum_{t=1}^{L}\\frac{\\left(\\bar{x}_{tz}\\hat{N}_z-\\left(\\sum_{u=1}^{L}\\bar{x}_{uz}\\hat{N}_{uz}\\right)\\right)^2}{\\hat{N}_z^4}\\hat{var}[\\hat{N}_{tz}]$$
")
        }
      } else {  # no proportions
        doCasella <- TRUE
        doGoodman <- TRUE
        cat("The mean length for each sampling stratum ", verb1, " estimated as follows, where $x_{ti}$ is the length of the *i*th ", species, " sampled within sampling stratum $t$, and $n_t$ is the number of ", species, " in stratum *t* sampled for length:

$$\\bar{x}_{t}=\\frac{\\sum_{i=1}^{n_{t}}x_{ti}}{n_{t}}$$

The sampling variance of $\\bar{x}_{t}$ ", verb1, " estimated as")
        if(doFPCmean) {
          cat("the following, in which $\\hat{N}_t$ represents the estimated abundance associated with stratum *t*:

$$\\hat{var}[\\bar{x}_{t}]=\\frac{\\sum_{i=1}^{n_{t}}(x_{ti}-\\bar{x}_{t})^2}{n_{t}(n_{t}-1)}\\left(\\frac{\\hat{N}_{t}-n_{t}}{\\hat{N}_{t}-1}\\right)$$

")
        } else {
          cat("

$$\\hat{var}[\\bar{x}_{t}]=\\frac{\\sum_{i=1}^{n_{t}}(x_{ti}-\\bar{x}_{t})^2}{n_{t}(n_{t}-1)}$$

")
        }
        cat("Stratified estimates of mean length ", verb3, " calculated as follows, in which ")
        if(doFPCmean) {
          cat("$\\bar{x}_t$ represents the mean length associated with stratum *t*")
        } else {
          cat("$\\hat{N}_t$ and $\\bar{x}_t$ represent the estimated abundance and mean length associated with stratum *t*, respectively")
        }
        cat(":

$$\\bar{x}=\\frac{\\sum_{t=1}^L \\hat{N}_t\\bar{x}_t}{\\sum_{t=1}^L \\hat{N}_t}$$

and

$$\\hat{var}[\\bar{x}] \\approx \\left(\\frac{\\sum_{t=1}^L \\hat{N}_t\\bar{x}_t}{\\sum_{t=1}^L \\hat{N}_t}\\right)^2\\left(\\frac{\\hat{var}[\\sum_{t=1}^L \\hat{N}_t\\bar{x}_t]}{\\left(\\sum_{t=1}^L \\hat{N}_t\\bar{x}_t\\right)^2}+\\frac{\\sum_{t=1}^L \\hat{var}[\\hat{N}_t]}{\\left(\\sum_{t=1}^L \\hat{N}_t\\right)^2}-2\\frac{\\hat{cov}[\\sum_{t=1}^L \\hat{N}_t,\\sum_{t=1}^L N_t\\bar{x}_t]}{\\left(\\sum_{t=1}^L \\hat{N}_t\\right)\\left(\\sum_{t=1}^L \\hat{N}_t\\bar{x}_t\\right)}\\right)$$

in which

$$\\hat{cov}[\\sum_{t=1}^L \\hat{N}_t,\\sum_{t=1}^L \\hat{N}_t\\bar{x}_t]=\\sum_{t=1}^L \\bar{x}_t\\hat{var}[\\hat{N}_t]$$

and

$$\\hat{var}[\\sum_{t=1}^L \\hat{N}_t\\bar{x}_t]=\\sum_{t=1}^L\\hat{N}_t^2\\hat{var}[\\bar{x}_t] + \\bar{x}_t^2\\hat{var}[\\hat{N}_t]-\\hat{var}[\\hat{N}_t]\\hat{var}[\\bar{x}_t]$$

by means of the delta method (Casella & Berger 2002) and Goodman (1960), respectively.")
      }

    }

    if(abundance == "unknown") {
      if(doProp) {
        cat("The proportion of each ", agesex, " category *z* ", verb1, " estimated for each sampling stratum *t* as follows:

$$\\hat{p}_{tz}=\\frac{n_{tz}}{n_t}$$

in which $n_{tz}$ equals the number of ", species," sampled during sampling stratum $t$ classified as ", agesex, " category $z$, and $n_t$ equals the number of ", species," sampled for ", agesex, " determination within sampling stratum $t$.

The sampling variance of $\\hat{p}_{tz}$ ", verb1, " estimated as the following (Cochran 1977):

$$\\hat{var}[\\hat{p}_{tz}]=\\frac{\\hat{p}_{tz}(1-\\hat{p}_{tz})}{n_t-1}$$

The overall proportion by ", agesex, " category and its variance ", verb3, " estimated as follows, in which $w_t$ represents the sampling weight associated with stratum *t* and *L* equals the number of strata.  It is worth noting that weights $w_t$ are treated as constant (i.e. known without error), therefore all variance estimates must be interpreted as minima without further assumptions.

$$\\hat{p}_z=\\frac{\\sum_{t=1}^Lw_t\\hat{p}_{tz}}{\\sum_{t=1}^Lw_t}$$

and

$$\\hat{var}[\\hat{p}_z]=\\frac{\\sum_{t=1}^Lw_t^2\\hat{var}[\\hat{p}_{tz}]}{\\left(\\sum_{t=1}^Lw_t\\right)^2}$$

")
        if(doLength) {
          doMood <- TRUE
          cat("The mean length by ", agesex, " for each sampling stratum ", verb1, " estimated as follows:

$$\\bar{x}_{tz}=\\frac{\\sum_{i=1}^{n_{tz}}x_{tzi}}{n_{tz}}$$

where $x_{tzi}$ is the length of the *i*th ", species," sampled of ", agesex, " $z$ during sampling stratum $t$.

The sampling variance of $\\bar{x}_{tz}$ ", verb1, " estimated as:

$$\\hat{var}[\\bar{x}_{tz}]=\\frac{\\sum_{i=1}^{n_{tz}}(x_{tzi}-\\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}$$

The mean length by ", agesex, " category ", verb2, " estimated as follows:

$$\\bar{x}_z=\\frac{\\sum_{t=1}^{L}w_t\\hat{p}_{tz}\\bar{x}_{tz}}{\\sum_{t=1}^{L}w_t\\hat{p}_{tz}}$$

with its variance approximated using a Taylor's series expansion (Mood et al. 1974):

$$\\hat{var}[\\bar{x}_z]\\approx\\frac{\\sum_{t=1}^{L}w_t\\hat{p}_{tz}^2\\hat{var}[\\bar{x}_{tz}]}{\\left(\\sum_{t=1}^Lw_t\\hat{p}_{tz}\\right)^2}+\\frac{\\sum_{t=1}^{L}\\left(\\bar{x}_{tz}\\sum_{u=1}^Lw_u\\hat{p}_{uz}-\\left(\\sum_{u=1}^{L}\\bar{x}_{uz}w_u\\hat{p}_{uz}\\right)\\right)^2w_t^2\\hat{var}[\\hat{p}_{tz}]}{\\left(\\sum_{t=1}^Lw_t\\hat{p}_{tz}\\right)^4}$$")
        }
      } else {  # no proportions
        cat("The mean length for each sampling stratum ", verb1, " estimated as follows, where $x_{ti}$ is the length of the *i*th ", species," sampled within sampling stratum $t$, and $n_t$ is the number of ", species," in stratum *t* sampled for length:

$$\\bar{x}_{t}=\\frac{\\sum_{i=1}^{n_{t}}x_{ti}}{n_{t}}$$

The sampling variance of $\\bar{x}_{t}$ ", verb1, " estimated as:

$$\\hat{var}[\\bar{x}_{t}]=\\frac{\\sum_{i=1}^{n_{t}}(x_{ti}-\\bar{x}_{t})^2}{n_{t}(n_{t}-1)}$$

Stratified estimates of mean length ", verb3, " calculated as follows, in which $w_t$ and $\\bar{x}_t$ represent the sampling weight and average length associated with sampling stratum $t$, respectively.  It is worth noting that weights $w_t$ are treated as constant (i.e. known without error), therefore all variance estimates must be interpreted as minima without further assumptions.

$$\\bar{x}=\\frac{\\sum_{t=1}^L w_t\\bar{x}_t}{\\sum_{t=1}^L w_t}$$

and

$$\\hat{var}[\\bar{x}]=\\frac{\\sum_{t=1}^L w_t^2\\hat{var}[\\bar{x}_t]}{\\left(\\sum_{t=1}^L w_t\\right)^2}$$")
      }
    }
  } else {  # if !stratified
    if(abundance == "known") {
      if(doProp) {
        cat("Proportions of each ", agesex, " category $z$ ", verb3, " estimated as follows (Cochran 1977):

$$\\hat{p}_z=\\frac{n_z}{n}$$

and

")
        if(doFPCprop) {
          cat("$$\\hat{var}[\\hat{p}_z]=\\frac{\\hat{p}(1-\\hat{p})}{n-1}\\left(\\frac{N-n}{N-1}\\right)$$

")
        } else {
          cat("$$\\hat{var}[\\hat{p}_z]=\\frac{\\hat{p}_z(1-\\hat{p}_z)}{n-1}$$

")
        }
        cat("in which ")
        if(doFPCprop) cat("*N* denotes the total abundance, ")
        cat("$n_z$ denotes the number of ", species, " sampled in ", agesex, " category $z$, and $n$ denotes the total number of ", species, " sampled.

Total abundance for ", agesex, " category $z$ ", verb1, " estimated as")
        if(!doFPCprop) cat(" the following, in which *N* denotes the total abundance:")
        cat("

$$\\hat{N}_z=N\\hat{p}_z$$

and

$$\\hat{var}[\\hat{N}_z]=N^2\\hat{var}[\\hat{p}_z]$$

")
        if(doLength) {
          cat("The mean length associated with ", agesex, " category $z$ ", verb1, " estimated as the following, in which $x_{zi}$ represents the length of the *i*th ", species, " in ", agesex, " category *z* and $n_z$ represents the number of ", species, " in ", agesex, " category *z* with an associated length measurement:

$$\\bar{x}_z=\\frac{\\sum_{i=1}^{n_z}x_{zi}}{n_z}$$

and

")
          if(doFPCmean) {
            cat("$$\\hat{var}[\\bar{x}_z]=\\frac{\\sum_{i=1}^{n_z}(x_{zi}-\\bar{x_z})^2}{n_z(n_z-1)}\\left(\\frac{\\hat{N}_z-n_z}{\\hat{N}_z-1}\\right)$$

")
          } else {
            cat("$$\\hat{var}[\\bar{x}_z]=\\frac{\\sum_{i=1}^{n_z}(x_{zi}-\\bar{x_z})^2}{n_z(n_z-1)}$$")
          }
        }
      } else {  # no proportions
        cat("The mean length of all ", species, " ", verb1, " estimated as the following, in which ")
        if(doFPCmean) cat("*N* represents the total abundance, ")
        cat("$x_{i}$ represents the length of the *i*th ", species, ", and $n$ represents the number of ", species, " with an associated length measurement:")
        cat("

$$\\bar{x}=\\frac{\\sum_{i=1}^{n}x_{i}}{n}$$

and

")
        if(doFPCmean) {
          cat("$$\\hat{var}[\\bar{x}]=\\frac{\\sum_{i=1}^{n}(x_{i}-\\bar{x})^2}{n(n-1)}\\left(\\frac{N-n}{N-1}\\right)$$

")
        } else {
          cat("$$\\hat{var}[\\bar{x}]=\\frac{\\sum_{i=1}^{n}(x_{i}-\\bar{x})^2}{n(n-1)}$$")
        }
      }
    }

    if(abundance == "estimated") {
      if(doProp) {
        cat("Proportions of each ", agesex, " category $z$ ", verb3, " estimated as follows (Cochran 1977):

$$\\hat{p}_z=\\frac{n_z}{n}$$

and

")
        if(doFPCprop) {
          cat("$$\\hat{var}[\\hat{p}]=\\frac{\\hat{p}(1-\\hat{p})}{n-1}\\left(\\frac{\\hat{N}-n}{\\hat{N}-1}\\right)$$

")
        } else {
          cat("$$\\hat{var}[\\hat{p}]=\\frac{\\hat{p}(1-\\hat{p})}{n-1}$$

")
        }
        cat("in which ")
        if(doFPCprop) cat("$\\hat{N}$ denotes the estimated abundance, ")
        cat("$n_z$ denotes the number of ", species, " sampled in ", agesex, " category $z$, and $n$ denotes the total number of ", species, " sampled.

Total abundance for ", agesex, " category $z$ ", verb1, " estimated as follows (Goodman 1960)")
        if(!doFPCprop) cat(" in which $\\hat{N}$ denotes the estimated abundance:")
        cat("

$$\\hat{N}_z=\\hat{N}\\hat{p}_z$$

and

$$\\hat{var}[\\hat{N}_z]=\\hat{N}^2\\hat{var}[\\hat{p}_z] + \\hat{p}_z^2\\hat{var}[\\hat{N}]-\\hat{var}[\\hat{N}]\\hat{var}[\\hat{p}_z]$$

")
        if(doLength) {
          cat("The mean length associated with ", agesex, " category $z$ ", verb1, " estimated as the following, in which $x_{zi}$ represents the length of ", species, " *i* within ", agesex, " category category *z*, and $n_z$ denotes the number of ", species, " within ", agesex, " category *z* with an associated length measurement:

$$\\bar{x}_z=\\frac{\\sum_{i=1}^{n_z}x_{zi}}{n_z}$$

and

")
          if(doFPCmean) {
            cat("$$\\hat{var}[\\bar{x}_z]=\\frac{\\sum_{i=1}^{n_z}(x_{zi}-\\bar{x_z})^2}{n_z(n_z-1)}\\left(\\frac{\\hat{N}_z-n_z}{\\hat{N}_z-1}\\right)$$

")
          } else {
            cat("$$\\hat{var}[\\bar{x}_z]=\\frac{\\sum_{i=1}^{n_z}(x_{zi}-\\bar{x_z})^2}{n_z(n_z-1)}$$")
          }
        }
      } else {  # no proportions
        cat("The mean length of all ", species, " ", verb1, " estimated as the following, in which ")
        if(doFPCmean) cat("$\\hat{N}$ represents the estimated total abundance, ")
        cat("$x_{i}$ represents the length of the *i*th ", species, ", and $n$ represents the number of ", species, " with an associated length measurement:

$$\\bar{x}=\\frac{\\sum_{i=1}^{n}x_{i}}{n}$$

and

")
        if(doFPCmean) {
          cat("$$\\hat{var}[\\bar{x}]=\\frac{\\sum_{i=1}^{n}(x_{i}-\\bar{x})^2}{n(n-1)}\\left(\\frac{\\hat{N}-n}{\\hat{N}-1}\\right)$$

")
        } else {
          cat("$$\\hat{var}[\\bar{x}]=\\frac{\\sum_{i=1}^{n}(x_{i}-\\bar{x})^2}{n(n-1)}$$")
        }
      }
    }

    if(abundance == "unknown") {
      if(doProp) {
        cat("Proportions of each ", agesex, " category $z$ ", verb3, " estimated as follows (Cochran 1977):

$$\\hat{p}_z=\\frac{n_z}{n}$$

and

$$\\hat{var}[\\hat{p}_z]=\\frac{\\hat{p}_z(1-\\hat{p}_z)}{n-1}$$

in which $n_z$ denotes the number of ", species, " sampled in ", agesex, " category $z$, and $n$ denotes the total number of ", species, " sampled.

")
        if(doLength) {
          cat("The mean length associated with ", agesex, " category $z$ ", verb1, " estimated as the following, in which $x_{zi}$ denotes the length measurement associated with the *i*th ", species, " in ", agesex, " category *z*:

$$\\bar{x}_z=\\frac{\\sum_{i=1}^{n_z}x_{zi}}{n_z}$$

and

$$\\hat{var}[\\bar{x}_z]=\\frac{\\sum_{i=1}^{n_z}(x_{zi}-\\bar{x_z})^2}{n_z(n_z-1)}$$
")
        }
      } else {
        cat("The mean length ", verb1, " estimated as the following, in which $x_{i}$ denotes the length measurement associated with the *i*th ", species, " and $n$ denotes the number of ", species, " with an associated length measurement:

$$\\bar{x}=\\frac{\\sum_{i=1}^{n}x_{i}}{n}$$

and

$$\\hat{var}[\\bar{x}]=\\frac{\\sum_{i=1}^{n}(x_{i}-\\bar{x})^2}{n(n-1)}$$")
      }
    }
  }

  ### Printing references

  #   cat("
  #
  # All calculations were performed in R^[",
  #       capture.output(print(citation(), style="text")),
  #       "] using the dsftools package.^[",
  #       capture.output(print(citation("dsftools"), style="text")),
  #       "]", sep="")
  cat("

All calculations were performed in R^[")
  cat(capture.output(print(citation(), style="text")))
  cat("] using the dsftools package.^[")
  cat(capture.output(print(citation("dsftools"), style="text")))
  cat("]")

  cat("

## References

")
  if(doCasella) {
    cat("Casella, George and Roger L. Berger. 2002. *Statistical Inference*. Australia ; Pacific Grove, CA, Thomson Learning

")
  }
  if(doCochran) {
    cat("Cochran, W. G.  1977.  *Sampling techniques*. 3rd edition.  John Wiley and Sons, New York.

")
  }
  if(doGoodman) {
    cat("Goodman, L.A., 1960. On the exact variance of products. *Journal of the American statistical association, 55*(292), pp.708-713.

")
  }
  if(doMood) {
    cat("Mood, A. M., F. A. Graybill, and D. C. Boes.  1974.  *Introduction to the theory of statistics*. 3rd edition.  McGraw-Hill Book Co., New York.

")
  }
  # print(citation(), style="text")
  # cat("\n")
  # print(citation("dsftools"), style="text")
  # cat("\n")
  # print(citation("knitr"), style="text")
}
