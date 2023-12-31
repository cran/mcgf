#' Ireland wind data, 1961-1978
#'
#' Daily average wind speeds for 1961-1978 at 11 synoptic meteorological
#' stations in the Republic of Ireland (Haslett and raftery 1989). Wind speeds
#' are in m/s. De-trended data sets are also provided.
#'
#' @name wind
#' @format `wind`: a list containing a data.frame with 6574 rows and 12 columns,
#' and a list of locations.
#'
#' @details
#' The data were obtained from the **gstat** package, and were modified so that
#' the first column is the time stamps. Locations of the 11 stations are given
#' in `wind_loc`. `wind_train` and `wind_test` contain de-trended and
#' square-root transformed train (1961-1970) and test (1971-1978) data sets.
#' See Gneiting et al. (2006) for de-trending details. `wind_trend` contains
#' the estimated annual trend and station-wise mean from the training dataset.
#'
#' @references
#' Haslett, J. and Raftery, A. E. (1989). Space-time Modelling with Long-memory
#' Dependence: Assessing Ireland's Wind Power Resource (with Discussion).
#' Applied Statistics 38, 1-50.
#'
#' Gneiting, T., Genton, M., & Guttorp, P. (2006). Geostatistical Space-Time
#' Models, Stationarity, Separability, and Full Symmetry. In C&amp;H/CRC
#' Monographs on Statistics &amp; Applied Probability (pp. 151–175).
#' Chapman and Hall/CRC.
#' @family (simulated) datasets
"wind"
