#' Describe Observation-Specific Censoring
#'
#' Use this response on the left of a [cgp()] formula. Every row may have a
#' different censoring type and limits. Bounds describe noisy measurements,
#' rather than hard constraints on the underlying smooth curve.
#'
#' @param y Numeric values: the exact measurement for `"none"`, the upper
#'   detection limit for `"left"`, the lower limit for `"right"`, or the lower
#'   interval endpoint for `"interval"`.
#' @param type Character vector of `"none"`, `"left"`, `"right"`, or
#'   `"interval"`. A single value is recycled to all rows. A logical vector is
#'   shorthand for left censoring (`TRUE`) versus exact observation (`FALSE`).
#' @param upper Upper endpoints for interval-censored rows; ignored for other
#'   types. May be a scalar or a vector the same length as `y`.
#' @return A two-column matrix of lower and upper bounds, of class `cgp_response`.
#' @examples
#' # Exact, left, right, and interval censoring in one response.
#' cens <- cens(c(1, 2, 3, 4), c("none", "left", "right", "interval"),
#'   upper = c(NA, NA, NA, 5)
#' )
#' cens
#' # Observation-specific detection limits are supplied directly in y.
#' cens(c(1.4, 1.7, 2.1), c(FALSE, TRUE, TRUE))
#' @export
cens <- function(y, type = "none", upper = NULL) {
  if (!is.numeric(y) || !length(y) || any(!is.finite(y))) {
    stop("`y` must contain finite numeric measurements or censoring limits.", call. = FALSE)
  }
  n <- length(y)
  if (is.logical(type)) type <- ifelse(type, "left", "none")
  if (is.factor(type)) type <- as.character(type)
  if (length(type) == 1L) type <- rep(type, n)
  allowed <- c("none", "left", "right", "interval")
  if (!is.character(type) || length(type) != n || anyNA(type) || any(!type %in% allowed)) {
    stop("`type` must give none, left, right, or interval for every observation.", call. = FALSE)
  }
  lower_bound <- upper_bound <- as.numeric(y)
  lower_bound[type == "left"] <- -Inf
  upper_bound[type == "right"] <- Inf
  intervals <- type == "interval"
  if (any(intervals)) {
    if (length(upper) == 1L) upper <- rep(upper, n)
    if (!is.numeric(upper) || length(upper) != n || any(!is.finite(upper[intervals]))) {
      stop("Supply finite `upper` endpoints for every interval-censored row.", call. = FALSE)
    }
    if (any(upper[intervals] <= y[intervals])) {
      stop("Interval upper endpoints must exceed their lower endpoints.", call. = FALSE)
    }
    upper_bound[intervals] <- upper[intervals]
  }
  structure(cbind(lower = lower_bound, upper = upper_bound),
    class = c("cgp_response", "matrix", "array")
  )
}

#' @export
print.cgp_response <- function(x, ...) {
  print(as_censoring_data(x), row.names = FALSE, ...)
  invisible(x)
}

# Convert bounds to a readable table without relying on attributes that could
# be lost when a response matrix is subsetted by standard R tools.
as_censoring_data <- function(response) {
  lower <- response[, 1L]
  upper <- response[, 2L]
  type <- rep("interval", length(lower))
  type[lower == upper] <- "none"
  type[is.infinite(lower)] <- "left"
  type[is.infinite(upper)] <- "right"
  data.frame(lower = lower, upper = upper, type = type)
}
