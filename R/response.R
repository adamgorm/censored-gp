#' Describe Observation-Specific Censoring
#'
#' Use this response on the left of a [cgp()] formula. Every row may have a
#' different censoring type and limits. Bounds describe noisy measurements,
#' rather than hard constraints on the underlying smooth curve.
#'
#' @param y Numeric values: the exact measurement for `"none"`, the upper
#'   detection limit for `"left"`, or the lower limit for `"right"`.
#' @param type Character vector of `"none"`, `"left"`, or `"right"`.
#'   A single value is recycled to all rows. A logical vector is
#'   shorthand for left censoring (`TRUE`) versus exact observation (`FALSE`).
#' @return A two-column matrix of lower and upper bounds, of class `cgp_response`.
#' @examples
#' # Exact, left-censored, and right-censored observations in one response.
#' cens(c(1, 2, 3), c("none", "left", "right"))
#' # Observation-specific detection limits are supplied directly in y.
#' cens(c(1.4, 1.7, 2.1), c(FALSE, TRUE, TRUE))
#' @export
cens <- function(y, type = "none") {
  if (!is.numeric(y) || !length(y) || any(!is.finite(y))) {
    stop("`y` must contain finite numeric measurements or censoring limits.", call. = FALSE)
  }
  n <- length(y)
  if (is.logical(type)) type <- ifelse(type, "left", "none")
  if (is.factor(type)) type <- as.character(type)
  if (length(type) == 1L) type <- rep(type, n)
  allowed <- c("none", "left", "right")
  if (!is.character(type) || length(type) != n || anyNA(type) || any(!type %in% allowed)) {
    stop("`type` must give none, left, or right for every observation.", call. = FALSE)
  }
  lower_bound <- upper_bound <- as.numeric(y)
  lower_bound[type == "left"] <- -Inf
  upper_bound[type == "right"] <- Inf
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
  exact <- is.finite(lower) & is.finite(upper) & lower == upper
  left <- is.infinite(lower) & lower < 0 & is.finite(upper)
  right <- is.finite(lower) & is.infinite(upper) & upper > 0
  # Reject manually constructed responses outside the public censoring types.
  if (anyNA(c(lower, upper)) || !all(exact | left | right)) {
    stop("Responses must be exact, left censored, or right censored.", call. = FALSE)
  }
  type <- ifelse(exact, "none", ifelse(left, "left", "right"))
  data.frame(lower = lower, upper = upper, type = type)
}
