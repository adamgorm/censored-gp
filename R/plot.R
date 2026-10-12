#' Plot Posterior Functional Curves and Censored Observations
#'
#' Left and right censoring are shown by downward and upward triangles.
#' Ribbons are pointwise credible intervals for the smooth latent function.
#'
#' @param x A fitted [cgp()] model.
#' @param newdata Optional plotting locations. A grid is constructed automatically
#'   for formulas with one numeric input variable.
#' @param type `"population"` or `"latent"` (population plus group deviation).
#' @param groups Group labels for automatically generated latent-curve grids.
#'   Defaults to the first six groups; supply all labels to show every group.
#' @param n Number of grid locations per curve.
#' @param level Pointwise credible interval coverage.
#' @param ndraws Number of posterior draws; defaults to at most 400.
#' @param seed Seed for prediction simulation.
#' @param ... Reserved for S3 compatibility.
#' @return A `ggplot2` plot, which may be customized with standard ggplot layers.
#' @importFrom rlang .data
#' @export
plot.cgp_fit <- function(x, newdata = NULL, type = c("population", "latent"),
                         groups = NULL, n = 100L, level = 0.95,
                         ndraws = NULL, seed = 1L, ...) {
  type <- match.arg(type)
  prepared <- x$prepared
  if (ncol(prepared$coordinates) != 1L) {
    stop("Automatic curve plots require one-dimensional GP inputs.", call. = FALSE)
  }
  if (is.null(newdata)) {
    variables <- all.vars(prepared$formula_parts$coordinates)
    if (length(variables) != 1L || !variables %in% names(x$data)) {
      stop("Supply `newdata` when the coordinate expression uses multiple variables.", call. = FALSE)
    }
    n <- check_count(n, "n", 2L)
    variable <- variables[[1L]]
    newdata <- x$data[rep(1L, n), , drop = FALSE]
    newdata[[variable]] <- seq(min(x$data[[variable]]), max(x$data[[variable]]), length.out = n)
    if (type == "latent" && !is.null(prepared$formula_parts$group)) {
      group_variables <- all.vars(prepared$formula_parts$group)
      if (length(group_variables) != 1L) stop("Supply newdata for a compound grouping expression.", call. = FALSE)
      if (is.null(groups)) groups <- utils::head(prepared$group_levels, 6L)
      newdata <- do.call(rbind, lapply(as.character(groups), function(group) {
        result <- newdata
        # Convert factors to character before adding independently modeled groups.
        result[[group_variables]] <- rep(group, n)
        result
      }))
    }
  }
  if (is.null(ndraws)) ndraws <- min(400L, posterior::ndraws(x$draws))
  predictions <- predict.cgp_fit(x, newdata, type = type, level = level, ndraws = ndraws, seed = seed)
  design <- prediction_design(x, newdata, if (type == "latent") "response" else "population",
    joint = FALSE
  )
  predictions$coordinate <- design$coordinates[, 1L]
  predictions$group <- if (type == "latent") design$group else ".population"
  predictions <- predictions[order(predictions$group, predictions$coordinate), ]
  observations <- prepared$bounds
  observations$value <- ifelse(observations$type == "left", observations$upper,
    observations$lower
  )
  observations$coordinate <- prepared$coordinates[, 1L]
  observations$group <- if (type == "latent") prepared$group else ".population"
  observations <- observations[observations$group %in% predictions$group, , drop = FALSE]
  plot <- ggplot2::ggplot(predictions, ggplot2::aes(x = .data$coordinate, y = .data$mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "#2878A0", alpha = 0.2
    ) +
    ggplot2::geom_line(color = "#176487", linewidth = 0.7) +
    ggplot2::geom_point(
      data = observations,
      ggplot2::aes(y = .data$value, shape = .data$type), alpha = 0.65
    ) +
    ggplot2::scale_shape_manual(
      values = c(none = 16, left = 25, right = 24),
      name = "Observation", labels = c(
        none = "Exact", left = "Left censored",
        right = "Right censored"
      )
    ) +
    ggplot2::labs(x = prepared$coordinate_label, y = "Response") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  if (type == "latent" && length(unique(predictions$group)) > 1L) {
    plot <- plot + ggplot2::facet_wrap(ggplot2::vars(.data$group))
  }
  plot
}
