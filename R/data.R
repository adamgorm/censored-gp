#' Load the HIV-1 RNA Example from the Paper
#'
#' Measurements below 50 copies/mL are left-censored at the assay limit.
#' `log_rna` contains the observed log10 value or this limit, rather than an
#' unobserved value below the limit. The original recorded measurements are
#' retained in `rna` for provenance.
#'
#' @return Data frame with 365 measurements from 44 subjects: `subject`, `day`,
#'   `rna`, `log_rna`, `censoring`, and `limit`.
#' @references Hoffmann et al. (2025), \doi{10.1002/sim.70277}.
#' @examples
#' hiv <- hiv_a5055()
#' table(hiv$censoring)
#' @export
hiv_a5055 <- function() {
  original <- utils::read.delim(system.file("extdata", "A5055data.txt", package = "censoredGP"))
  data.frame(
    subject = factor(original$Subject), day = original$day, rna = original$rna,
    log_rna = log10(pmax(original$rna, 50)),
    censoring = ifelse(original$rna < 50, "left", "none"), limit = log10(50)
  )
}
