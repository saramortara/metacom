#' Function to fit the Trait-mediated Selection & Drift model
#'
#' Each observation must me relative to a particular species at a particular sampling site.
#'
#' @param ab vector of species abundance
#' @param trait vector of species traits. If `NULL` a model with a random effect to correct for the absence of traits is fitted
#' @param grad vector of environmental variable
#' @param site vector of localities
#' @param region vector of regions (localities should be nested within regions)
#' @param spp vector of species identities
#' @param ... any parameter from glmer
#'
#' @return A linear mixed model of class glmer
#'
#' @export
#'
#' @importFrom lme4 glmer
#'
fit_selection_drift <- function(ab, trait = NULL, grad, site, region, spp, ...){
  if (is.null(trait)) {
    m_full <- glmer(ab ~ grad + I(grad^2) +
                      (1|spp:region) + (1|spp:site) +
                      (1+grad|spp) + (1+grad|site), ...)
  } else {
    m_full <- glmer(ab ~ trait + grad + I(grad^2) +
                      trait:grad + trait:I(grad^2) +
                      (1|spp) +
                      (1|spp:region) + (1|spp:site) + (1+grad|site), ...)

  }
  return(m_full)
}
