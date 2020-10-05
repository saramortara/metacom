#' Function to fit the Drift model
#'
#' Each observation must me relative to a particular species at a particular sampling site.
#'
#' @inheritParams fit_selection_drift
#'
#' @return A linear mixed model of class glmer
#'
#' @export
#'
#' @importFrom lme4 glmer
#'
fit_null <- function(ab, site, region, spp, ...){
  m_null <- glmer(ab ~ 1 + (1|region) + (1|spp) + (1|site), ...)
  return(m_null)
}