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
neutral <- function(ab, site, region, spp, ...){
  m_neutral <- glmer(ab ~ 1 +
                       + (1|site) +
                       (1|spp:region) + (1|spp:site), ...)
  return(m_neutral)
}