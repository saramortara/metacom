#' Function to fit the Selection model
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
selection <- function(ab, trait = NULL, grad, site, spp, ...){
  if (is.null) {
    m_selection <- glmer(ab ~ grad + I(grad^2) +
                           (1+grad|spp) + (1+grad|site), ...)
  } else {
    m_selection <- glmer(ab ~ trait + grad + I(grad^2) +
                           trait:grad + trait:I(grad^2) +
                           (1|spp) + (1+grad|site), ...)
  }
  return(m_selection)
}