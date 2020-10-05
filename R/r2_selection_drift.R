#' Calculates partinioned R-squared for selection-drift model
#'
#' @param model drift model
#' @param null.model the correspondent null model
#'
#' @return
#' @export
#'
r2_selection_drift <- function(model, null_model){
  ## Calculates null model
  if (missing(null_model)) {
    m0 <- best_null(model)
  } else {
    m0 <- null_model
  }

  ## Variance for fixed effects
  VarF <- var(as.vector(fixef(model) %*% t(model@pp$X)))

  ## Variance for slope/intercept random effects
  X <- model.matrix(model)
  n <- nrow(X)
  Z <- X[, c("(Intercept)","grad")]

  # Check if model has trait or not

  if (trait) {
    sigma <- VarCorr(model)$site
    VarSite <- sum(diag(Z %*% sigma %*% t(Z)))/n

    ## Variance for each component of random effects
    VarR <-  c(VarSite,
               VarCorr(model)$`spp:site`[1],
               VarCorr(model)$`spp:region`[1],
               VarCorr(model)$`spp`[1])
    names(VarR) <- c("(1+grad|site)", "(1|spp:site)", "(1|spp:region)", "(1|spp)")

    ## Denominator for R2GLMM formula works for Poisson distribution only
    deno <- (VarF + sum(VarR) +
               log(1 + 1/exp(as.numeric(fixef(m0)))))

    ## R2GLMM(m) - marginal R2GLMM
    r2f <- VarF/deno

    ## R2GLMM(c) - conditional R2GLMM for full model
    r2t <- (VarF + sum(VarR))/deno

    # R2 random effects only
    r2rand <- r2t - r2f

    ## R2 Residuals
    r2res <- 1 - r2t

    ## Partitioning R2 GLMM for each random effect
    r2rand.part <- VarR/deno
    r2.tab <- data.frame(component = c("conditional", "fixed", "random",
                                       names(VarR)),
                         R2 = c(r2t,r2f,r2rand, r2rand.part),
                         type = c("all", "niche", "all.random",
                                  "niche", "neutral", "neutral", "niche"))
    r2.partition <- aggregate(r2.tab$R2, list(type = r2.tab$type), sum)
    names(r2.partition)[2] <- "R2.partition"
    R2 <- r2.tab$R2
    names(R2) <- r2.tab$component
  } else {
    sigma.site <- VarCorr(model)$site
    VarSite <- sum(diag(Z %*% sigma.site %*% t(Z)))/n
    ## Variance for slope/intercept random effects
    ## spp
    sigma.spp <- VarCorr(model)$spp
    VarSpp <- sum(diag(Z %*% sigma.spp %*% t(Z)))/n
    ## Variance for each component of random effects
    VarR <-  c(VarSite, VarSpp,
               VarCorr(model)$`spp:site`[1],
               VarCorr(model)$`spp:region`[1])
    names(VarR) <- c("(1+grad|site)", "(1+grad|spp)", "(1|spp:site)", "(1|spp:region)")
    # Denominator for R2GLMM formula works for Poisson distribution only
    deno <- (VarF + sum(VarR) + #pi^2/3
               log(1 + 1/exp(as.numeric(fixef(m0)))))
    # R2GLMM(m) - marginal R2GLMM
    r2f <- VarF/deno
    # R2GLMM(c) - conditional R2GLMM for full model
    r2t <- (VarF + sum(VarR))/deno
    # R2 random effects only
    r2rand <- r2t - r2f
    ## R2 Residuals
    r2res <- 1 - r2t
    ## Partitioning R2 GLMM for each random effect
    r2rand.part <- VarR/deno
    r2.tab <- data.frame(component = c("conditional", "fixed", "random",
                                     names(VarR)),
                         R2 = c(r2t,r2f,r2rand, r2rand.part),
                         type = c("all", "niche", "all.random",
                                  "niche",  "niche", "neutral", "neutral"))
    r2.partition <- aggregate(r2.tab$R2, list(type = r2.tab$type), sum)
    names(r2.partition)[2] <- "R2.partition"
    R2 <- r2.tab$R2
    names(R2) <- r2.tab$component
  }

  return(list(full.table = r2.tab, part.table = r2.partition, R2 = R2))

}