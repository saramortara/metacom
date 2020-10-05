## Function to calculate the null model, model with all random terms
best_null <- function(model) {
  parens <- function(x) paste0("(", x, ")")
  onlyBars <- function(form) reformulate(sapply(findbars(form),
                                                function(x)  parens(deparse(x))),
                                         response = ".")
  onlyBars(formula(model))
  best_null <- update(model, onlyBars(formula(model)), ...)
  return(best_null)
}