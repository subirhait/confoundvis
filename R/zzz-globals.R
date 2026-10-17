# Suppress R CMD CHECK notes for column names used in ggplot2 non-standard
# evaluation via .data$x pronoun. Every name referenced as .data$name in
# aes() calls must be listed here.
utils::globalVariables(c(
  "lambda", "y", "lo", "hi", "level",
  "covariate", "impact",
  "r_yu", "r_du", "label",
  "delta", "value", "series",
  "group", "q", "p", "zone", "theta_adj"
))
