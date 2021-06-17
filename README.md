# CodingModels

Fit one of three coding models --- BACE, MACE, or DS --- on long-format coder x object labels.

Example code:
```R
# if not already installed:
# devtools::install_github("matthewtyler/CodingModels")

library(CodingModels)

set.seed(123)

n_objects <- 200
n_cats <- 3
n_coders <- 5

true_labels <- data.frame(
  ii = 1:n_objects,
  truth = sample(n_cats, n_objects, replace = TRUE)
)

coder_acc <- data.frame(
  jj = 1:n_coders,
  acc = runif(n_coders, 0.6, 0.99)
)

test_labels <- expand.grid(
  ii = 1:n_objects,
  jj = 1:n_coders
)

test_labels <- merge(test_labels, true_labels,
  by = "ii", all.x = TRUE)
test_labels <- merge(test_labels, coder_acc,
  by = "jj", all.x = TRUE)

guess_labels <- sample(n_cats, nrow(test_labels), replace = TRUE)
correct <- rbinom(nrow(test_labels), 1, test_labels$acc)


test_labels$yy <- correct * test_labels$truth +
  (1-correct) * guess_labels

fit <- fit_em(test_labels, model = "BACE")

round(cbind(fit$pars$beta, coder_acc$acc), 2)
```