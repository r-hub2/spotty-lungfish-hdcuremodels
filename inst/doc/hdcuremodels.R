## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(hdcuremodels)
library(survival)

## -----------------------------------------------------------------------------
withr::local_seed(23)
data <- generate_cure_data(n = 200, j = 10, n_true = 5, a = 1.8, rho = 0.2)
training <- data$training
testing <- data$testing
head(training)

## -----------------------------------------------------------------------------
data$parameters

## -----------------------------------------------------------------------------
names(training)[grep("^X", names(training))][data$parameters$nonzero_b]

## -----------------------------------------------------------------------------
names(training)[grep("^X", names(training))][data$parameters$nonzero_beta]

## -----------------------------------------------------------------------------
km_train <- survfit(Surv(cryr, relapse.death) ~ 1, data = amltrain)

## ----echo=FALSE, fig=TRUE-----------------------------------------------------
plot(km_train, mark.time = TRUE, xlab = "Time (years)", ylab = "Relapse-free survival")

## ----echo=TRUE----------------------------------------------------------------
nonzerocure_test(km_train)

## -----------------------------------------------------------------------------
cure_estimate(km_train)

## ----echo=TRUE----------------------------------------------------------------
sufficient_fu_test(km_train)

## ----args---------------------------------------------------------------------
args(curegmifs)

## ----eval=FALSE---------------------------------------------------------------
# fitgmifs <- curegmifs(Surv(cryr, relapse.death) ~ .,
#   data = amltrain,
#   x_latency = amltrain, model = "weibull"
# )

## ----args2--------------------------------------------------------------------
args(cureem)

## -----------------------------------------------------------------------------
fitem <- cureem(Surv(cryr, relapse.death) ~ .,
  data = amltrain,
  x_latency = amltrain, model = "cox",
  lambda_inc = 0.009993, lambda_lat = 0.02655
)

## ----args3--------------------------------------------------------------------
args(cv_cureem)

## -----------------------------------------------------------------------------
fit_cv <- cv_cureem(Surv(Time, Censor) ~ .,
  data = training,
  x_latency = training, fdr_control = FALSE,
  grid_tuning = FALSE, nlambda_inc = 10,
  nlambda_lat = 10, n_folds = 2, seed = 23,
  verbose = TRUE
)

## ----eval = FALSE-------------------------------------------------------------
# lambda_inc <- lambda_lat <- rep(0, 100)
# for (k in 1:100) {
#   print(k)
#   coxem_auc_k <- cv_cureem(Surv(cryr, relapse.death) ~ .,
#     data = amltrain, x_latency = amltrain,
#     model = "cox", penalty = "lasso",
#     scale = TRUE, grid_tuning = TRUE,
#     nfolds = 10, nlambda_inc = 20,
#     nlambda_lat = 20, verbose = FALSE,
#     parallel = TRUE, measure_inc = "auc"
#   )
#   lambda_inc[k] <- coxem_auc_k$selected_lambda_inc
#   coxem_c_k <- cv_cureem(Surv(cryr, relapse.death) ~ .,
#     data = amltrain,
#     x_latency = amltrain, model = "cox",
#     penalty = "lasso", scale = TRUE,
#     grid_tuning = TRUE, nfolds = 10,
#     nlambda_inc = 20, nlambda_lat = 20,
#     verbose = FALSE, parallel = TRUE,
#     measure_inc = "c"
#   )
#   lambda_lat[k] <- coxem_c_k$selected_lambda_lat
# }
# table(lambda_inc)
# table(lambda_lat)

## ----args4--------------------------------------------------------------------
args(cv_curegmifs)

## ----print, eval = FALSE------------------------------------------------------
# print(fitem)

## ----summary------------------------------------------------------------------
summary(fitem)

## ----cv-----------------------------------------------------------------------
summary(fit_cv)

## ----plot---------------------------------------------------------------------
plot(fitem)

## ----plotAIC------------------------------------------------------------------
plot(fitem, type = "cAIC")

## ----cv2----------------------------------------------------------------------
plot(fit_cv)

## ----cAIC---------------------------------------------------------------------
coef_cAIC <- coef(fitem, model_select = "cAIC")

## ----m12----------------------------------------------------------------------
coef_12 <- coef(fitem, model_select = 12)

## ----compareAIC---------------------------------------------------------------
names(coef_cAIC)
all.equal(coef_cAIC$rate, coef_12$rate)
all.equal(coef_cAIC$alpha, coef_12$alpha)
all.equal(coef_cAIC$b0, coef_12$b0)
all.equal(coef_cAIC$beta_inc, coef_12$beta_inc)
all.equal(coef_cAIC$beta_lat, coef_12$beta_lat)

## ----pred---------------------------------------------------------------------
train_predict <- predict(fitem, model_select = "cAIC")

## ----group--------------------------------------------------------------------
p_group <- ifelse(train_predict$p_uncured < 0.50, "Cured", "Susceptible")

## ----km-----------------------------------------------------------------------
km_cured <- survfit(Surv(cryr, relapse.death) ~ p_group, data = amltrain)

## ----plotkm, echo = FALSE, fig = TRUE-----------------------------------------
plot(km_cured, mark.time = TRUE, lty = c(1, 2), xlab = "Time (years)", ylab = "Relapse-free survival")
legend(c(.9, .1), legend = c("Cured", "Susceptible"), lty = c(1, 2), bty = "n")

## ----suscept------------------------------------------------------------------
km_suscept <- survfit(Surv(cryr, relapse.death) ~ train_predict$latency_risk, data = amltrain, 
                      subset = (p_group == "Susceptible"))

## ----plotsusc, echo = FALSE, fig = TRUE---------------------------------------
plot(km_suscept, mark.time = TRUE, lty = c(1, 2), xlab = "Time (years)", ylab = "Relapse-free survival")
legend(c(.9, .1), legend = c("Higher risk", "Lower risk"), lty = c(1, 2), bty = "n")

## ----testpred-----------------------------------------------------------------
test_predict <- predict(fitem, newdata = amltest, model_select = "cAIC")

## ----pgroup-------------------------------------------------------------------
test_p_group <- ifelse(test_predict$p_uncured < 0.50, "Cured", "Susceptible")

## ----kmtest-------------------------------------------------------------------
km_cured_test <- survfit(Surv(cryr, relapse.death) ~ test_p_group, data = amltest)

## ----plottest, echo = FALSE, fig = TRUE---------------------------------------
plot(km_cured_test, mark.time = TRUE, lty = c(1, 2), xlab = "Time (years)", ylab = "Relapse-free survival")
legend(c(.4, .1), legend = c("Cured", "Susceptible"), lty = c(1, 2), bty = "n")

## -----------------------------------------------------------------------------
km_suscept_test <- survfit(Surv(cryr, relapse.death) ~ test_predict$latency_risk, 
                          data = amltest, subset = (test_p_group == "Susceptible"))

## ----echo = FALSE, fig = TRUE-------------------------------------------------
plot(km_suscept_test, mark.time = TRUE, lty = c(1, 2), xlab = "Time (years)", ylab = "Relapse-free survival")
legend(c(.4, .1), legend = c("Higher risk", "Lower risk"), lty = c(1, 2), bty = "n")

## -----------------------------------------------------------------------------
auc_mcm(fitem, model_select = "cAIC")
auc_mcm(fitem, newdata = amltest, model_select = "cAIC")
concordance_mcm(fitem, model_select = "cAIC")
concordance_mcm(fitem, newdata = amltest, model_select = "cAIC")

