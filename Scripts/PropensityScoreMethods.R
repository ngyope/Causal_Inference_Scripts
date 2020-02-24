library(tidyverse)
library(MatchIt)
library(WeightIt)
library(cobalt)

# get the lalonde data
data(lalonde)
summary(lalonde)

# define formulas
# formula of outcome model
frml_out <- as.formula(re78 ~ treat + age + educ + race + married + nodegree + re74 + I(re74^2))
# formula of treatment model
frml_treat <- as.formula(treat ~ age + educ + race + married + nodegree + re74 + I(re74^2))

# matching methods
# nearest neighbors
m_out_nn_1 <- matchit(frml_treat, lalonde, method="nearest", replace=TRUE, ratio=1, discard="control")
m_out_nn_4 <- matchit(frml_treat, lalonde, method="nearest", replace=TRUE, ratio=4, discard="control")
# subclassification
m_out_subclass <- matchit(frml_treat, lalonde, method="subclass",subclass=5, discard="control")
# caliper
m_out_caliper <- matchit(frml_treat, lalonde, method="nearest", replace=TRUE, caliper=.025, ratio=1, discard="control")

# weighting methods
w_out_ps <- weightit(frml_treat, lalonde, method="ps", estimand="ATT")
w_out_gbm <- weightit(frml_treat, lalonde, method="gbm", stop.method="es.max", estimand="ATT")
w_out_cbps <- weightit(frml_treat, lalonde, method="ps", estimand="ATT")

# check the covariate balance
love.plot(frml_treat, data=lalonde,
          weights=data.frame(Nearest_Neighbor_1=get.w(m_out_nn_1),
                             Nearest_Neighbor_4=get.w(m_out_nn_4),
                             Subclassification=get.w(m_out_subclass),
                             Caliper=get.w(m_out_caliper),
                             PS_Weighting=get.w(w_out_ps)),
          method=c("matching", "matching", "matching", "matching", "weighting"),
          binary="std", s.d.denom="treated", grid=FALSE, threshold=.1, limits=c(0,1))

# check the covariate balance per covariate
bal.plot(m_out_nn_1, var.name="distance", which="both", type="histogram", mirror=TRUE)

# get the ATT
matched_data <- list(match.data(m_out_nn_1),
                     match.data(m_out_nn_4),
                     match.data(m_out_subclass),
                     match.data(m_out_caliper))
model_label <- c("nearest neighbor ratio 1", "nearest neighbor ratio 4", "subclassification", "caliper")

for(i in 1:length(matched_data)){
  print(model_label[[i]])
  print(lm(re78 ~ treat, data=matched_data[[i]]) %>% broom::tidy())
}
