library(tidyverse)
library(haven)
library(MatchIt)
library(WeightIt)
library(cobalt)

# This scripts uses the famous Lalonde data. So, let's get it.
cps1_data <- read_dta("https://users.nber.org/~rdehejia/data/cps_controls.dta")
cps3_data <- read_dta("https://users.nber.org/~rdehejia/data/cps_controls3.dta")
nswdw_data <- read_dta("https://users.nber.org/~rdehejia/data/nsw_dw.dta")

# Extract the treatment group from NSW Data and treat it as a treatment group in CPS1.
cps1_nsw_data <- nswdw_data %>% filter(treat==1) %>% rbind(cps1_data)
# Extract the treatment group from NSW Data and treat it as a treatment group in CPS3.
cps3_nsw_data <- nswdw_data %>% filter(treat==1) %>% rbind(cps3_data)



# Define the formulas previously.
# Formula of outcome model
frml_out <- as.formula(re78 ~ treat + age + education + black + hispanic + married + nodegree + re74 + re75 + I(re74^2) + I(re75^2))

# Formula of treatment model
frml_treat <- as.formula(treat ~ age + education + black + hispanic + married + nodegree + re74 + re75 + I(re74^2) + I(re75^2))



# Create matching models. We're going to estimate the ATT.
# Therefore, discard the control units whose propensity score are smaller than the min of treatment units' ps.

# nearest neighbors
m_out_nn_1 <- matchit(frml_treat, cps1_nsw_data, method="nearest", ratio=1, replace=TRUE)
m_out_nn_4 <- matchit(frml_treat, cps1_nsw_data, method="nearest", ratio=4, replace=TRUE)

# subclassification
m_out_subclass <- matchit(frml_treat, cps1_nsw_data, method="subclass", subclass=10)

# caliper matching
m_out_caliper <- matchit(frml_treat, cps1_nsw_data, method="nearest", replace=TRUE, caliper=.01, ratio=1)

# weighting methods
w_out_ps <- weightit(frml_treat, cps1_nsw_data, method="ps", estimand="ATT")
w_out_gbm <- weightit(frml_treat, cps1_nsw_data, method="gbm", stop.method="es.mean", estimand="ATT", bag.fraction=0.7, cv.folds=5)
w_out_cbps <- weightit(frml_treat, cps1_nsw_data, method="cbps", estimand="ATT")
summary(w_out_ps)
summary(w_out_gbm)
summary(w_out_cbps)
w_out_ps_trm <- w_out_ps %>% trim(at=5, lower=TRUE)
w_out_gbm_trm <- w_out_gbm %>% trim(at=2, lower=TRUE)
w_out_cbps_trm <- w_out_cbps %>% trim(at=5, lower=TRUE)
summary(w_out_ps_trm)
summary(w_out_gbm_trm)
summary(w_out_cbps_trm)

# Check the covariate balance by love.plot() in "cobalt" package.
# matching models
love.plot(frml_treat, data=cps1_nsw_data,
          weights=data.frame(Nearest_Neighbor_1=get.w(m_out_nn_1),
                             Nearest_Neighbor_4=get.w(m_out_nn_4),
                             Subclassification=get.w(m_out_subclass),
                             Caliper=get.w(m_out_caliper)),
          method=c("matching", "matching", "matching", "matching"),
          binary="std", s.d.denom="treated", grid=FALSE, threshold=.1, limits=c(0,1))

# Weighting models
love.plot(frml_treat, data=cps1_nsw_data,
          weights=data.frame(PS=get.w(w_out_ps_trm),
                             GBM=get.w(w_out_gbm_trm),
                             CBPS=get.w(w_out_cbps_trm)),
          method=c("weighting", "weighting", "weighting"),
          binary="std", s.d.denom="treated", grid=FALSE, threshold=.1, limits=c(0,1))

# check the covariate balance per covariate
bal.plot(m_out_nn_1, var.name="distance", which="both", type="histogram", mirror=TRUE)



# Get the estimated ATT via the linear model.
# matching methods
matched_data <- list(match.data(m_out_nn_1),
                     match.data(m_out_nn_4),
                     match.data(m_out_subclass),
                     match.data(m_out_caliper))
model_label1 <- c("nearest neighbor ratio 1", "nearest neighbor ratio 4", "subclassification", "caliper")

for(i in 1:length(matched_data)){
  print(model_label1[[i]])
  print(lm(re78 ~ treat, data=matched_data[[i]]) %>% broom::tidy())
}

# weighting methods
weighted_data <- list(w_out_ps_trm, w_out_gbm_trm, w_out_cbps_trm)
model_label2 <- c("PS", "GBM", "CBPS")

for(i in 1:length(weighted_data)){
  print(model_label2[[i]])
  print(lm(re78 ~ treat, data=cps1_nsw_data, weights=weighted_data[[i]]$weights) %>% broom::tidy())
}
