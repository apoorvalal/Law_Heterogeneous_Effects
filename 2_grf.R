# %%
library(LalRUtils)
LalRUtils::libreq(tidyverse, data.table, fst, glue, fastDummies,
                  janitor, tictoc, lubridate,
                  IRdisplay, doMC, RPushbullet,
                  # ml/ CI packages
                  grf, DoubleML, rlearner, glmnet
)

options(scipen=999)
set.seed(42)
root = "~/aa_scratch/motivatedreasoning-main"
pbnot = function(x) pbPost("note", x)

# %%
######   ########  ########
##    ##  ##     ## ##
##        ##     ## ##
##   #### ########  ######
##    ##  ##   ##   ##
##    ##  ##    ##  ##
######   ##     ## ##

# load data

load(file.path(root, "scratch/residualised_data.RData"))
dt_samp = data_main[sample(1:nrow(data_main), 1e5)]
Y = dt_samp[['favors_plaintiff']]; W = dt_samp[['rep']]
# propensity score is share of republican judges in that district X year
W_hat = dt_samp$x_republican/dt_samp$n_judges
district_id = dt_samp[['district']] %>% as.factor %>%  as.numeric
X = dt_samp %>% .[, 9:ncol(.)] %>% remove_constant %>% as.matrix()
dyDummies = model.matrix(~ -1 + as.factor(district_year), dt_samp)

# %%
# residualise on district X year FEs
Y.forest = regression_forest(dyDummies, Y, clusters = district_id, tune.parameters = 'all')
Y.hat = predict(Y.forest)$predictions
pbnot("Partialling out done")

# %% estimate pscore - no longer necessary given pscore is available
# W.forest = regression_forest(dyDummies, W, clusters = district_id)
# W.hat = predict(W.forest)$predictions

# %% # causal forest takes orthogonalized fitted values from above and pscore
cf0 = causal_forest(X, Y, W, tune.parameters = 'all',
                       Y.hat = Y.hat, W.hat = W_hat,
                       clusters = district_id)
cat("ATE basic")
average_treatment_effect(cf0, target.sample='overlap')
pbnot("CF0 done")

# %%
varimp = variable_importance(cf0)
selected= which(varimp > mean(varimp))
cf1 = causal_forest(X[, selected], Y, W, tune.parameters = 'all',
                       Y.hat = Y.hat, W.hat = W_hat,
                       clusters = district_id)
cat("ATE parsimonious")
average_treatment_effect(cf1, target.sample='overlap')
pbnot("CF1 done")

################################################################################

tc <- test_calibration(cf1)
cat("Calibration")
tc

save(cf0, cf1, file = file.path(root, "scratch/causal_forest_fit.rdata"))
pbnot("wrote forests to disk")
# %%

