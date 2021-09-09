# %%
rm(list = ls())
library(LalRUtils)
LalRUtils::libreq(tidyverse, data.table, fst, fixest, rio, glue, fastDummies,
                  janitor, tictoc, patchwork, lubridate, cobalt, modelsummary,
                  IRdisplay, knitr, rmarkdown, kableExtra, tictoc, doMC,
                  # ml/ CI packages
                  grf, DoubleML, npcausal, rlearner, glmnet
)
theme_set(lal_plot_theme()) # add _d() for dark
options(repr.plot.width=12, repr.plot.height=9)
options(ggplot2.discrete.fill = RColorBrewer::brewer.pal(9, "Set1"))
options(ggplot2.discrete.colour = RColorBrewer::brewer.pal(9, "Set1"))
options(ggplot2.continuous.fill = "viridis"); options(ggplot2.continuous.colour = "viridis")
options(scipen=999)

set.seed(42)
chr = function(...) as.character(...) %>% display_html()
root = "/home/alal/Desktop/WorldBank/projects/motivatedreasoning-main"

# %%
########  ########  ######## ########
##     ## ##     ## ##       ##     ##
##     ## ##     ## ##       ##     ##
########  ########  ######   ########
##        ##   ##   ##       ##
##        ##    ##  ##       ##
##        ##     ## ######## ##
tic()
df = fread(file.path(root, "scratch/new_cases.csv.tar.gz")) %>% clean_names
toc()

dist_yr_pscore = fread(file.path(root, "scratch/district_year_r_d_counts.csv"))

# %%
df[, courtType := case_when(
  str_detect(tolower(court_name), "circuit") ~  "Circuit",
  str_detect(tolower(court_name), "appeals") ~  "Appeals",
  str_detect(tolower(court_name), "district") ~  "District",
  TRUE ~ court_name
)]

# %% prep
df = df[courtType == "District"]
dropcols = c('decision_date', 'termdate', 'tapeyear')
df[, (dropcols) := NULL]
# subset on outcome
df = df[judgment %in% c(1, 2)]
df[, favors_plaintiff := ifelse(judgment == 2, 1, 0)]
# drop marginal parties
setnames(df, 'party_of_appointing_president', 'party')
df = df[party %in% c('Republican', 'Democratic')]
df[, rep := ifelse(party == "Republican", 1, 0)]
# timestamp
df[, filing_year := year(filing_date)]
cat(glue("Missing share: {format(sum(is.na(df$filing_year))/nrow(df), 3)}"))
# district year FE
df[, district_year := .GRP, by =  .(district, filing_year)]
# %% merge in pscore by district
df = merge(df, dist_yr_pscore, by.x = c('district', 'filing_year'), by.y = c('dist', 'yr'))
# %% # number of unique values for covariates
X = c('classact', 'juris', 'origin', 'office', 'nos', 'residenc')
df[, lapply(.SD, nunique), .SDcols = X]
# coerce to factor before constructing dummies
df[, (X) := lapply(.SD, function(x) as.factor(x)), .SDcols = X]
# drop rows with missing
kvs = c('favors_plaintiff', 'rep',  X, 'filing_year', 'district_year',
  'district', 'x_republican', 'x_dem', 'n_judges')
df = df[, ..kvs] %>% na.omit()

# %% dummy columns
Xmat = dummy_cols(df[, ..X], select_columns = X, remove_first_dummy = TRUE,
  remove_selected_columns = TRUE)
data_main = cbind(df[, .(favors_plaintiff, rep, district, district_year,
  filing_year, x_republican, x_dem, n_judges)], Xmat)
fwrite(data_main, file.path(root, "scratch/prepped_data.csv"))

# %%
cat("Tabulations of categorical variables")
summary(df[, ..X])
data_main %>% glimpse
# %%
partialer = function(v, fitted = F){
  # subset dataframe
  kv = c(v, 'district_year'); dat = data_main[, ..kv]
  f = formula_fixest(v, X = "1", D = "district_year")
  m = feols(f, data = dat)
  if (fitted == F){
   return(m$residuals)
  } else{
   return(m$fitted.values)
  }
}
# residualise treatment and outcome
resid_y  = partialer('favors_plaintiff')
resid_d  = partialer('rep')
# residualise covariates
tic()
Xvars = colnames(Xmat)
resid_X = map(Xvars, partialer)
toc()
# construct dataframe partialled out
resid_X = resid_X  %>% setDT
resid_df = cbind(resid_y, resid_d, resid_X, data_main$district)
setnames(resid_df, c("y", "d", c(colnames(Xmat), 'district')))
save(resid_df, data_main, file = file.path(root, "scratch/residualised_data.RData"))

# %% # prep for Rlearner - fit nuisance functions
fitted_y  = partialer('favors_plaintiff', T)
fitted_d  = partialer('rep', T)
big_df2 = data.table(
 y = data_main[['favors_plaintiff']],
 w = data_main[['rep']],
 fitted_y,
 fitted_d,
 data_main[, ..Xvars]
)
save(big_df2, file = file.path(root, "scratch/wide_df_rlearner.RData"))


# %%
load(file.path(root, "scratch/residualised_data.RData"))
m0 = feols(y ~ d-1, data = resid_df, cluster = ~district)
etable(m0)


# %%
########  ##     ## ##
##     ## ###   ### ##
##     ## #### #### ##
##     ## ## ### ## ##
##     ## ##     ## ##
##     ## ##     ## ##
########  ##     ## ########


load(file.path(root, "scratch/residualised_data.RData"))
dml_data = DoubleMLData$new(resid_df,
                            y_col = "y",
                            d_cols = "d",
                            x_cols = colnames(Xmat))
dml_data
library(mlr3)
library(mlr3learners)
# surpress messages from mlr3 package during fitting
lgr::get_logger("mlr3")$set_threshold("warn")
# random forest for nuisance functions takes too long
# learner = lrn("regr.ranger", num.trees=500, max.depth=5, min.node.size=2)
# lasso
learner = lrn("regr.glmnet")
# set to use 4 CPUs
set_threads(learner, n = 4)
ml_g_bonus = learner$clone()
ml_m_bonus = learner$clone()
set.seed(94305)
dml_1 = DoubleMLPLR$new(dml_data, ml_g=ml_g_bonus, ml_m=ml_m_bonus)
dml_1$fit()
dml_1$summary()
save(dml_1, file = file.path(root, "scratch/plr1.RData"))


# %%
######   ########  ########
##    ##  ##     ## ##
##        ##     ## ##
##   #### ########  ######
##    ##  ##   ##   ##
##    ##  ##    ##  ##
######   ##     ## ##


load(file.path(root, "scratch/residualised_data.RData"))
dt_samp = data_main[sample(1:nrow(data_main), 1e5)]

Y = dt_samp[['favors_plaintiff']]
W = dt_samp[['rep']]
W_hat = dt_samp$x_republican/dt_samp$n_judges
district_id = dt_samp[['district']] %>% as.factor %>%  as.numeric
X = dt_samp %>% .[, 9:ncol(.)] %>% remove_constant %>% as.matrix()
dyDummies = model.matrix(~ -1 + as.factor(district_year), dt_samp)


# %%
# residualise on district X year FEs
Y.forest = regression_forest(dyDummies, Y, clusters = district_id,
  tune.parameters = 'all')
Y.hat = predict(Y.forest)$predictions
# %% estimate pscore - no longer necessary given pscore is available
# W.forest = regression_forest(dyDummies, W, clusters = district_id)
# W.hat = predict(W.forest)$predictions

# %% # causal forest takes orthogonalized fitted values from above and pscore
cf0 = causal_forest(X, Y, W,
                       Y.hat = Y.hat, W.hat = W_hat,
                       clusters = district_id)
average_treatment_effect(cf0, target.sample='overlap')
# %%
varimp = variable_importance(cf0)
selected= which(varimp > mean(varimp))
cf1 = causal_forest(X[, selected], Y, W,
                       Y.hat = Y.hat, W.hat = W_hat,
                       clusters = district_id)
average_treatment_effect(cf1, target.sample='overlap')
# %%
save(cf0, cf1, file = file.path(root, "scratch/causal_forest_fit.rdata"))

# %%
#######  ##     ## ######## ########
##     ## ##     ##    ##    ##     ##
##     ## ##     ##    ##    ##     ##
##     ## ##     ##    ##    ########
##     ## ##     ##    ##    ##
##     ## ##     ##    ##    ##
#######   #######     ##    ##

load(file.path(root, "scratch/causal_forest_fit.rdata"))

# %%
oob_pred <- predict(cf1)
oob_tauhat_cf <- oob_pred$predictions
hist(oob_tauhat_cf, main="Causal forests: out-of-bag CATE",
  col = "cornflowerblue", las = 1)

# %%
average_treatment_effect(cf1, target.sample = 'overlap')

# %%
tic()
tc <- test_calibration(cf1)
toc()

tc
# %%
save(tc, file = file.path(root, "scratch/omnibus.RData"))


# %%
var_imp <- c(variable_importance(cf1))
names(var_imp) <- colnames(cf1$X.orig)
sorted_var_imp <- sort(var_imp, decreasing=TRUE)
as.data.frame(sorted_var_imp, row.names = names(sorted_var_imp)) %>%
  slice(1:10) %>%
  kable("html", digits = 4, row.names = T) %>%
  chr

# %%
# 440 Other Civil Rights
# 441 Voting
# 442 Employment
# 443 Housing/Accommodations
# 444 Welfare 1
# 445 Amer w/Disabilities-Employment
# 446 Amer w/Disabilities - Other
# 448 Education
# LABOR
# 710 Fair Labor Standards Act
# 720 Labor/Management Relations
# 730 Labor/Management Reporting & Disclosure Act 1
# 740 Railway Labor Act
# 751 Family and Medical Leave Act
# 790 Other Labor Litigation
# 791 Employee Retirement Income Security Act
# PROPERTY RIGHTS
# 820 Copyrights
# 830 Patent
# 840 Trademark
# OTHER STATUTES
# 410 Antitrust
# 892 Economic Stabilization Act 1
# 893 Environmental Matters
# 950 Constitutionality of State Statutes
# %%
best_linear_projection(cf0, subset = W_hat < 1 & W_hat > 0, cf0$X.orig)


# %%
load(file.path(root, "scratch/wide_df_rlearner.RData"))
df2 = big_df2 #[sample(.N, 1e5)]
y    = df2[['y']]
mhat = df2[['fitted_y']]
w    = df2[['w']]
what = df2[['fitted_d']]
X = df2[, 5:ncol(df2)] %>% remove_constant() %>% as.matrix()


# %%
## tic()
## tau_lasso = rlasso(X, w, y, p_hat = what, m_hat = mhat)
## toc()
##
## tic()
## tau_lasso2 = rlasso(X, w, y, p_hat = what, m_hat = mhat, rs = T)
## toc()
##
## save(tau_lasso, tau_lasso2, file = file.path(root, "scratch/rlasso_fit.RData"))


# %%
load(file.path(root, "scratch/rlasso_fit.RData"))
# fit_df = data.frame(tau1 = predict(tau_lasso, X),  tau2 = predict(tau_lasso2, X))
cat("ATE:", mean(tau_lasso$tau_hat))
cat("ATE:", mean(tau_lasso2$tau_hat))

coefs = data.table(var = c("intercept", colnames(X)), est1 = round(tau_lasso$tau_beta, 3),
  est2 = round(tau_lasso2$tau_beta, 3))

coefs[order(-abs(est1))][est1>0]
