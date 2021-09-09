# %%
rm(list = ls())
library(LalRUtils)
LalRUtils::libreq(tidyverse, data.table, fst, fixest, rio, glue, fastDummies,
                  janitor, tictoc, patchwork, lubridate, cobalt,
                  IRdisplay, tictoc, doMC, RPushbullet,
                  # ml/ CI packages
                  grf, DoubleML, rlearner, glmnet
)
theme_set(lal_plot_theme()) # add _d() for dark
options(repr.plot.width=12, repr.plot.height=9)
options(scipen=999)
set.seed(42)

root = "~/aa_scratch/motivatedreasoning-main"
pbnot = function(x) pbPost("note", x)


# %%
########  ########  ######## ########
##     ## ##     ## ##       ##     ##
##     ## ##     ## ##       ##     ##
########  ########  ######   ########
##        ##   ##   ##       ##
##        ##    ##  ##       ##
##        ##     ## ######## ##

tic()
df = fread(file.path(root, "scratch/new_cases.csv")) %>% clean_names
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
# propensity score is share of republican judges in that district X year
W_hat = data_main$x_republican/data_main$n_judges

big_df2 = data.table(
 y = data_main[['favors_plaintiff']],
 w = data_main[['rep']],
 W_hat,
 fitted_y,
 fitted_d,
 data_main[, ..Xvars]
)


save(big_df2, file = file.path(root, "scratch/wide_df_rlearner.RData"))

# %%

m0 = feols(y ~ d-1, data = resid_df, cluster = ~district)
etable(m0)


# %%

