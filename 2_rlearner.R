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

######################################################################

load(file.path(root, "scratch/wide_df_rlearner.RData"))

df2 = big_df2[W_hat %between% c(0.02, 0.98)] # subset to overlap sample
y    = df2[['y']]
mhat = df2[['fitted_y']]
w    = df2[['w']]
W_hat = df2[['W_hat']]
X = df2[, 6:ncol(df2)] %>% remove_constant() %>% as.matrix()


# %%
tic()
tau_lasso = rlasso(X, w, y, p_hat = W_hat, m_hat = mhat)
toc()

pbnot("rlasso 1 done")

tic()
tau_lasso2 = rlasso(X, w, y, p_hat = W_hat, m_hat = mhat, rs = T)
toc()
##

pbnot("rlasso 2 done")


save(tau_lasso, tau_lasso2, file = file.path(root, "scratch/rlasso_fit.RData"))


# %%
load(file.path(root, "scratch/rlasso_fit.RData"))
# fit_df = data.frame(tau1 = predict(tau_lasso, X),  tau2 = predict(tau_lasso2, X))

cat("ATE:", mean(tau_lasso$tau_hat))
cat("ATE:", mean(tau_lasso2$tau_hat))

coefs = data.table(var = c("intercept", colnames(X)), est1 = round(tau_lasso$tau_beta, 3),
  est2 = round(tau_lasso2$tau_beta, 3))

coefs[order(-abs(est1))][est1>0]
