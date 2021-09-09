# %%
library(LalRUtils)
LalRUtils::libreq(tidyverse, data.table, fst, glue, fastDummies,
                  janitor, tictoc, lubridate,
                  IRdisplay, doMC, RPushbullet,
                  grf, DoubleML, rlearner, glmnet, fixest
)
options(scipen=999)
set.seed(42)
root = "~/aa_scratch/motivatedreasoning-main"
pbnot = function(x) pbPost("note", x)

load(file.path(root, "scratch/causal_forest_fit.rdata"))

################################################################################
# ATE
average_treatment_effect(cf0, target = "overlap")
average_treatment_effect(cf1, target = "overlap")
oob_pred <- predict(cf1)
oob_tauhat_cf <- oob_pred$predictions
hist(oob_tauhat_cf, main="Causal forests: out-of-bag CATE",
  col = "cornflowerblue", las = 1)

################################################################################
# overlap subpop
#trimmed_ind = which(cf0$W.hat > 0.05 & cf0$W.hat < 0.95)

#tauhat = oob_tauhat_cf[trimmed_ind]
#Xorig =  cf0$X.orig[trimmed_ind, ]
#dfnew = data.frame(tauhat, Xorig)

var_imp <- c(variable_importance(cf0))
names(var_imp) <- colnames(cf1$X.orig)
sorted_var_imp <- sort(var_imp, decreasing=TRUE)
as.data.frame(sorted_var_imp, row.names = names(sorted_var_imp)) %>%
  slice(1:10)

var_imp <- c(variable_importance(cf1))
names(var_imp) <- colnames(cf1$X.orig)
sorted_var_imp <- sort(var_imp, decreasing=TRUE)
as.data.frame(sorted_var_imp, row.names = names(sorted_var_imp)) %>%
  slice(1:10)




