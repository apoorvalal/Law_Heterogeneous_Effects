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

