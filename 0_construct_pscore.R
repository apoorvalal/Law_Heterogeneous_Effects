# %% ####################################################
rm(list = ls())
library(LalRUtils)
LalRUtils::libreq(tidyverse, data.table, fst, fixest, rio, foreach, lubridate,
                  janitor, tictoc, RColorBrewer, patchwork, RPushbullet, IRdisplay)
theme_set(lal_plot_theme()) # add _d() for dark
options(repr.plot.width=12, repr.plot.height=9) ; options(scipen=999)
set.seed(42)
root = "/home/alal/Desktop/WorldBank/projects/motivatedreasoning-main"
# %% ####################################################
df = import(file.path(root, "scratch/DistrictJudgesBIO.dta")) %>% setDT
df[, `:=`(
  ApptYear = year(AppointmentDate),
  TermYear = year(TerminationDate)
)]
setorder(df, courtname, ApptYear)

df[, courtName := str_replace_all(courtname, "u. s. district court+, ", "")]
df[, courtState := str_replace_all(courtName, ".*district of ", "")]
df[, dist := str_sub(district, -2)]

# %% dist indicator checks
df$dist%>% unique %>% sort %>% print
df[, .N, .(dist, courtName)][order(dist)]
# check years of missings
df[dist == "", unique(ApptYear)] %>% sort %>% print
# %% create grid to populate
df2 = df[dist != ""]
dist_year_grid = expand.grid(dist = sort(unique(df2$dist)), yr = 1978:2014) %>% setDT
df3 = df2[, .(dist, songername, ApptYear, TermYear, x_dem, x_republican)]
# %%
tic()
res = list()
for(i in 1:nrow(dist_year_grid)){
  dycomb = dist_year_grid[i]; d = dycomb$dist; y = dycomb$yr
  matches = df3[eval(parse(text = glue::glue("dist == '{d}' & ApptYear < {y} & (TermYear >= {y} | is.na(TermYear))")))]
  d_r_counts = matches[, lapply(.SD, sum), .SDcols = c('x_dem', 'x_republican')]
  res[[i]] = cbind(dycomb, d_r_counts[, n_judges := nrow(matches)])
}
toc()
# %%
dist_yr_counts = rbindlist(res)
dist_yr_counts = dist_yr_counts[n_judges>0]
# %%
dist_yr_counts %>% fwrite(file.path(root, "scratch/district_year_r_d_counts.csv"))
