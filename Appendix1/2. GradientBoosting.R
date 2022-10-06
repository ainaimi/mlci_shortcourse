## ----setup, include=FALSE---------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=40),tidy=TRUE)

packages <- c( "data.table","tidyverse","ggplot2","ggExtra","formatR","broom",
               "gridExtra","skimr","here","Hmisc","RColorBrewer")#,"gmm")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN',dependencies=T)
  }
}

for (package in packages) {
  library(package, character.only=T)
}

remotes::install_github("rstudio/fontawesome")

library(fontawesome)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)


## ---- warning = F, message = F----------------------------------------------------------------------------

library(tidyverse)
library(xgboost)
library(caret)
library(DiagrammeR)

nhefs <- read_csv(here("data","nhefs.csv")) %>% 
  mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)))

head(nhefs)

# create tuning grid
grid_default <- expand.grid(nrounds = c(250),
                            max_depth = c(4),
                            eta = c(0.01),
                            gamma = c(0,1),
                            min_child_weight = c(10,25),
                            colsample_bytree = c(0.7),
                            subsample = c(0.6))
# set random seed
set.seed(123)
# train XGBoost model
xgboost1 <- train(factor(wt_delta) ~ qsmk + sex + age + income + sbp + dbp + price71 + tax71 +race, 
                      data = nhefs,
                      tuneGrid = grid_default,
                      method = "xgbTree",
                      metric = "Kappa")

tree_plot <- xgb.plot.tree(model = xgboost1$finalModel, trees = 1:2, plot_width=2000, plot_height=2000)


## ----xgbeg, out.width = "10cm", fig.align="center", fig.cap="Two trees from the example extreme gradient boosting algorithm fit to the NHEFS data", echo=F----
knitr::include_graphics(here("figures", "xgboost_example.pdf"))

