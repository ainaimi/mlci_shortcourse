## ----setup, include=FALSE---------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=40),tidy=TRUE)

packages <- c( "data.table","tidyverse","ggplot2","ggExtra","formatR","broom",
               "gridExtra","skimr","here","Hmisc","RColorBrewer")

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


## ---- out.width = "10cm", fig.align="center", fig.cap="Classification and Regression Tree Algorithms fit to the NHEFS Data with a complexity parameter of 0.01 (left panel) and 0.005 (right panel).",echo=F----
knitr::include_graphics(here("figures", "cart_tuning.pdf"))


## ---- out.width = "10cm", fig.align="center", fig.cap="Nine Classification and Regression Tree Algorithms fit to bootstrap resamples of the NHEFS Data under default tuning parameters.",echo=F----
knitr::include_graphics(here("figures","rf_example.pdf"))


## ---------------------------------------------------------------------------------------------------------

library(tidyverse)
library(ranger)

nhefs <- read_csv(here("data","nhefs.csv")) %>% 
  mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)),
         age = scale(age),
         sbp = scale(sbp),
         dbp = scale(dbp),
         price71 = scale(price71), 
         tax71 = scale(tax71)) %>% 
  select(-wt82_71) 

head(nhefs)

# implement ranger
# set random seed
set.seed(123)
rf_example <- ranger(factor(wt_delta) ~ qsmk + sex + 
                       age + income + sbp + dbp + 
                       price71 + tax71 +race,
                     num.trees=500, 
                     mtry=3, 
                     min.node.size = 50, 
                     data=nhefs)

rf_example

