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

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)


## ---- warning = F, message = F----------------------------------------------------------------------------

library(tmle3)
library(sl3)
library(tidyverse)
library(here)

scale_ <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
nhefs <- read_csv(here("data","nhefs.csv")) %>% 
  mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)),
         age = scale_(age),
         sbp = scale_(sbp),
         dbp = scale_(dbp),
         price71 = scale_(price71), 
         tax71 = scale_(tax71)) %>% 
  select(-wt82_71) 

head(nhefs)

# CREATE SUPERLEARNER LIBRARY
# choose base learners

sl3_list_learners("binomial") 

lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)

# ranger learner
grid_params <- list(num.trees = c(250, 500, 1000, 2000),
                    mtry = c(2,4,6),
                    min.node.size = c(50,100))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
lrnr_ranger <- vector("list", length = nrow(grid))
for(i in 1:nrow(grid)){
  lrnr_ranger[[i]] <- make_learner(Lrnr_ranger, 
                                   num.trees=grid[i,]$num.trees, 
                                   mtry=grid[i,]$mtry,
                                   min.node.size=grid[i,]$min.node.size)
}
########################################################################
 lrnr_ranger <- make_learner(Lrnr_ranger)  ###########  FLAG! ###########
########################################################################

# glmnet learner
grid_params <- seq(0,1,by=.25)
lrnr_glmnet <- vector("list", length = length(grid_params))
for(i in 1:length(grid_params)){
  lrnr_glmnet[[i]] <- make_learner(Lrnr_glmnet, alpha = grid_params[i])
}
########################################################################
 lrnr_glmnet <- make_learner(Lrnr_glmnet)  ###########  FLAG! ###########
########################################################################

# xgboost learner
grid_params <- list(max_depth = c(2, 4, 6, 8),
                    eta = c(0.01, 0.1, 0.2),
                    nrounds = c(50, 100, 500)
                    )
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
lrnr_xgboost <- vector("list", length = nrow(grid))
for(i in 1:nrow(grid)){
  lrnr_xgboost[[i]] <- make_learner(Lrnr_xgboost, max_depth=grid[i,]$max_depth, eta=grid[i,]$eta)
}
########################################################################
 lrnr_xgboost <- make_learner(Lrnr_xgboost)  ##########  FLAG! ##########
########################################################################

# earth learner
grid_params <- c(2,3,4,5,6)
lrnr_earth <- vector("list", length = length(grid_params))
for(i in 1:length(grid_params)){
  lrnr_earth[[i]] <- make_learner(Lrnr_earth, degree = grid_params[i])
}
########################################################################
 lrnr_earth <- make_learner(Lrnr_earth)  ############  FLAG! ############
########################################################################

sl_ <- make_learner(Stack, unlist(list(lrnr_mean, 
                                       lrnr_glm,
                                       lrnr_ranger, 
                                       lrnr_glmnet,
                                       lrnr_xgboost,
                                       lrnr_earth), 
                                      recursive = TRUE))

# DEFINE SL_Y AND SL_A 
# We only need one, because they're the same

Q_learner <- Lrnr_sl$new(learners = sl_, 
                         metalearner = Lrnr_nnls$new(convex=T))
g_learner <- Lrnr_sl$new(learners = sl_, 
                         metalearner = Lrnr_nnls$new(convex=T))
learner_list <- list(Y = Q_learner,
                     A = g_learner)

######################################################################

# PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
ate_spec <- tmle_ATE(treatment_level = 1, control_level = 0)

nodes_ <- list(W = c("age", 
                     "sbp", 
                     "dbp", 
                     "price71", 
                     "tax71", 
                     "sex", 
                     "income",
                     "race"), 
               A = "qsmk", 
               Y = "wt_delta")

# RUN TMLE3 
set.seed(123)
tmle_fit_ <- tmle3(ate_spec, nhefs, nodes_, learner_list)

print(tmle_fit_)

saveRDS(tmle_fit_, here("misc","tmle_fit-preliminary.rds"))

tmle_task <- ate_spec$make_tmle_task(nhefs, nodes_)

initial_likelihood <- ate_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

## save propensity score for diagnosis
propensity_score <- initial_likelihood$get_likelihoods(tmle_task)$A
propensity_score <- propensity_score*nhefs$qsmk + (1-propensity_score)*(1-nhefs$qsmk)

# min and max
print(min(propensity_score))
print(max(propensity_score))

plap_ <- tibble(exposure=nhefs$qsmk,
                pscore=propensity_score)
plap_ <- plap_ %>% mutate(sw = exposure*(mean(exposure)/propensity_score) + 
                               (1-exposure)*((1-mean(exposure))/(1-propensity_score)))

# distribution of PS and stabilized weights
summary(plap_$sw)
summary(propensity_score)

# ps overlap plot
ggplot(plap_) + geom_histogram(aes(pscore,fill=factor(exposure)),
                                   colour="grey50", 
                               alpha=0.75, bins=50,
                               position="identity") +
  scale_fill_manual(values=c("blue","orange")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

ggsave(here("figures","ps_overlap_hist-2022_06_02.pdf"),
       width = 15,
       height = 15, 
       units = "cm")

ggplot(plap_) + geom_density(aes(pscore,color=factor(exposure))) +
  scale_fill_manual(values=c("blue","orange")) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

ggsave(here("figures","ps_overlap_dens-2022_06_02.pdf"),
       width = 15,
       height = 15, 
       units = "cm")

# save outcome predictions for diagnosis
outcome_preds <- initial_likelihood$get_likelihoods(tmle_task)$Y

# super learner coefficients for PS model
g_fit <- tmle_fit_$likelihood$factor_list[["A"]]$learner
g_fit$fit_object$full_fit$learner_fits$Lrnr_nnls_TRUE

# super learner coefficients for outcome model
Q_fit <- tmle_fit_$likelihood$factor_list[["Y"]]$learner
Q_fit$fit_object$full_fit$learner_fits$Lrnr_nnls_TRUE




## ---- warning=F, message=F--------------------------------------------------------------------------------

# remotes::install_github("yqzhong7/AIPW")
library(AIPW)

stacklearner <- Stack$new(lrnr_mean, 
                          lrnr_glm,
                          lrnr_ranger, 
                          lrnr_glmnet,
                          lrnr_xgboost,
                          lrnr_earth) 

metalearner <- Lrnr_nnls$new(convex=T)

sl.lib <- Lrnr_sl$new(learners = stacklearner,
                      metalearner = metalearner)

outcome <- nhefs$wt_delta
exposure <- nhefs$qsmk
covariates <- nhefs[,c("age", 
                     "sbp", 
                     "dbp", 
                     "price71", 
                     "tax71", 
                     "sex", 
                     "income",
                     "race")]

set.seed(123)
AIPW_SL <- AIPW$new(Y = outcome,
                    A = exposure,
                    W = covariates, 
                    Q.SL.library = sl.lib,
                    g.SL.library = sl.lib,
                    k_split = 10,
                    verbose=FALSE)$
  fit()$
  summary(g.bound = 0.025)$ 
  plot.p_score()



## ---- warning= F, message = F-----------------------------------------------------------------------------

library(SuperLearner)

listWrappers()

lrnr_mean <- "SL.mean"
lrnr_glm <- "SL.glm"

# ranger learner
lrnr_ranger = create.Learner("SL.ranger", 
                             tune = list(num.trees = c(250, 500, 1000, 2000),
                                         mtry = c(2,4,6),
                                         min.node.size = c(50,100)))
########################################################################
lrnr_ranger <- NULL; lrnr_ranger$names <- "SL.ranger"  ###########  FLAG! ###########
########################################################################

# glmnet learner
grid_params <- seq(0,1,by=.25)
lrnr_glmnet = create.Learner("SL.glmnet", 
                             tune = list(alpha = grid_params))
########################################################################
lrnr_glmnet <- NULL; lrnr_glmnet$names <- "SL.glmnet"  ###########  FLAG! ###########
########################################################################

# xgboost learner
lrnr_xgboost= create.Learner("SL.xgboost", 
                             tune = list(max_depth = c(2, 4, 6, 8),
                                         eta = c(0.01, 0.1, 0.2),
                                         nrounds = c(50, 100, 500)))
########################################################################
lrnr_xgboost <- NULL; lrnr_xgboost$names <- "SL.xgboost"  ###########  FLAG! ###########
########################################################################

# earth learner
grid_params <- c(2,3,4,5,6)
lrnr_earth = create.Learner("SL.earth", 
                             tune = list(degree = grid_params))
########################################################################
lrnr_earth <- NULL; lrnr_earth$names <- "SL.earth"  ###########  FLAG! ###########
########################################################################

sl.lib <- c(lrnr_mean,lrnr_glm,
            lrnr_ranger$names,lrnr_glmnet$names,
            lrnr_xgboost$names,lrnr_earth$names)

outcome <- nhefs$wt_delta
exposure <- nhefs$qsmk
covariates <- nhefs[,c("age", 
                       "sbp", 
                       "dbp", 
                       "price71", 
                       "tax71", 
                       "sex", 
                       "income",
                       "race")]

set.seed(123)
AIPW_SL <- AIPW$new(Y = outcome,
                    A = exposure,
                    W = covariates, 
                    Q.SL.library = sl.lib,
                    g.SL.library = sl.lib,
                    k_split = 10,
                    save.sl.fit=T,
                    verbose=FALSE)$
  fit()$
  summary(g.bound = 0.025)$ 
  plot.p_score()

#stratified_fit()$

# AIPW_SL$libs$Q.fit
# AIPW_SL$libs$g.fit

print(AIPW_SL$result, digits = 2)


