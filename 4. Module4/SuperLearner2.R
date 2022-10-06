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


## ---- warning=F, message=F, tidy=T------------------------------------------------------------------------

library(sl3)

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

# Begin modeling nhefs with super learner
# create the prediction task (i.e., use nhefs to predict outcome)
task <- make_sl3_Task(
  data = nhefs,
  outcome = "wt_delta",
  covariates = c("qsmk","age","sbp","dbp","price71","tax71","sex","income","race"),
  folds=5
)

# let's look at the task
task



## ---- warning=F, message=F, tidy=T------------------------------------------------------------------------

sl3_list_properties()



## ---- warning=F, message=F, tidy=T------------------------------------------------------------------------

sl3_list_learners(properties = "binomial")



## ---- warning=F, message=F, tidy=T------------------------------------------------------------------------

# create a simple glm learner and a simple mean learner
## note: no change in any tuning parameters
lrn_glm <- Lrnr_glm$new()
lrn_mean <- Lrnr_mean$new()

# create ridge and lasso regression:
## note: rigde and lasso are defined by setting alpha tuning parameter in the glmnet function
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_lasso <- Lrnr_glmnet$new(alpha = 1)

# create xgboost and ranger:
## note: using default tuning parameters
lrn_ranger <- Lrnr_ranger$new()
lrn_xgb <- Lrnr_xgboost$new()



## ---- warning=F, message=F, tidy=T------------------------------------------------------------------------
stack <- Stack$new(
  lrn_glm, lrn_mean, lrn_ridge, lrn_lasso, lrn_ranger, lrn_xgb
)
stack


## ---- warning=F, message=F--------------------------------------------------------------------------------

lrn_ridge



## ---- warning = F, message = F----------------------------------------------------------------------------

sl <- Lrnr_sl$new(
  learners = stack, 
  metalearner = Lrnr_nnls$new(convex=T)
  )



## ---- warning = F, message = F----------------------------------------------------------------------------

set.seed(123)
sl_fit <- sl$train(task = task)

sl_fit$fit_object$cv_meta_fit



## ---- warning=F, message=F, eval=F------------------------------------------------------------------------
## lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
## lrn_lasso <- Lrnr_glmnet$new(alpha = 1)


## ---- warning=F, message=F--------------------------------------------------------------------------------

# glmnet learner
grid_params <- seq(0,1,by=.1)
lrnr_glmnet <- vector("list", length = length(grid_params))
for(i in 1:length(grid_params)){
  lrnr_glmnet[[i]] <- make_learner(Lrnr_glmnet, alpha = grid_params[i])
}



## ---------------------------------------------------------------------------------------------------------

# first four elements of the new lrnr_glmnet object
lrnr_glmnet[1:4]



## ---- warning=F, message=F--------------------------------------------------------------------------------
# create the learning stack with the new lrnr_glmnet library
stack <- Stack$new(
  lrnr_glmnet
)

# define the metalearner
sl <- Lrnr_sl$new(
  learners = stack, 
  metalearner = Lrnr_nnls$new(convex=T)
  )

# run the super learner and print the results
set.seed(123)
sl_fit <- sl$train(task = task)

sl_fit$fit_object$cv_meta_fit


## ---------------------------------------------------------------------------------------------------------
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

# create the learning stack with the new lrnr_ranger library
stack <- Stack$new(
  lrnr_ranger
)
stack

# define the metalearner
sl <- Lrnr_sl$new(
  learners = stack, 
  metalearner = Lrnr_nnls$new(convex=T)
  )

# run the super learner and print the results
set.seed(123)
sl_fit <- sl$train(task = task)

sl_fit$fit_object$cv_meta_fit

