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


## ---------------------------------------------------------------------------------------------------------
set.seed(123)
n = 250

c1 <- factor(round(rnorm(n),2))
c2 <- factor(round(rnorm(n),2))
c3 <- factor(round(rnorm(n),2))
c4 <- factor(round(rnorm(n),2))
x <- factor(rbinom(n, 1, .5))

dat_ <- data.frame(x,c1,c2,c3,c4)

head(dat_, 10)

mod_mat <- model.matrix(~., data=dat_)

dim(mod_mat)

mod_mat_int <- model.matrix(~.^2, data=dat_)

dim(mod_mat_int)



## ----mlresults1, out.width="12cm", fig.align='center', fig.cap="Absolute bias of inverse probability weighted, g-computation, and doubly robust estimators for sample sizes of  N = 200, N = 1200, and N = 5000. Bar color intensity, from black to light gray, represent IPW, g Computation, AIPW, and TMLE estimators, respectively. Plot labels refer to the following scenarios: Nonpar Complex = nonparametric method fit to the transformed confounders; Nonpar Simple = nonparametric method fit to the untransformed confounders; Par Misspec = parametric method fit to transformed confounders; Par Correct = parametric method fit to untransformed confounders. Parametric regression included logistic regression for the exposure model, and linear regression for the outcome model. Nonparametric method consisted of a stacked generalization with random forests and extreme gradient boosting algorithms, and no sample splitting.", echo=F----
knitr::include_graphics(here("figures", "AJE-00517-2020_Naimi_Figure1_v2.pdf"))


## ---- warning = F, message = F, include = F---------------------------------------------------------------
library(boot)
library(ranger)

nhefs <- read_csv(here("data","nhefs.csv")) %>% 
  mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)),
         age = scale(age),
         sbp = scale(sbp),
         dbp = scale(dbp),
         price71 = scale(price71), 
         tax71 = scale(tax71)) %>% 
  select(-wt82_71) 

#' Marginal Standardization
formulaVars <- "sex + age + income + sbp + dbp + price71 + tax71 + race"
modelForm <- as.formula(paste0("wt_delta ~", formulaVars))

model0 <- glm(modelForm,data=subset(nhefs,qsmk==0),family=binomial("logit"))
model1 <- glm(modelForm,data=subset(nhefs,qsmk==1),family=binomial("logit"))
mu1 <- predict(model1,newdata=nhefs,type="response")
mu0 <- predict(model0,newdata=nhefs,type="response")

marg_stand_RD <- mean(mu1)-mean(mu0)

bootfunc <- function(data,index){
  boot_dat <- data[index,]
  model0 <- glm(modelForm,data=subset(boot_dat,qsmk==0),family=binomial("logit"))
  model1 <- glm(modelForm,data=subset(boot_dat,qsmk==1),family=binomial("logit"))
  mu1 <- predict(model1,newdata=boot_dat,type="response")
  mu0 <- predict(model0,newdata=boot_dat,type="response")
  
  marg_stand_RD_ <- mean(mu1)-mean(mu0)
  return(marg_stand_RD_)
}

#' Run the boot function. Set a seed to obtain reproducibility
set.seed(123)
boot_res <- boot(nhefs,bootfunc,R=2000)
boot_RD <- boot.ci(boot_res)

marg_stand_RD
boot_RD


## ---- warning=F, message=F--------------------------------------------------------------------------------
library(ranger)
#' Marginal Standardization with Random Forest
model0 <- ranger(modelForm, num.trees=500, mtry=3, min.node.size = 50, data=subset(nhefs, qsmk==0))
model1 <- ranger(modelForm, num.trees=500, mtry=3, min.node.size = 50, data=subset(nhefs, qsmk==1))
mu1 <- predict(model1, data=nhefs, type="response")$pred
mu0 <- predict(model0, data=nhefs, type="response")$pred

marg_stand_RDrf <- mean(mu1) - mean(mu0)

bootfunc <- function(data,index){
  boot_dat <- data[index,]
  model0 <- ranger(modelForm, num.trees=500, mtry=5, data=subset(boot_dat, qsmk==0))
  model1 <- ranger(modelForm, num.trees=500, mtry=5, data=subset(boot_dat, qsmk==1))
  mu1 <- predict(model1,data=boot_dat,type="response")$pred
  mu0 <- predict(model0,data=boot_dat,type="response")$pred
  
  marg_stand_RD_ <- mean(mu1)-mean(mu0)
  return(marg_stand_RD_)
}

#' Run the boot function. Set a seed to obtain reproducibility
set.seed(123)
boot_res <- boot(nhefs,bootfunc,R=2000)
boot_RDrf <- boot.ci(boot_res)

marg_stand_RDrf
boot_RDrf$bca

