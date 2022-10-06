## ----setup, include=FALSE---------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=40),tidy=TRUE)

packages <- c( "data.table","tidyverse","ggplot2","ggExtra","formatR",
               "gridExtra","skimr","here","Hmisc","RColorBrewer",
               "broom","boot", "lmtest")

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


## ---- warning = F, message = F, echo = F, include = F-----------------------------------------------------
file_loc <- url("https://cdn1.sph.harvard.edu/wp-content/uploads/sites/1268/1268/20/nhefs.csv")

nhefs <- read_csv(file_loc) %>% 
  select(qsmk,wt82_71,wt82, wt71, exercise,sex,age,
         race,income, marital,school,
         asthma,bronch, 
         starts_with("alcohol"),
         starts_with("price"),
         starts_with("tax"), 
         starts_with("smoke"),
         smkintensity82_71) %>% 
  mutate(income=as.numeric(income>15),
         marital=as.numeric(marital>2)) %>% 
  na.omit(.)

mod <- lm(wt82~wt71,data=nhefs)

plot1 <- nhefs %>% 
  sample_n(100) %>% 
  select(wt82, wt71) %>% 
  ggplot(.) +
  geom_abline(intercept = coef(mod)[1], slope = coef(mod)[2], col = "blue") +
  geom_segment(aes(y = coef(mod)[1] + coef(mod)[2] * wt71, 
                   yend = wt82, 
                   x = wt71, 
                   xend = wt71), col = "red") +
  geom_point(aes(x = wt71, y = wt82), size = 1, shape = 20) +
  scale_y_continuous(expand=c(0,0), limits=c(40,120)) +
  scale_x_continuous(expand=c(0,1), limits=c(40,120)) +
  xlab("x") + ylab("y")

ggsave(here("figures","2022_02_21-ssr_plot.pdf"), plot=plot1)

x <- seq(-3,3,.5)
plot2 <- ggplot() + 
  xlim(-3, 3) + 
  geom_function(fun = function(x) 7 + (0.57 - x)^2) +
  xlab("f(X)") + ylab("MSE[f(X)]")

ggsave(here("figures","2022_02_21-mse_optimization_plot.pdf"), plot=plot2)

plot3 <- gridExtra::grid.arrange(plot1,plot2,nrow=2)

ggsave(here("figures","2022_02_21-mse_ssr_plot.pdf"), 
       width = 7,
       height = 14,
       units = "cm",
       plot=plot3)



## ----ssrmseplot, out.width="5cm", fig.align='center', fig.margin=TRUE, echo=F, fig.cap="Line of 'best fit' (blue line) defined on the basis of minimizing the sum of squared residuals (red lines) displayed in the top panel; Partial representation of the mean squared error as a function of f(X) in the bottom panel."----
knitr::include_graphics(here("figures","2022_02_21-mse_ssr_plot.pdf"))


## ---- message=F, warning=F--------------------------------------------------------------------------------
nhefs <- read_csv(here("data","nhefs.csv")) %>% 
  mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)))

#' Quick view of data
dim(nhefs)

names(nhefs)


## ---- warning = F, message = F----------------------------------------------------------------------------
#' Here, we start fitting relevant regression models to the data.
#' modelForm is a regression argument that one can use to regress the 
#' outcome (wt_delta) against the exposure (qsmk) and selected confounders.

formulaVars <- "qsmk + sex + age + income + sbp + dbp + price71 + tax71 + race"
modelForm <- as.formula(paste0("wt_delta ~", formulaVars))
modelForm

#' This model can be used to quantify a conditionally adjusted 
#' odds ratio with correct standard error
modelOR <- glm(modelForm,data=nhefs,family = binomial("logit"))

summary(modelOR)

tidy(modelOR)[2,]



## ---- warning = F, message = F----------------------------------------------------------------------------
#' Regress the outcome against the confounders with interaction
ms_model <- glm(modelForm,data=nhefs,family=binomial("logit"))
##' Generate predictions for everyone in the sample to obtain 
##' unexposed (mu0 predictions) and exposed (mu1 predictions) risks.
mu1 <- predict(ms_model,newdata=transform(nhefs,qsmk=1),type="response")
mu0 <- predict(ms_model,newdata=transform(nhefs,qsmk=0),type="response")

#' Marginally adjusted odds ratio
marg_stand_OR <- (mean(mu1)/mean(1-mu1))/(mean(mu0)/mean(1-mu0))
#' Marginally adjusted risk ratio
marg_stand_RR <- mean(mu1)/mean(mu0)
#' Marginally adjusted risk difference
marg_stand_RD <- mean(mu1)-mean(mu0)

#' Using the bootstrap to obtain confidence intervals for the marginally adjusted 
#' risk ratio and risk difference.
bootfunc <- function(data,index){
  boot_dat <- data[index,]
  ms_model <- glm(modelForm,data=boot_dat,family=binomial("logit"))
  mu1 <- predict(ms_model,newdata=transform(boot_dat,qsmk=1),type="response")
  mu0 <- predict(ms_model,newdata=transform(boot_dat,qsmk=0),type="response")
  
  marg_stand_OR_ <- (mean(mu1)/mean(1-mu1))/(mean(mu0)/mean(1-mu0))
  marg_stand_RR_ <- mean(mu1)/mean(mu0)
  marg_stand_RD_ <- mean(mu1)-mean(mu0)
  res <- c(marg_stand_RD_,marg_stand_RR_,marg_stand_OR_)
  return(res)
}

#' Run the boot function. Set a seed to obtain reproducibility
set.seed(123)
boot_res <- boot(nhefs,bootfunc,R=2000)

boot_RD <- boot.ci(boot_res,index=1)
boot_RR <- boot.ci(boot_res,index=2)
boot_OR <- boot.ci(boot_res,index=3)

marg_stand_OR
marg_stand_RR
marg_stand_RD

boot_RD
boot_RR
boot_OR



## ---- warning = F, message = F----------------------------------------------------------------------------
#' Marginal Standardization
##' To avoid assuming no interaction between 
##' quitting smoking and any of the other variables
##' in the model, we subset modeling among 
##' exposed/unexposed. This code removes qsmk from the model,
##' which will allow us to regress the outcome 
##' against the confounders among the exposed and 
##' the unexposed separately. Doing so will allow us 
##' to account for any potential exposure-covariate interactions
##' that may be present. 
formulaVars <- "sex + age + income + sbp + dbp + price71 + tax71 + race"
modelForm <- as.formula(paste0("wt_delta ~", formulaVars))
modelForm

#' Regress the outcome against the confounders 
#' among the unexposed (model0) and then among the exposed (model1)
model0 <- glm(modelForm,data=subset(nhefs,qsmk==0),family=binomial("logit"))
model1 <- glm(modelForm,data=subset(nhefs,qsmk==1),family=binomial("logit"))
##' Generate predictions for everyone in the sample using the model fit to only the 
##' unexposed (mu0 predictions) and only the exposed (mu1 predictions).
mu1 <- predict(model1,newdata=nhefs,type="response")
mu0 <- predict(model0,newdata=nhefs,type="response")

#' Marginally adjusted odds ratio
marg_stand_OR <- (mean(mu1)/mean(1-mu1))/(mean(mu0)/mean(1-mu0))
#' Marginally adjusted risk ratio
marg_stand_RR <- mean(mu1)/mean(mu0)
#' Marginally adjusted risk difference
marg_stand_RD <- mean(mu1)-mean(mu0)

#' Using the bootstrap to obtain confidence intervals for the marginally adjusted 
#' risk ratio and risk difference.
bootfunc <- function(data,index){
  boot_dat <- data[index,]
  model0 <- glm(modelForm,data=subset(boot_dat,qsmk==0),family=binomial("logit"))
  model1 <- glm(modelForm,data=subset(boot_dat,qsmk==1),family=binomial("logit"))
  mu1 <- predict(model1,newdata=boot_dat,type="response")
  mu0 <- predict(model0,newdata=boot_dat,type="response")
  
  marg_stand_OR_ <- (mean(mu1)/mean(1-mu1))/(mean(mu0)/mean(1-mu0))
  marg_stand_RR_ <- mean(mu1)/mean(mu0)
  marg_stand_RD_ <- mean(mu1)-mean(mu0)
  res <- c(marg_stand_RD_,marg_stand_RR_,marg_stand_OR_)
  return(res)
}

#' Run the boot function. Set a seed to obtain reproducibility
set.seed(123)
boot_res <- boot(nhefs,bootfunc,R=2000)

boot_RD <- boot.ci(boot_res,index=1)
boot_RR <- boot.ci(boot_res,index=2)
boot_OR <- boot.ci(boot_res,index=3)

marg_stand_OR
marg_stand_RR
marg_stand_RD

boot_RD
boot_RR
boot_OR



## ----figure1, out.width="10cm", fig.align='center', fig.margin=FALSE, echo=F, fig.cap="Directed acyclic graph depicting confounder adjustment using an outcome modelling approach. In this approach, information from the confounder to the outcome is 'blocked' (blue arrow). With this adjustment, the exposure $X$ os $d$-separated from the outcome $Y$, rendering the estimate of the exposure-outcome association unconfounded."----
knitr::include_graphics(here("figures","2022_03_21-Section3_Figure1.pdf"))


## ----figure3, out.width="10cm", fig.align='center', fig.margin=FALSE, echo=F, fig.cap="Directed acyclic graph depicting the causal relations between variables in a 'pseudo-population' obtained via inverse probability weighting. Again, with this adjustment, the exposure $X$ os $d$-separated from the outcome $Y$, rendering the estimate of the exposure-outcome association unconfounded."----
knitr::include_graphics(here("figures","2022_03_21-Section3_Figure3.pdf"))


## ---- warning = F, message = F----------------------------------------------------------------------------

c <- c(0,0,1,1)
x <- c(1,0,1,0)
n <- c(2,8,9,1)

tibble(x,c,n)



## ---- warning = F, message = F----------------------------------------------------------------------------

# create the propensity score in the dataset
nhefs$propensity_score <- glm(qsmk ~ sex + age + income + sbp + dbp + price71 + tax71 + race, data = nhefs, family = binomial("logit"))$fitted.values

# stabilized inverse probability weights
nhefs$sw <- (mean(nhefs$qsmk)/nhefs$propensity_score)*nhefs$qsmk + 
  ((1-mean(nhefs$qsmk))/(1-nhefs$propensity_score))*(1-nhefs$qsmk)

summary(nhefs$sw)

nhefs %>% select(seqn, qsmk, wt_delta, propensity_score, sw) %>% print(n = 5)



## ---- warning = F, message = F----------------------------------------------------------------------------

model_RD_weighted <- glm(wt_delta ~ qsmk, data = nhefs, weights=sw, family = quasibinomial("identity"))

summary(model_RD_weighted)$coefficients



## ---- warning = F, message = F----------------------------------------------------------------------------
library(lmtest)
library(sandwich)
coeftest(model_RD_weighted, vcov. = vcovHC)

