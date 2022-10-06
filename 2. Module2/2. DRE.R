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


## ----dag1, fig.margin=TRUE, fig.cap="Causal diagram representing the structure from which the simple simulated data were generated.", echo=F----
knitr::include_graphics(here("figures", "F1.pdf"))


## ---- message=F, warning=F--------------------------------------------------------------------------------

# the inverse logit function, to simulate from logistic model
# (NB: the inverse logit function is the softmax function)
expit <- function(x){(1/(1 + exp(-x)))}
# sample size
n <- 2500

# simulation loop
sim_dat <- function(index){
  set.seed(index)
  # confounders
  c1 <- rbinom(n,1,.5)
  c2 <- rbinom(n,1,.5)
  
  # propensity score model
  pi_x <- expit(-1.5 + log(2.5)*c1 + log(2.5)*c2)
  x <- rbinom(n,1,pi_x)
  
  # outcome model
  mu_y <- expit(-1.5 + log(1.75)*x + log(2.5)*c1 + log(2.5)*c2)
  y <- rbinom(n,1,mu_y)
  
  # simulated data
  a <- data.frame(x,y,c1,c2)
  
  # regression 1: correctly specified 
  mod_true <- glm(y ~ x + c1 + c2, data=a, family=binomial("logit"))
  mu1 <- mean(predict(mod_true, newdata=transform(a,x=1), type="response"))
  mu0 <- mean(predict(mod_true, newdata=transform(a,x=0), type="response"))
  ate_true <- mu1 - mu0
  ate_true
  
  # regression 2: misspecified, one confounder left out
  mod_noconf <- glm(y ~ x + c1, data=a, family=binomial("logit"))
  mu1 <- mean(predict(mod_noconf, newdata=transform(a,x=1), type="response"))
  mu0 <- mean(predict(mod_noconf, newdata=transform(a,x=0), type="response"))
  ate_noconf <- mu1 - mu0
  ate_noconf
  
  # regression 3: double robust estimation
  # correctly specified PS model
  a$propensity_score <- glm(x ~ c1 + c2, data = a, family = binomial("logit"))$fitted.values
  # stabilized inverse probability weights
  a$sw <- (mean(a$x)/a$propensity_score)*a$x + ((1-mean(a$x))/(1-a$propensity_score))*(1-a$x)
  
  # misspecified outcome model
  mod_dr1 <- glm(y ~ x + c1, weights=sw, data=a, family=binomial("logit"))
  mu1 <- mean(predict(mod_dr1, newdata=transform(a,x=1), type="response"))
  mu0 <- mean(predict(mod_dr1, newdata=transform(a,x=0), type="response"))
  ate_dr1 <- mu1 - mu0
  ate_dr1

  # regression 5: double robust estimation  
  # misspecified PS model
  a$propensity_score <- glm(x ~ c1, data = a, family = binomial("logit"))$fitted.values
  # stabilized inverse probability weights
  a$sw <- (mean(a$x)/a$propensity_score)*a$x + ((1-mean(a$x))/(1-a$propensity_score))*(1-a$x)

  # correctly specified outcome model
  mod_dr2 <- glm(y ~ x + c1 + c2, weights=sw, data=a, family=binomial("logit"))
  mu1 <- mean(predict(mod_dr1, newdata=transform(a,x=1), type="response"))
  mu0 <- mean(predict(mod_dr1, newdata=transform(a,x=0), type="response"))
  ate_dr2 <- mu1 - mu0
  ate_dr2
  
  res <- data.frame(
    rbind(c(Method="True", Estimate = ate_true),
          c(Method="Confounded",      Estimate = ate_noconf),
          c(Method="Double Robust 1",   Estimate = ate_dr1),
          c(Method="Double Robust 2",   Estimate = ate_dr2))
    )
  
  return(res)
}

res <- lapply(1:500, function(x) sim_dat(index=x))

res <- do.call(rbind, res)

res[,2] <- as.numeric(res[,2])

plot1 <- ggplot(res) +
  geom_density(aes(x = Estimate,
                     group = Method,
                     fill = Method),
                 alpha=.2) +
  scale_x_continuous(expand=c(0,0), 
                     limits=c(0,.3)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_discrete(name="") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(here("figures","simulation_results.pdf"),
       plot = plot1,
       width = 10,
       height = 10,
       units = "cm")



## ----simfigure, out.width="6cm", fig.align='center', fig.cap="Distribution of 500 Estimates from the Simple Simulation Illustration of DR Estimation", echo=F----
knitr::include_graphics(here("figures","simulation_results.pdf"))


## ----drs, fig.fullwidth = TRUE, fig.margin=FALSE, fig.cap=" Comparisons of R Packages for Implementing Doubly Robust Estimators.", echo=F----
knitr::include_graphics(here("figures", "DRSoftware.pdf"))


## ---- warning = F, message = F, eval = F------------------------------------------------------------------
## 
## remotes::install_github("tlverse/tlverse")
## library(tmle3)
## library(sl3)
## 


## ---- warning = F, message = F----------------------------------------------------------------------------
library(tidyverse)

remotes::install_github("tlverse/tlverse")
library(tmle3)
library(sl3)

# Import NHEFS data
nhefs <- read_csv(here("data","nhefs.csv")) %>% 
  mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)))

# Quick view of data
dim(nhefs)
names(nhefs)

# using the sl3 function create the super learner, that contains a simple
# average estimator and a glm estimator.
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_mean <- make_learner(Lrnr_mean)
sl <- Lrnr_sl$new(learners = list(lrnr_glm,lrnr_mean))

learner_list <- list(Y = sl, A = sl)

######################################################################

# PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
# ate_spec defines the effect we want to target, and the exposure levels we want 
# to compare
ate_spec <- tmle_ATE(treatment_level = 1, control_level = 0)
# nodes define the variables we want to adjsut for, and the exposure and the outcome
nodes <- list(W = c("sex", "age", "income", "sbp", "dbp", "price71", "tax71", "race"), 
              A = "qsmk", 
              Y = "wt_delta")

# Run the tmle3 function using the pieces defined above 
tmle_fit <- tmle3(ate_spec, nhefs, nodes, learner_list)

tmle_fit$summary$psi_transformed

tmle_fit$summary$lower_transformed

tmle_fit$summary$upper_transformed



## ---- warning = F, messages = F, results='hide'-----------------------------------------------------------

install.packages("AIPW",
                 repos='http://lib.stat.cmu.edu/R/CRAN',
                 dependencies=T)

library(AIPW)

install.packages("SuperLearner",
                 repos='http://lib.stat.cmu.edu/R/CRAN',
                 dependencies=T)

library(SuperLearner)



## ---- warning = F, message = F----------------------------------------------------------------------------

# we will see, the AIPW function includes important visual diagnostics
# which can be implemented with the ggplot library
library(ggplot2)

# Import NHEFS data, again
nhefs <- read_csv(here("data","nhefs.csv")) %>% 
  mutate(wt_delta = as.numeric(wt82_71>median(wt82_71)))

# using the SuperLearner function with simple average estimator and a glm estimator.

sl <- c("SL.mean","SL.glm")

outcome <- nhefs$wt_delta
exposure <- nhefs$qsmk
covariates <- nhefs[,c("sex", "age", "income", "sbp", "dbp", "price71", "tax71", "race")]

AIPW_SL <- AIPW$new(Y = outcome,
                    A = exposure,
                    W = covariates,
                    Q.SL.library = sl,
                    g.SL.library = sl,
                    k_split = 3,
                    verbose=FALSE)$
  fit()$
  summary(g.bound = 0.025)$ 
  plot.p_score()$
  plot.ip_weights()

print(AIPW_SL$result[3,], digits=3)


