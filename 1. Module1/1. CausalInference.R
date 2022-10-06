## ---------------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=40),tidy=TRUE)

packages <- c( "data.table","tidyverse","ggplot2","ggExtra","formatR",
               "gridExtra","skimr","here","Hmisc","RColorBrewer", "MatchIt")#,"gmm")

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


## --------------------------------------------------
remotes::install_github("yqzhong7/AIPW")
library(AIPW)

install.packages("SuperLearner",repos = "https://cloud.r-project.org/", dependencies=TRUE)
library(SuperLearner)

# define the expit function
expit<-function(z){1/(1+exp(-(z)))}
set.seed(123)
n<-1e6
confounder<-rbinom(n,1,.5)
smoking<-rbinom(n,1,expit(-2+log(2)*confounder))
CVD<-rbinom(n,1,.1+.05*smoking+.05*confounder)

# the data
head(data.frame(CVD,smoking,confounder))

round(mean(confounder),3)
round(mean(smoking),3)
round(mean(CVD),3)

#OLS
round(coef(lm(CVD~smoking+confounder)),4)

#ML1
round(coef(glm(CVD~smoking+confounder,family=poisson("identity"))),4)

#ML2
round(coef(glm(CVD~smoking+confounder,family=binomial("identity"))),4)

#GMM
# round(gmm(CVD~smoking+confounder,x=cbind(smoking, confounder))$coefficients,4)

#AIPW
AIPW_SL <- AIPW$new(Y = CVD,
                    A = smoking,
                    W = confounder, 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 3,
                    verbose=FALSE)$
  fit()$
  summary()

round(AIPW_SL$result[3,1],4)


## ----------------------------------------------------------------------------------------------------
## for calculations below
set.seed(123)
n<-1e6;confounder<-rbinom(n,1,.5)
smoking<-rbinom(n,1,expit(-2+log(2)*confounder))
CVD<-rbinom(n,1,.1+.05*smoking+.05*confounder)
python_data <- data.frame(cbind(CVD,smoking,confounder))
python_data_x <- data.frame(cbind(smoking,confounder))
python_data_y <- data.frame(CVD)
pC<-round(mean(confounder),3)
ols_RD<-round(coef(lm(CVD~smoking+confounder)),4)


## ------------------------------------------------------------------------------

library(MatchIt)
data("lalonde")

head(lalonde)

propensity_score <- glm(treat ~ age + educ + re74 + re75, data=lalonde, family=binomial(link="logit"))$fitted.values

head(propensity_score)

## by appending a "$fitted.values" to the end
## of this glm function, we are 
## keeping the predicted values from the model 
## under the observed data settings.



## ------------------------------------------------------------------------------

set.seed(123)

exposure <- lalonde$treat

plot_data <- data.frame(propensity_score,Exposure=as.factor(lalonde$treat))

p1 <- ggplot(data=plot_data) + 
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  ylab("Density") +
  xlab("Propensity Score") +
  scale_color_manual(values=c("#000000","#D55E00")) +
  geom_density(aes(x=propensity_score,
                   group=Exposure,
                   color=Exposure)) + 
  geom_histogram(aes(y = ..density.., 
                     x=propensity_score,
                     alpha=.25,
                   group=Exposure,
                   color=Exposure)) +
  xlim(0,1)

ggsave(here("figures", "2022_01_10-ps_overlap.pdf"), plot=p1)



## ----psfigure, out.width="10cm", fig.align='center', fig.cap="Propensity score overlap plot for the training intervention in 614 individuals in the Lalonde dataset.", echo=F----
knitr::include_graphics(here("figures", "2022_01_10-ps_overlap.pdf"))


## ----------------------------------------------------------------------------------------------------

sw <- (mean(exposure)/propensity_score)*exposure + 
  ((1 - mean(exposure))/(1 - propensity_score))*(1 - exposure)

summary(sw)


