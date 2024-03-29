## ----setup, include=FALSE---------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=40),tidy=TRUE)

packages <- c( "data.table","tidyverse","ggplot2","ggExtra","formatR","broom",
               "gridExtra","skimr","here","Hmisc","RColorBrewer", "rmutil")

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


## ----F1, out.width = "10cm", fig.align="center", fig.cap="Example dose-response relation between a continuous variable X and an outcome variable Y. The black dashed line represents the true dose-response relation, with all additional lines representing various machine learning estimators of the fit, including the Super Learner (Stacked Regression).",echo=F----
knitr::include_graphics(here("figures", "promo_fig1.pdf"))


## ---- warning = F, message = F----------------------------------------------------------------------------

library(rmutil)
library(here)
library(SuperLearner)

# set the seed for reproducibility
set.seed(123)

# generate the observed data
n=1000
x = runif(n,0,8)
y = 5 + 4*sqrt(9 * x)*as.numeric(x<2) + as.numeric(x>=2)*(abs(x-6)^(2)) + rlaplace(n)

# to plot the true dose-response curve, generate sequence of 'doses' from 0 to 8 at every 0.1,
#	then generate the true outcome
xl<-seq(0,8,.1)
yl<-5 + 4 * sqrt(9 * xl)*as.numeric(xl<2) + as.numeric(xl>=2)*(abs(xl-6)^(2))

D<-data.frame(x,y)  	# observed data
Dl<-data.frame(xl,yl)   # for plotting the true dose-response curve

head(D)



## ---- warning = F, message = F----------------------------------------------------------------------------
# Specify the number of folds for V-fold cross-validation
folds= 5
## split data into 5 groups for 5-fold cross-validation 
## we do this here so that the exact same folds will be used in 
## both the SL fit with the R package, and the hand coded SL
index<-split(1:1000,1:folds)
splt<-lapply(1:folds,function(ind) D[index[[ind]],])
# view the first 6 observations in the first [[1]] and second [[2]] folds
head(splt[[1]])
head(splt[[2]])

#-------------------------------------------------------------------------------
# Fit using the SuperLearner Package
#-------------------------------------------------------------------------------
# Create the 5 df GAMs using functions called from programs "sourced" above
SL.gam.5 <- create.Learner("SL.gam", params = list(deg.gam = 5))

# Specifying the SuperLearner library of candidate algorithms
sl.lib <- c(SL.gam.5$names,"SL.earth")

# Fit using the SuperLearner package, specify 
#		outcome-for-prediction (y), the predictors (x), the loss function (L2),
#		the library (sl.lib), and number of folds 
fitY<-SuperLearner(Y=y,
                   X=data.frame(x), 
                   method="method.NNLS", 
                   SL.library=sl.lib,
                   cvControl=list(V= folds,
                                  validRows=index))

# View the output: 'Risk' column returns the CV-MSE estimates
#		'Coef' column gives the weights for the final SuperLearner (meta-learner)
fitY
# Now predict the outcome for all possible x 'doses'
yS<-predict(fitY,newdata=data.frame(x=xl),onlySL=T)$pred

# Create a dataframe of all x 'doses' and predicted SL responses
Dl1<-data.frame(xl,yS)


## ---- warning=F, message=F--------------------------------------------------------------------------------

library(nnls)
library(earth)
library(gam)

#-------------------------------------------------------------------------------
# Hand-coding Super Learner
#-------------------------------------------------------------------------------

## 2: the lapply() function is an efficient way to rotate through the folds to execute the following:
#	(a) set the ii-th fold to be the validation set; (b) fit each algorithm on the training set, 
#     which is what's left after taking the ii-th fold out (i.e., splt[-ii]); 
# (c) obtain the predicted outcomes for observations in the validation set;
# (d) estimate the estimated risk (CV-MSE) for each fold
#
## 2b: fit each algorithm on the training set (but not the ii-th validation set)
m1<-lapply(1:folds,function(ii) gam(y~s(x,5),family="gaussian",data=rbindlist(splt[-ii])))
m2<-lapply(1:folds,function(ii) earth(y~x,data=rbindlist(splt[-ii]),degree = 2, penalty = 3, 
                                      nk = 21, pmethod = "backward", nfold = 0, 
                                      ncross = 1, minspan = 0, endspan = 0))

## 2c: predict the outcomes for observation in the ii-th validation set
p1<-lapply(1:folds,function(ii) predict(m1[[ii]],newdata=rbindlist(splt[ii]),type="response"))
p2<-lapply(1:folds,function(ii) predict(m2[[ii]],newdata=rbindlist(splt[ii]),type="response"))

# add the predictions to grouped dataset 'splt'
for(i in 1:folds){
	splt[[i]]<-cbind(splt[[i]],p1[[i]],p2[[i]])
}
# view the first 6 observations in the first fold 
#   column2 (y) is the observed outcome;  
#   column3 is the CV-predictions from gam
#   column4 is the CV-predictions from earth; 
head(splt[[1]])

## 2d: calculate CV risk for each method for with ii-th validation set
# our loss function is L2-squared error; so our risk is mean squared error
risk1<-lapply(1:folds,function(ii) mean((splt[[ii]][,2]-splt[[ii]][,3])^2))
risk2<-lapply(1:folds,function(ii) mean((splt[[ii]][,2]-splt[[ii]][,4])^2))

## 3: average the estimated risks across the 5 folds to obtain 1 measure of performance for each algorithm
a<-rbind(cbind("gam",mean(do.call(rbind,risk1),na.rm=T)),
      cbind("earth",mean(do.call(rbind,risk2),na.rm=T)))

# checking to see match with SL output 
fitY
a

#------------
## 4: estimate SL weights using nnls (for convex combination) and normalize

# create a new datafame with the observed outcome (y) and CV-predictions from the 3 algorithms
X<-data.frame(do.call(rbind,splt))[,-1];names(X)<-c("y","gam","earth")
head(X)
SL.r<-nnls(cbind(X[,2],X[,3]),X[,1])$x
alpha<-as.matrix(SL.r/sum(SL.r))
round(alpha,3)
# compare to the package's coefficients 
fitY

#----------------
## 5a: fit all algorithms to original data and generate predictions
m1<-gam(y~s(x,5),family="gaussian",data=D)
m2<-earth(y~x,data=D,degree = 2, penalty = 3, 
          nk = 21, pmethod = "backward", nfold = 0, 
          ncross = 1, minspan = 0, endspan = 0)

## 5b: predict from each fit using all data
p1<-predict(m1,newdata=D,type="response")
p2<-predict(m2,newdata=D,type="response")

predictions<-cbind(p1,p2)

## 5c: for the observed data take a weighted combination of predictions using nnls coeficients as weights
y_pred <- predictions%*%alpha


#--------------------------------------------
## now apply to new dataset (all doses (xl)) to  verify that our work predicts similar results as actual SL function
p1<-predict(m1,newdata=data.frame(x=xl),type="response")
p2<-predict(m2,newdata=data.frame(x=xl),type="response")

predictions<-cbind(p1,p2)
yS2 <- predictions%*%alpha

# now we have a new dataframe of doses (xl) and SL manual predicted outcome
Dl2<-data.frame(xl,yS2)

head(Dl2)

# plotting these data can give us a visual of our estimate 
# of the dose-response function


