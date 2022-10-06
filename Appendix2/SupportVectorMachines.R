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


## ---------------------------------------------------------------------------------------------------------

spam <- spam %>% mutate(type=if_else(type=="spam",1,-1))

spam %>% count(type)



## ---- message=F, warning=F, include=T---------------------------------------------------------------------

set.seed(123)
library(e1071)
library(kernlab)

# let's create some data that look like those in Figure 1.1
dat <- as_tibble(data.frame(x1=c(rnorm(9,8,1),rnorm(13,2,1)),
                            x2=c(rnorm(9,3,1),rnorm(13,7,2.5)),
                            y=c(rep(-1,9),rep(1,13))))

dat

# here's a simple svm representation of figure 1.1

svm_simple = svm(as.factor(y) ~ ., data = dat, kernel = "linear", scale = FALSE)

svm_simple

# Plot Results

png("../figures/svm_plot2.png", width=450, height=450, units="px")
plot(svm_simple, dat, symbolPalette = c("red","blue"), color.palette = cm.colors)
dev.off()



## ----fig:svm-simple, out.width = "275px", fig.align='center', fig.cap="SVM Example in Two Dimensions, replicated from page 6 of The 100 Page ML Book, replicated via plot.svm", echo=F----
knitr::include_graphics("../figures/svm_plot2.png")


## ----message=F,warning=F,include=T------------------------------------------------------------------------

p1 <- spam %>% 
  ggplot(.) + 
  theme(text = element_text(size=25)) +
  geom_point(aes(y=you,
                 x=our,
                 group=as.factor(type),
                 color=as.factor(type)),
                 size=5, alpha=.25)

png("../figures/svm_plot3.png", width=450, height=450, units="px")
p1
dev.off()



## ----fig:svm-simple2, out.width = "275px", fig.align='center', fig.cap="More complex data for fitting SVMs, showing spam and nonspam emails as a function of the frequency of `our' and `you' in the email.",echo=F----
knitr::include_graphics("../figures/svm_plot3.png")


## ----message=F,warning=F,include=T------------------------------------------------------------------------

spam2 <- spam %>% select(type,you,our)

svm_polynomial = svm(as.factor(type) ~ ., data = spam2, kernel = "polynomial", gamma = 1, cost = 1)

svm_radial = svm(as.factor(type) ~ ., data = spam2, kernel = "radial", gamma = 1, cost = 1)

# Plot Results
png("../figures/svm_plot4.png", width=450, height=450, units="px")
plot(svm_polynomial, spam2, symbolPalette = c("red","blue"), color.palette = cm.colors)
dev.off()

png("../figures/svm_plot5.png", width=450, height=450, units="px")
plot(svm_radial, spam2, symbolPalette = c("red","blue"), color.palette = cm.colors)
dev.off()



## ---- out.width = "275px", fig.align='center', fig.cap="More complex data for fitting SVMs, showing spam and nonspam emails as a function of the frequency of `our' and `you' in the email, classified with a polynomial separator in SVM",echo=F----
knitr::include_graphics("../figures/svm_plot4.png")


## ---- out.width = "275px", fig.align='center', fig.cap="More complex data for fitting SVMs, showing spam and nonspam emails as a function of the frequency of `our' and `you' in the email, classified with a radial basis function separator in SVM",echo=F----
knitr::include_graphics("../figures/svm_plot5.png")


## ----message=F,warning=F,include=T------------------------------------------------------------------------

svm_radial = svm(as.factor(type) ~ ., data = spam2, kernel = "radial", gamma = 1, cost = 100)

png("../figures/svm_plot6.png", width=450, height=450, units="px")
plot(svm_radial, spam2, symbolPalette = c("red","blue"), color.palette = cm.colors)
dev.off()



## ---- out.width = "275px", fig.align='center', fig.cap="More complex data for fitting SVMs, showing spam and nonspam emails as a function of the frequency of ``our'' and ``you'' in the email, classified with a radial basis function separator in SVM with a cost of 100",echo=F----
knitr::include_graphics("../figures/svm_plot6.png")


## ----message=F,warning=F,include=T------------------------------------------------------------------------

svm_radial = svm(as.factor(type) ~ ., data = spam2, kernel = "radial", gamma = 100, cost = 1)

png("../figures/svm_plot7.png", width=450, height=450, units="px")
plot(svm_radial, spam2, symbolPalette = c("red","blue"), color.palette = cm.colors)
dev.off()



## ---- out.width = "275px", fig.align='center', fig.cap="More complex data for fitting SVMs, showing spam and nonspam emails as a function of the frequency of ``our'' and ``you'' in the email, classified with a radial basis function separator in SVM with a gamma value of 100",echo=F----
knitr::include_graphics("../figures/svm_plot7.png")

