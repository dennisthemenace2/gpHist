## Overview


This package provides a Gaussian Process with histogram intersection kernel and variance approximations.
The package is design to provide a high speed for estimation of the Gaussian Process. Therefore, it only encompasses a very limited number of functions and is supposed to be light weight.

This source code is part of the publication:
Becker, D. (2017). Analysis of the Histogram Intersection Kernel for Use in Bayesian Optimization, 7(6), 337–345. doi:10.7763/IJMO.2017.V7.609

## Installation

```{r, eval = FALSE}

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("dennisthemenace2/gpHist")
```

If you encounter a clear bug, please file a minimal reproducible example.

## Usage

```{r, message = FALSE}
library(gpHist)

#define a function
testFn = function(x){
  y = sin(2*pi*x*2) 
}

#Get data
X = seq(0,1,0.1)
Y = testFn(X)

#Call gpHist function
gp_hist = gpHist(matrix(X),matrix(Y),sigma=0.01)


x_pred = matrix(seq(0,1,0.01))

prediction = gpHistPredict(gp_hist,matrix( X), x_pred)

vars = gpHistVariance(gp_hist,matrix( X), x_pred)

plot(X,Y)
lines(x_pred, prediction,col='red')

lines(x_pred, prediction+sqrt(vars),lty=2,col='red')
lines(x_pred, prediction-sqrt(vars),lty=2,col='red')

legend('bottomleft',legend=c('Data', 'Approximation','Coarse std. dev'), col=c('black','red','red') ,lty=c(NA,1,2),pch=c(1,NA,NA))


```

![Result](https://github.com/dennisthemenace2/gpHist/tree/master/vignettes/gpPlot.png)

