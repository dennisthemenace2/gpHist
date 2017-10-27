\name{estimateHyperParameters}
\alias{estimateHyperParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  estimateHyperParameters Function
}
\description{
 This function can be used to estimate the hyperparameters of the GP Hist. THe function provides to set bounds on the parameters and the specification of a data transformation function. The data transformation function can be used to make use of the generalized HIK. Where this function provided the g() function that transforms the data. 
}
\usage{
estimateHyperParameters(X, Y, paramLower=0,paramUpper=1,datatransform=NULL,
nParams=0,it = 100,tol=0.001,k=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	 Data used to estimate the GP
  }
  \item{Y}{
	Y values of the data
  }

  \item{paramLower}{
	 lower bound for each parameter
  }

  \item{paramUpper}{
	 upper bound for each parameter
  }
  \item{datatransform}{
	Function that is used to transform the data. This function can make use of hyperparameters, but you have to preserve the orders of the values, because the orders object will be reused to save time
  }

  \item{nParams}{
	Number of hyperparameters to estimate. This will be number of parameters that are used in the datatransform function.
  }

  \item{it}{
	Maximal number of iterations in the downhillsimplex algorithm
  }

  \item{tol}{
	Tollerance till convergence of the downhillsimplex algorithm
  }


  \item{k}{
	Number of estimatied hyperparameters to approximate the log-likelihood. Note that number of eigenvectors scales quadratically in terms of runtime.
  }


}
\details{
 This function is a wrapper for the downhillsimplex algoirthm and can be used to estimate hyperparameters of the GP HIST.
}
\value{
Function will return a matrix of estimated hypervalues.
First item will be the estimated sigma value. The other paramters will be the on you used in the datatransform function.
}


%%\references{
  %% ~put references to the literature/web site here ~
%%}
\author{
  Dennis Becker
}
\note{
	Provides you a wrapper for the downhillsimplex method for hyperparameter estimation  
}

\examples{

#We apply a shift to the function
testFn = function(x){
  y = sin(2*pi*x*2+0.2) 
}

#Get data
X = seq(0,1,0.1)
Y = testFn(X)

#Call gpHist function
#gp_hist = gpHist(matrix(X),matrix(Y),sigma=0.01)


#add a positive shift to the data
transform= function(X,p){
  X+p
}

##estimate hyperparameters
res = estimateHyperParameters(matrix(X),matrix(Y),paramLower = c(0.001,0.0001),
paramUpper = c(0.01,1),datatransform = transform,nParams = 1)

##data need to be trannsformed for the training
X_trans = matrix(transform(X,res[2] ) )
gp_hist = gpHist(X_trans,matrix(Y),sigma=res[1])

##new data to be predicted has also be to transformed
x_pred = matrix(seq(0,1,0.01))
x_pred_trans = transform(x_pred,res[2])

prediction = gpHistPredict(gp_hist,X_trans, x_pred_trans)

plot(X,Y)
## note that we need to use the data before transformation for real x values
lines(x_pred, prediction,col='red')


legend('topleft',legend=c('Data', 'Approximation'), col=c('black','red'),
lty=c(NA,1),pch=c(1,NA))

}
%% ~Make other sections like Warning with \section{Warning }{....} ~
  
\seealso{
	 Package Overview:
	\code{\link[gpHist:gpHist-package]{gpHist-Package}}

Implementaion of the downhillsimplex algortihm:
\code{\link[gpHist:downhillsimplex]{downhillsimplex}}
  }

