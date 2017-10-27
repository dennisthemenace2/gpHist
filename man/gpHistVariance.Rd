\name{gpHistVariance}
\alias{gpHistVariance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  gpHistVariance Function
}
\description{
  This function is used to estimate (approximate) the variance of new sample from a trained GP.
}
\usage{
gpHistVariance( GP,X,x_pred)
}

\arguments{
  \item{GP}{
    Gaussian Process object that has been trained with gpHist
  }
  \item{X}{
    Data that has been used to train the GP object
  }
  \item{x_pred}{
    New data whos variance has to be predicted
  }
}
\details{
Based on the number of estimaited eigenvalues/eigenvetors, the variance is approximated. When only 1 eigenvalue is used, the so called coarse approximation is used. If more eigenvectors are present, the so called fine approximation is used. For deatils look in to paper of refrenced under references
}
\value{
 Will return NAN on error or the estimated approximaited variances of the samples

}


\references{
The utilized approximations are described in the following paper:

Rodner, E., Freytag, A., Bodesheim, P., Froehlich, B., & Denzler, J. (2016). Large-Scale Gaussian Process Inference with Generalized Histogram Intersection Kernels for Visual Recognition Tasks. International Journal of Computer Vision, pp. 1-28. Springer US. doi:10.1007/s11263-016-0929-y

}
\author{
  Dennis Becker
}
\note{
  A negative variance can be retunred. This can happen when one eigenvalue is not the highest, rather a lower one. Might also happen when eigenvalue and/or eigenvector is not estimated accurate enough. 
As workaround consider using less eigenvectors for the estimation of the Gaussian process. Or remove the lowest eigenvalue and vector from the object and call the variance function again. This of course will you provide a worse approximation, but give you a valid variance that you might be able to use instead.
}

\section{Warning }{ A negative variance cen be retunred. For reasons and workarounds look into note section.}


  \seealso{

 Package Overview:
\code{\link[gpHist:gpHist-package]{gpHist-Package}}

Function for estimation of the GP: 
\code{\link[gpHist:gpHist]{gpHist}}

Function for prediction of the mean of new samples:
\code{\link[gpHist:gpHistPredict]{gpHistPredict}}

Function for hyperparameter estimation:
\code{\link[gpHist:estimateHyperParameters]{estimateHyperParameters}}

  }
\examples{

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

legend('topleft',legend=c('Data', 'Approximation','Coarse std. dev'),
col=c('black','red','red') ,lty=c(NA,1,2),pch=c(1,NA,NA))

}
