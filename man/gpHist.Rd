\name{gpHist}
\alias{gpHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gpHist function
}
\description{

This function is used to estimate the Gaussian Process regression
}
\usage{
gpHist(X,Y,sigma,orders=NULL,alphas = NULL,k=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	Matrix of X values. Make sure that this is a matrix in the right dimension. It should be a NxD matrix. So that each row is one example of D dimensions.
}
  \item{Y}{
	Matrix of Y values. Make sure that this is a matrix in the right dimension. Should be a 1xN matrix.
}

\item{sigma}{
	Noise in the data.	
}

\item{orders}{
	If the correct orders of the samples in each dimension is already known, this can be used to avoid sorting of each dimension. If not set, order will be created and returned.
}

 \item{alphas}{
 With this parameters the initial value of the alphas can be set. If not set, a random vector will be used. By setting the alpha vector the convergence can potentially be speed up.
}
 \item{k}{
 The number of estimated eigenvalues and eigenvectors. This effects the approximation of the log-likelyhood and the predicted variance. It has no influence on the prediction of the mean of new sample values.
Estimation of eigenvectors scales quadratically, so be awere of this additional runtime.
}
}
\details{
	This function will estimate the Gaussian process with HIK kernen using conjugate gradient descent and return all required estimates as lists.
If only one eigenvector is estimated the power method will be used. For multible eigenvectors, Lanczos algorihm is applied. The eigenvalues will be bounded using sturm sequences and then estimatied using inverse iteration, wich also makes use of the conjugate gradient descent. Estimated time of multible eigenvectors greatly varies.
}
\value{
If the functions fails or wrong parameters are provided NAN is returned.


%% list('logmarginal'=output$logmarginal,'orders'=orders,'alpha'=matrix(output$result,nrow=nrowA),'sigma'=sigma,'lambda'=output$lambda,'vectors'=matrix(output$vector,nrow=nrow(X), ncol=k))

Return value of the function is a List with following items:
\item{logmarginal }{Approximated log-likelihood of the GP}
\item{orders }{Matrix that describes the ordering of each dimension}
\item{alpha }{Estimatied alpha vector. alpha = K_xx^-1 y}
\item{lambda }{Esimated eigenvalues}
\item{vectors }{Estimated eigenvectors}

}

\author{
Dennis Becker
}
\note{
You need to keep the X data that you used for training because it is not saved within the returned GP object.
 In the future more parameters will be added to controll the internals of the function.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Package Overview:
\code{\link[gpHist:gpHist-package]{gpHist-Package}}

Function for prediction of new samples: 
\code{\link[gpHist:gpHistPredict]{gpHistPredict}}

Function for prediction of new sample variance:
\code{\link[gpHist:gpHistVariance]{gpHistVariance}}

Function for hyperparameter estimation:
\code{\link[gpHist:estimateHyperParameters]{estimateHyperParameters}}

}
\examples{

require('gpHist')


testFn = function(x){
  y = sin(2*pi*x*2) 
}

#Get data
X = seq(0,1,0.01)
Y = testFn(X)

#Call gpHist function
gp_hist = gpHist(matrix(X),matrix(Y),sigma=0.01)
}

