\name{downhillsimplex}
\alias{downhillsimplex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  downhillsimplex function
}
\description{
	Impementation of the downhillsimplex function for optimisation
}
\usage{

downhillsimplex(fc,nParams,bp=NULL,lower=0,upper=1,it=1000,tol=1e-10,
control= list(alpha=1,gamma=2, beta=0.5,sigma = 0.5) )
  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fc}{
    Function to be minimized
  }
  \item{nParams}{
    Number of parameters that the function to be optimized requires
  }
  \item{bp}{
	Matrix used to start the optimisation. If no parameters is given, this matrix will be randamly initialized.
  }

  \item{lower}{
	Vector for the lower bound of the parameters	
 }

 \item{upper}{
	Vector for the upper bound of the parameters	
 }

 \item{it}{
	Maximium number of iterations done by the downhillsimplex algortihm	
 }

 \item{tol}{
	tolearnce that is used to test for convergence	
 }

 \item{control}{
	list of variables that control the behavior of the downhillsimpley algorithm.	
	alpha = Reflexion
  	gamma = Expansion
  	beta = Contraktion
  	sigma = Compression
 }


}
\details{
	Implemented in C, will call you R function

}
\value{
Returns estimated parameters
}


\references{
Implementation according to describtion on wiki:
https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
}
\author{
    Dennis Becker
}
\note{
 Implementation according to Wikipedia page
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
Package Overview:
\code{\link[gpHist:gpHist-package]{gpHist-Package}}

Function for hyperparameter estimation:
\code{\link[gpHist:estimateHyperParameters]{estimateHyperParameters}}

  }
\examples{
#simple function to be optimized
fn = function(x){
 x^2 
}


#call downhillsimplex
res = downhillsimplex(fn,1,lower=-10,upper=10,it=1000,tol=1e-10)

#plot results
x = seq(-10,10,0.01)
plot(x,fn(x),type='l')
points(res[1],fn(res[1]),col='red')


legend('topleft',legend=c('Function', 'Estimated minimum'), col=c('black','red'),
lty=c(1,NA),pch=c(NA,1))


}
