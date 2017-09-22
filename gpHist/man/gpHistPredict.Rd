\name{gpHistPredict}
\alias{gpHistPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gpHistPredict Function
}
\description{

This function is used to predict new sample means. For prediction of the variance the function gpHistVariance is used. This has been done intentially, to provide the chance not to predict the variance if it is not required. Also keep in mind, that you need the data that you used for training of GP to make predictions.

}
\usage{
  gpHistPredict(GP, X,x_pred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GP}{
	Gaussian Process object returned by gpHist function
  }
  \item{X}{
    Original data that has been used for the training of the GP process.
  }

  \item{x_pred}{
    New data that is to  be predicted
  }


}
\details{
  X and x_pred need to be matrix in correct format. Each row is one example of D dimensions.
}
\value{
If the function fails or spotts an error in the parameters NAN is returned. Otherwise the predicted sample mean is returned.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
 Package Overview:
\code{\link[gpHist:gpHist-package]{gpHist-Package}}

Function for estimation of the GP: 
\code{\link[gpHist:gpHist]{gpHist}}

Function for prediction of new sample variance:
\code{\link[gpHist:gpHistVariance]{gpHistVariance}}

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

# New data to predict
x_pred = matrix(seq(0,1,0.01))

#Prediction
prediction = gpHistPredict(gp_hist,matrix( X), x_pred)

# Plot results
plot(X,Y)
lines(x_pred, prediction,col='red')

legend('topleft',legend=c('Data', 'Approximation'), col=c('black','red') ,lty=c(NA,1),pch=c(1,NA))

}


