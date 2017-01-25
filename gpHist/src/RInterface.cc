
#include "RMat.h"
#include "R.h"     // R functions
#include "Rmath.h" // Rmath

#include <R_ext/Utils.h>
#include <Rinternals.h>


#include <iostream> // For cout, cerr, endl, and flush

void fastMultiply(const RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){
  
//  RMat prev(Y.NumRows ,1);
  result.SetZero();
  unsigned int NR = X.NumRows();
    
  for(unsigned int i=1;i<=X.NumCols();++i){
  
    unsigned int idx = 0;
    double matcount = 0;
    double Ysum;
    
    for(unsigned int k=1;k<NR;++k ){ // estimate prev
      idx = orders(k,i);
      
      Ysum= 0;
      for(unsigned int x=k;x<=NR;++x ){ // estimate Y.sum(k,NR); 
        Ysum += Y(orders(x,i) ,1);
      }
      
      result(idx,1) +=  matcount + X(idx, i)*Ysum;
      matcount = matcount + X(idx,i) * Y(idx,1);

    }
    idx = orders(NR,i);
    
    result(idx,1) += matcount + X(idx,i) * Y(idx,1);
      
  }
 // for(unsigned i =1;i<=Y.NumRows();++i){
  //  result(i,1)  += sigma*Y(i,1);
//  }
  result += sigma*Y ;
  
}



bool conjgradFP(const RMat& A, RMat& b,RMat& x ,RMat& orders,double sigma, unsigned int max_iterations=1000 ) {

  RMat res( x.NumRows(),1 );

  fastMultiply(A,x,res,orders,sigma); //A%*%x;

  RMat r= b-res;// #A%*%x;
  RMat p=r;
  
  double rsold=r.SquareSum();//t(r)%*%r
    
  for (unsigned int i=1;i<=max_iterations;++i){
    fastMultiply(A,p,res,orders,sigma);// res =  A%*%p;
//double alpha= rsold / ( p.Transpose() *res);
//  std::cout<<alpha<<std::endl;
    double alpha= rsold / p.ScalarProd(res) ;

    x += alpha*p;
    
    if(i % 50 == 0){
      fastMultiply(A,x,res,orders,sigma);
      r=b-res; //#A%*%x;
  //    std::cout<<"update"<<std::endl;
  //    r.Print();
    }else{
      r -= alpha*res;
     // r.Print();
    }
    
    
    double rsnew= r.SquareSum();//t(r)%*%r;
    if(sqrt(rsnew) <1e-12){
//#  cat(c('break:', i, ' err:',sum(sqrt(rsnew)),'\n') )
      return true;
    }
    p*=(rsnew/rsold);
    p+=r;
    //p.Print();
    rsold=rsnew;
  }
    
  return false;
}
// power method
double getEigen(const RMat& A,RMat &b,RMat& orders,double sigma){
  double last = 0;
  RMat res( b.NumRows(),1 );
  
  for(unsigned int i=0;i<10000;++i){
    fastMultiply(A,b,res,orders,sigma); //A%*%b;
    double norm = sqrt(res.SquareSum()) ;
  //  res.Print();
   // std::cout<< norm<<std::endl;
    b = res;
    b /= norm;
  //  b.Print();
    if(abs(last-norm)<1e-10 ){
        break;
    }
    last= norm;
  }
  return last;
}

// two helper functions for predictions
double sumOrdered(RMat& alpha,RMat&orders,int start,int end,int col){
  double result = 0;
  for(int i= start;i<=end;++i){
    int idx = orders(i,col);
    result += alpha(idx,1);
  }
  return result;
}
  
double sumSquareOrdered(RMat& X,RMat&orders,int start,int end, int col){
  double result = 0;
  for(int i= start;i<=end;++i){
    int idx = orders(i,col);
    //  std::cout<<"sumordidx"<<idx<<std::endl;
    result += X(idx,col)*X(idx,col);
  }
  
  return result;
}
    
double sumOrdered(RMat& X,RMat& alpha,RMat&orders,int start,int end, int col){
  double result = 0;
  for(int i= start;i<=end;++i){
    int idx = orders(i,col);
  //  std::cout<<"sumordidx"<<idx<<std::endl;
    result += X(idx,col)*alpha(idx,1);
  }
  
  return result;
}

int findIdx(RMat& X,RMat& orders,double value,int col){
  int k =1;
  for(;k<=X.NumRows();++k){ // this is linear this suckz
    int idx = orders(k,col);
    if(value<X(idx,col)){
      break;
    }
  }
  
  return k-1;
}
void CppHistVariance(double* result,double numRows,double numCols,double numRows2,double numCols2,double* mat1,double* mat2,double *mat3,double lambda,double sigma,double* orders){
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatPred(mat2,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
  RMat  vMatAlpha(mat3,numRows,1);
  RMat  vMatOrders(orders,numRows,numCols);
  
  for(int i=1; i<=numRows2;++i ){
      double A=0;
      double B=0;
      
    for(int d=1; d<= numCols2;++d){
      int idx = findIdx(vMatX,vMatOrders, vMatPred(i,d),d);
     // std::cout<<"idx="<<idx<<std::endl;
      
       if(idx>0){
           A = sumSquareOrdered(vMatX,vMatOrders,1,idx,d);
      //   std::cout<<"A="<<A<<std::endl;
        }else{
          A=0;
          idx=0;
        }
        
        B = (numRows-idx )* ( vMatPred(i,d)*vMatPred(i,d) );
        //std::cout<<"B="<<B<<std::endl;
     // C = C+ (A+B);
    }
    vMatResult(i,1) = A+B;
  }
//  vMatResult.Print("vmatresult");
  vMatResult *= (1/lambda);

  
 // vMatPred.RowSums().Print("rowssums") ;
  
  vMatResult = vMatPred.RowSums() -vMatResult;
  vMatResult += sigma;

}

void CppHistPredict(double *result,
                    int numRows, int numCols,int numRows2 ,int numCols2,
                    double *mat1, double *mat2,double* mat3,double *orders){
  
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatPred(mat2,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
  RMat  vMatAlpha(mat3,numRows,1);
  RMat  vMatOrders(orders,numRows,numCols);
  
  //NRA = nrow(GP$alpha)
    
  //  pred = c()
    for(int i=1; i<=numRows2;++i ){
      double C=0;
      double A=0;
      double B=0;
      
      for(int d=1; d<= numCols2;++d){
        // find idx....
        
        int idx = findIdx(vMatX,vMatOrders, vMatPred(i,d),d);
      //  std::cout<<"idx"<<idx<<std::endl;
        
        if(idx>0){
          if(idx==numRows){//###largest
            C = C+ sumOrdered(vMatX,vMatAlpha,vMatOrders,1,numRows,d);
        //    std::cout<<"C"<<C<<std::endl;
            continue;
          }else{
            A = sumOrdered(vMatX,vMatAlpha,vMatOrders,1,idx,d);
      //      std::cout<<"A"<<A<<std::endl;
          }
        }else{ //##point even further left....
          A=0;
        }
        B = sumOrdered(vMatAlpha,vMatOrders,idx+1,numRows,d) * vMatPred(i,d);
    //    std::cout<<"B"<<B<<std::endl;
        C = C+A+B;
      }
     
     vMatResult(i,1)=C;
    }
    
  
}

void CppHist(double *result,
                   int numRows, int numCols,int numRows2 ,int numCols2,
                   double *mat1, double *mat2,double sigma , double *orders, double* logmarginal,double* lambda)
{
    RMat  vMatX(mat1,numRows,numCols);
    RMat  vMatY(mat2,numRows2,numCols2);
    RMat  vMatAlphas(result,numRows,1);
    RMat  vMatOrders(orders,numRows,numCols);
    
   
   if(!conjgradFP(vMatX, vMatY ,vMatAlphas,vMatOrders, sigma ) ){
     Rprintf("conjgradFP not converged\n");
   }
   //vMatAlphas.Print();
   
   double mu1 = numRows*sigma + vMatX.Sum();
   // nedd to get larges eigenvalue
   RMat eigenVec(numRows,1);
   eigenVec = vMatY;//1.0/(*numRows);
   
   double beta = getEigen(vMatX,eigenVec,vMatOrders, sigma);
   *lambda = beta;
//   std::cout<<beta<<std::endl;
 //  std::cout<<" "<<std::endl;
 //  eigenVec.Print();
   
   double mu2 =  beta *beta;
     
   double t_bar = (beta*mu1 - mu2) / (beta*numRows - mu1);
   double t_bar2 = t_bar*t_bar;
   double beta2 = beta*beta; /// CHECK THIS 
     
   RMat t1(1,2);
   t1(1,1) =log(beta);
   t1(1,2) =log(t_bar);
 
   RMat t2(2,2);
   t2(1,1) = t_bar2;
   t2(2,1) = -beta2;
   t2(1,2) = -t_bar;
   t2(2,2) = beta;  
   t2*=(1/(beta*t_bar2 - t_bar*beta2 ));
       
   RMat t3(2,1);
   t3(1,1)=mu1;
   t3(2,1)=mu2;
    
   double logdetM =  (t1 * t2).ScalarProd(t3);
   
 //  std::cout<<"logdet:"<<logdetM<<std::endl;
   
   *logmarginal= 0.5*vMatY.ScalarProd( vMatAlphas ) + logdetM/2 + (numRows/2.0)*log(2.0*PI);
//   std::cout<<"   vMatY.ScalarProd( vMatAlphas ) :"<<   vMatY.ScalarProd( vMatAlphas ) <<std::endl;
//   std::cout<<"    numRows/2*log(2*PI):"<<   (numRows/2.0)*log(2.0*PI)<<std::endl;
   
  // std::cout<<"approxmarginal:"<<*logmarginal<<std::endl;


}

void Cmtest(double *result,
            int *numRows, int *numCols,int *numRows2 ,int *numCols2,
            double *mat1, double *mat2,double *sigma , double *orders)
{
  RMat  vMatX(mat1,*numRows,*numCols);
  RMat  vMatY(mat2,*numRows2,*numCols2);
  RMat  vMatAlphas(result,*numRows,1);
  RMat  vMatOrders(orders,*numRows,*numCols);
  
  fastMultiply(vMatX, vMatY ,vMatAlphas,vMatOrders, *sigma );
  
}


extern "C" {
    void CgpHist(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *sigma,double *orders,double* logmarginal,double *lambda){
      CppHist(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,*sigma,orders,logmarginal,lambda);
    }

  void CgpHistPredict(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *mat3,double *orders){
    CppHistPredict(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,mat3,orders);
  }
  
  void CgpHistVariance(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *mat3,double*lambda,double* sigma, double *orders){
    CppHistVariance(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,mat3,*lambda,*sigma,orders);
  }

  void Cdownhillsimplex(SEXP func,SEXP sys,int* nParams,double* bp){
    SEXP arg;
    SEXP environment;
    
    PROTECT(environment = sys);
    PROTECT(arg = allocVector(REALSXP, *nParams) );
 // memcpy(REAL(arg), 5, (*nParams) * sizeof(double));
    double *input = REAL(arg);
 
    for( int i = 0;i<*nParams;++i){
      input[i] =i;
    }
    
    SEXP function_call;
    PROTECT(function_call = lang2(func, arg));
    
    SEXP s_fval;
    PROTECT(s_fval = eval(function_call, environment));
    
    if(length(s_fval) != 1) {
      error("Function should return error to minimize");
    }
    double *err = REAL(s_fval);
    std::cout<<"res:"<< *err <<std::endl;
    
    UNPROTECT(4); /* s_fval */
        
        
    //Cdownhillsimplex SEXP
    std::cout<<"OK"<<std::endl;
  }
  
  
  void CgpfMult(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *sigma,double *orders)
  {
    Cmtest(result,numRows,numCols,numRows2,numCols2,mat1,mat2,sigma,orders);  
  }
  
}
