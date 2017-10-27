#include "RMat.h"   // For RMat class
#include <limits>   // For NaN
#include <iostream> // For cout, cerr, endl, and flush
#include <assert.h> // For assert

#if COMPILE_WITH_R
#include "R.h"      // For Rprintf
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::flush;

//- Constructor #1
//
RMat::RMat() : mValues(0), mNumRows(0),
mNumCols(0), mNaN(std::numeric_limits<double>::quiet_NaN()),allocated(false)
{
}

//- Constructor #2
RMat::RMat(int x, int y) : mNumRows(x),mNumCols(y),allocated(false),mValues(0),
 mNaN(std::numeric_limits<double>::quiet_NaN())
{
    AllocateMemory(x,y);
}

//- Constructor #3
//  This constructor will be used to point an RMat object to the memory
//  location provided by R, via the '.C()' function.
//
RMat::RMat(double *vElements, int x, int y) :
mNaN(std::numeric_limits<double>::quiet_NaN()),allocated(false)
{
    mValues  = vElements;
    mNumRows = x;
    mNumCols = y;
}

//- Constructor #4 (copy constructor)
//
RMat::RMat(const RMat &vRhs) : mValues(0), mNumRows(0),
mNumCols(0), mNaN(std::numeric_limits<double>::quiet_NaN())
{
  AllocateMemory(vRhs.NumRows(),vRhs.NumCols());
  
  memcpy(mValues,vRhs.mValues,mNumRows*mNumCols*sizeof(double));
}

//- Destructor
//
RMat::~RMat(){
  DeallocateMemory();
}

void RMat::SetZero(){
  if(mValues!= NULL){
    memset( mValues,0,mNumCols*mNumRows*sizeof(double) );
  }
}
//- Assignment operator #1
//
RMat& RMat::operator= (const RMat &vRhs){
 
    if(mNumCols != vRhs.mNumCols || mNumRows!=vRhs.mNumRows){
      DeallocateMemory();
      AllocateMemory(vRhs.NumRows(),vRhs.NumCols());
    }
    for ( int i=0; i<mNumRows*mNumCols; i++){
      mValues[i] = vRhs.mValues[i];
    }
    return (*this);
}

//- Assignment operator #2
//
RMat& RMat::operator= (double vValue){
    for (int i=0; i<mNumRows*mNumCols; i++){
        mValues[i] = vValue;
    }
    return (*this);
}


RMat& RMat::operator-= (const RMat &vRhs){
  assert(mNumCols==vRhs.mNumCols);
  assert(mNumRows==vRhs.mNumRows);
    
  for (int i=0; i<mNumRows*mNumCols; i++){
    mValues[i] -= vRhs.mValues[i];
  }
  return (*this);
}

//- Subtraction operator
//

RMat RMat::operator- (const RMat &vRhs){
  assert(mNumCols==vRhs.mNumCols);
  assert(mNumRows==vRhs.mNumRows);
  
  RMat vProduct(mNumRows,mNumCols);
  for (unsigned int i=0; i<mNumRows*mNumCols; i++){
    vProduct.mValues[i] = mValues[i] - vRhs.mValues[i];
  }
  return vProduct;
}

RMat RMat::operator/ (const double vRhs){
  RMat vProduct(mNumRows,mNumCols);
  for (unsigned int i=0; i<mNumRows*mNumCols; i++){
    vProduct.mValues[i] = mValues[i] / vRhs;
  }
  return vProduct;
}

//Add
RMat RMat::operator+ (const RMat &vRhs){
    assert(mNumCols==vRhs.mNumCols);
    assert(mNumRows==vRhs.mNumRows);
    
    RMat vProduct(mNumRows,mNumCols);
    for (unsigned int i=0; i<mNumRows*mNumCols; i++){
      vProduct.mValues[i]= mValues[i]+vRhs.mValues[i]; 
    }
    return vProduct;
  }

RMat& RMat::operator+= (const double vRhs){
  for (unsigned int i=0; i<mNumRows*mNumCols; i++){
    mValues[i] +=  vRhs;
  }
    
  return (*this);
}

RMat& RMat::operator*= (const double vRhs){
  for (unsigned int i=0; i<mNumRows*mNumCols; i++){
    mValues[i] *=  vRhs;
  }
  
  return (*this);
}

RMat& RMat::operator/= (const double vRhs){
  for (unsigned int i=0; i<mNumRows*mNumCols; i++){
    mValues[i] /=  vRhs;
  }
    
  return (*this);
}

RMat& RMat::operator+= (const RMat &vRhs){
    assert(mNumCols==vRhs.mNumCols);
    assert(mNumRows==vRhs.mNumRows);
    
    for (int i=0; i<mNumRows*mNumCols; i++){
      mValues[i] += vRhs.mValues[i];
    }
    return (*this);
}

double RMat::Sum(){
    double result=0;
    for ( int i=0; i<mNumRows*mNumCols; ++i){
      result+= mValues[i]; 
    }
    return result;
}


double RMat::ScalarProd(const RMat &vRhs){

  if (mNumCols!= vRhs.mNumCols || mNumRows!= vRhs.mNumRows){
   // cerr << "Rmat::ScalarProd: matrices have not the same dimension!" << endl;
    Rprintf("Rmat::ScalarProd: matrices have not the same dimension!\n");
    return (double &) mNaN;
  }
  double result=0;
  for ( int i=0; i<mNumRows*mNumCols; ++i){
    result+= mValues[i]*vRhs.mValues[i]; 
  }
  
  return result;
}

double RMat::SquareSum(){// or call this square norm, could also call ScalarProd
  double result=0;
  for ( int i=0; i<mNumRows*mNumCols; ++i){
    result+= mValues[i]*mValues[i]; 
  }
  
  return result;
}

// THIS IS a hack
double* RMat::getColPtr( int col ){
  if ( col<1 || col>mNumCols){
    //cerr << "Rmat::getColPtr: range error!" << endl;
    Rprintf("Rmat::getColPtr: range error!\n");
    return NULL;
  }
  return &mValues[ (col-1)*mNumRows];
}

/*
//RMat RMat::operator/(double vLhs, const RMat &vRhs){
  RMat vProduct(vRhs.mNumRows,vRhs.mNumCols);
  for ( int i=0; i<vRhs.mNumRows*vRhs.mNumCols; i++){
    vProduct.mValues[i] = vLhs/vRhs.mValues[i];
  }
  
  return vProduct;
}
*/

//- Multiplication operator
//
RMat
RMat::operator* (const RMat &vRhs) // maybe make this faster..
{
    assert(mNumCols==vRhs.mNumRows);
    RMat vProduct(mNumRows,vRhs.NumCols());
    for ( int i=1; i<=mNumRows; i++)
    {
        for ( int j=1; j<=vRhs.NumCols(); j++)
        {
            for ( int k=1;k<=mNumCols;k++)
            { 
                vProduct(i,j)+=(*this)(i,k)*vRhs(k,j); 
            } 
        }
    }
    return vProduct;
}

RMat RMat::tMultiply (RMat &vRhs){
  
  assert(mNumRows==vRhs.mNumRows);
  
  RMat vProduct(mNumCols,vRhs.NumCols());
  for ( int i=1; i<=mNumCols; ++i){
    double *p1=getColPtr(i);
    for ( int j=1; j<=vRhs.NumCols(); ++j){
      double *p2= vRhs.getColPtr(j);
      double sum=0;
      for(int z=0;z<mNumRows;++z){
        sum+= p1[z]*p2[z];
      }
      vProduct(i,j)=sum;
    }
  }
  return vProduct;
}


RMat RMat::ColSums(){
  RMat vProduct(mNumCols,1);
  for(int j=0;j<mNumCols;++j){
    for(int i=0;i<mNumRows;++i){
      vProduct.mValues[j] += mValues[ (j)*mNumRows+(i)];
    }
  }
  return vProduct;
}

RMat RMat::RowSums(){
  RMat vProduct(mNumRows,1);

  for(int j=0;j<mNumCols;++j){
    for(int i=0;i<mNumRows;++i){
      vProduct.mValues[i] += mValues[ (j)*mNumRows+(i)];
    }
  }
  return vProduct;
}

RMat RMat::operator* (const double vRhs){
    
  RMat vProduct(mNumRows,mNumCols);
  for ( int i=0; i<mNumRows*mNumCols; i++){
       vProduct.mValues[i] =  mValues[i]* vRhs;
  }
  return vProduct;
}

//- Element access
//
double&
RMat::operator()( int x, int y) const
{
    //- Basic range checks.
    //
    if ( ( x < 1 ) || ( x > mNumRows )
      || ( y < 1 ) || ( y > mNumCols ) )
    {
        // cerr << "Rmat::operator(): range error!" << endl;
        Rprintf("Rmat::operator(): range error!");
        return (double &) mNaN;
    }
    
    //- Offset is ((y-1)*mNumRows)+(x-1)
    //  rather than ((x-1)*mNumCols)+(y-1)
    //  because we want to organize memory the same
    //  way as it is in an R matrix.  R is column-major.
    //  Reference:
    //  http://en.wikipedia.org/wiki/Row-major_order#Column-major_order
    //
    return mValues[((y-1)*mNumRows)+(x-1)];
}


//- Transpose
//
RMat RMat::Transpose() // this should be slow...
{
    RMat vTranspose(mNumCols,mNumRows);
    for (int i=1; i<=mNumRows; i++)
    {
        for (int j=1; j<=mNumCols; j++)
        {
            vTranspose(j,i) = (*this)(i,j); 
        }
    }
    return vTranspose;
}


#if COMPILE_WITH_R
//- RPrint #1
//
void RMat::RPrint(){
    for ( int i=1; i<=mNumRows; i++)
    {
        for ( int j=1; j<=mNumCols; j++)
        {
            Rprintf("%g  ",(*this)(i,j));
        }
        Rprintf("\n");
        R_FlushConsole();
        R_ProcessEvents();
    }
}

//- RPrint #2
//
void RMat::RPrint(const char *vString){
    Rprintf("%s",vString);
    RPrint();
}
#endif

//- Get number of rows
//
int RMat::NumRows() const{
    return mNumRows;
}

//- Get number of columns
//
int RMat::NumCols() const {
    return mNumCols;
}

//- Allocate memory
//  Initialize all entries to 0.
//
bool RMat::AllocateMemory( int x, int y)
{
    DeallocateMemory();
    try {
        mValues  = new double [x*y];
    }
    catch (...) {
        /*cerr << "Rmat::AllocateMemory(): unable to allocate memory for "
               << x << " by " << y << " matrix!" << endl;*/
        Rprintf("Rmat::AllocateMemory(): unable to allocate memory\n");
        return false;
    }

    mNumRows = x;
    mNumCols = y;
    allocated=true;
    SetZero(); // do i realy need to zero this ?
  
    return true;
}

//- Deallocate memory
//
void RMat::DeallocateMemory(){
    if ( mValues != 0 && allocated){ // make sure we dont delete R's memory
        delete [] mValues;
        mValues = 0;
    }
}
