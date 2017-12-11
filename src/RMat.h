#ifndef _RMAT_H
#define _RMAT_H

#define COMPILE_WITH_R 1

class RMat {
public:      
    //- Constructors and destructor.
    //
    RMat();                                // Constructor #1
    RMat( int x,   int y);                    // Constructor #2 for x-by-y matrix.
    RMat(double *vElements,   int x,  int y); // Constructor #3 for R interface
    RMat(const RMat &vRhs);                // Constructor #4 (copy constructor)
    ~RMat();                               // Destructor

    //- Overloaded operators.
    //
    RMat& operator= (const RMat &vRhs);     // Assignment #1
    RMat& operator= (double vValue);        // Assignment #2
    RMat operator* (const RMat &vRhs);      // Multiplication
    RMat operator* (const double vRhs);      // Multiplication
    RMat operator/ (const double vRhs);      // Div
    RMat operator- (const RMat &vRhs);      // Subtraction
    RMat& operator-= (const RMat &vRhs);      // Subtraction
    
    RMat operator+ (const RMat &vRhs);      // Add#
    RMat& operator+= (const RMat &vRhs);      // add in place
    RMat& operator+= (const double vRhs);      // add in place
    RMat& operator*= (const double);      // multi in place
    RMat& operator/= (const double);      // multi in place
    
    RMat RowSums();
    RMat ColSums();
    void SetZero();
    double* getColPtr( int col );
    double* getValuesPtr(){return mValues;};
    void setValuesPtr(double* ptr){mValues=ptr;};
    RMat tMultiply(RMat &vRhs);
      
    friend RMat operator/(double vLhs, const RMat &vRhs){
      RMat vProduct(vRhs.mNumRows,vRhs.mNumCols);
      for ( int i=0; i<vRhs.mNumRows*vRhs.mNumCols; i++){
        vProduct.mValues[i] = vLhs/vRhs.mValues[i];
      }
      return vProduct;
    }
    friend RMat operator*(double vLhs, const RMat &vRhs){
      RMat vProduct(vRhs.mNumRows,vRhs.mNumCols);
      for ( int i=0; i<vRhs.mNumRows*vRhs.mNumCols; i++){
        vProduct.mValues[i] = vLhs*vRhs.mValues[i];
      }
      return vProduct;
    }
    
        
   // operator double() const {return mValues[0]; }
    
    double& operator()( int x,  int y) const; // Element access

    double SquareSum();
    double Sum();
    double ScalarProd(const RMat &vRhs);
      
    //- Some utility functions.
    //
    RMat Transpose();                 // Transpose
     int NumRows() const;              // Returns number of rows.
     int NumCols() const;              // Returns number of columns.
  //  void Print();                     // Print matrix to stdout
//    void Print(const char *vString);  // Print matrix to stdout with string.
#if COMPILE_WITH_R
    void RPrint();                    // Print matrix to R console.
    void RPrint(const char *vString); // Print matrix to R console with string.
#endif

    //- A function for deallocating memory is made publicly available because
    //  I want to not only run test this code in a pure C++ test
    //  program outside of R, testRMat.cc, but I also want to run this code
    //  from within R.  In the latter case, I will let R handle memory
    //  management.  This means that I can't automatically deallocate memory
    //  in the destructor, which is the usual thing one might do.
    //  For symmetry, the function for memory allocation is also made publicly
    //  available, although it didn't have to be.  (Maybe it should be private!)
    //
    bool AllocateMemory( int x,  int y); // Allocate memory for x-by-y matrix.
    void DeallocateMemory();           // Deallocate memory

private:
    double  *mValues;           // Matrix values
    int     mNumRows, mNumCols; // Matrix dimensions
    bool    allocated;
    double  mNaN;               // Return NaN on bad input
};

#endif // _RMAT_H defined
