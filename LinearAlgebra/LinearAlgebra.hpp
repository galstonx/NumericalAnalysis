#ifndef __LINEAR_ALGEBRA_HPP
#define __LINEAR_ALGEBRA_HPP

#include <vector>
#include "../Containers/Point.hpp"


namespace NumericalAnalysis {

  /*
    
    ABSTRACT MATRIX CLASS
    
  */
  template<typename T>
  class AbstractMatrix {
  private:
  protected:
  public:
    virtual unsigned long long rows() const=0;
    virtual unsigned long long cols() const=0;
    virtual T get(unsigned long long,unsigned long long) const=0;
    virtual int set(unsigned long long,unsigned long long,T)=0;
  };
  


  /*

    MATRIX CLASS
    T is the number type

  */
  template<typename T>
  class Matrix : public AbstractMatrix<T> {
  private:
    unsigned long long n_rows;
    unsigned long long n_cols;
    std::vector<std::vector<T> > data;
  protected:
  public:
    Matrix(unsigned long long,unsigned long long);
    unsigned long long rows() const;
    unsigned long long cols() const;
    T get(unsigned long long,unsigned long long) const;
    int set(unsigned long long,unsigned long long,T);
    int rowReduce();
    int swapRows(unsigned long long,unsigned long long);
    int rowOperation(T,unsigned long long,unsigned long long);
    int GaussBackSolve(Point<T>&); // only call after row reduce has been called
    template<typename S>
    friend std::ostream& operator<<(std::ostream&,const Matrix<S>&);
  };
    
}


#include "LinearAlgebra.cpp"




#endif
