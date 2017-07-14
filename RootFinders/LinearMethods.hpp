#ifndef __LINEAR_METHODS_H
#define __LINEAR_METHODS_H

#include "../LinearAlgebra/LinearAlgebra.hpp"
#include "../Containers/Point.hpp"


namespace NumericalAnalysis {


  // assume the first arg is an augmented matrix
  // the second arg is the solution
  // if no soltn, returns 1 (no guaranterr on what 2nd arg is)
  // otherwise returns 0
  // if not a unique solution, no guarantee what the answer is -- it will use whatever
  // is currently in the free entries of the 2nd arg
  template<typename T>
  int GaussBackSolve(Matrix<T>&,Point<T>&);




















}








#include "LinearMethods.cpp"

#endif
