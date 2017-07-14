namespace NumericalAnalysis {


#include<iostream>


  template<typename T>
  int GaussBackSolve(Matrix<T>& aug_matrix,Point<T>& ans) {
    int rv=0;
    aug_matrix.rowReduce();
    rv=aug_matrix.GaussBackSolve(ans);
    return rv;
  }












}
