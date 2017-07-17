#ifndef __POISSON_SOLVER_HPP
#define __POISSON_SOLVER_HPP

#include <vector>

namespace NumericalAnalysis {


  template<typename T>
  class PoissonBasicMethod2d {
  private:
  protected:
    Matrix<T> p_matrix;
    virtual int init_matrix();
    virtual int init_nonhomog_vector();
    T xmin;
    T xmax;
    T ymin;
    T ymax;
    T n_x_pts; // does not include xmin,xmax
    T n_y_pts; // does not include ymin,ymax
    T x_stepsize;
    T y_stepsize;
    const std::vector<T>* ic_xmin;
    const std::vector<T>* ic_xmax;
    const std::vector<T>* ic_ymin;
    const std::vector<T>* ic_ymax;
    //std::vector<T> bc;
    const RealFunction<T>* rf_nonhomog;
    const RealFunction<T>* rf_bc;
    Point<T> soltn_without_bc;
  public:
    PoissonBasicMethod2d(const RealFunction<T>*,const RealFunction<T>*,T,T,T,T,unsigned long long,unsigned long long);
    virtual int solve(Point<T>&);
    virtual int printMatrix();
  };




}






#include "PoissonSolvers.cpp"
#endif
