#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include <vector>

namespace NumericalAnalysis {


  template<typename T>
  class PoissonBasicMethod2d {
  public:
    PoissonBasicMethod2d(const RealFunction<T>* const rf_nonhomog,const RealFunction<T>* const rf_bc,T xmin,T xmax,T ymin,T ymax,unsigned long long x_steps,unsigned long long y_steps) :
      rf_nonhomog(rf_nonhomog), rf_bc(rf_bc), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), n_x_pts(x_steps-1), n_y_pts(y_steps-1),  
									  x_stepsize((xmax-xmin)/(n_x_pts+1)), 
									  y_stepsize((ymax-ymin)/(n_y_pts+1)),
									  p_matrix(n_x_pts*n_y_pts,n_x_pts*n_y_pts+1),
									  soltn_without_bc(n_x_pts*n_y_pts) {
									    init_matrix();
									    init_nonhomog_vector();
									  }
    virtual void solve(Point<T>&);
    virtual void printMatrix();
  protected:
    virtual int init_matrix();
    virtual int init_nonhomog_vector();
    const RealFunction<T>* rf_nonhomog;
    const RealFunction<T>* rf_bc;
    T xmin;
    T xmax;
    T ymin;
    T ymax;
    unsigned long long n_x_pts; // does not include xmin,xmax
    unsigned long long n_y_pts; // does not include ymin,ymax
    T x_stepsize;
    T y_stepsize;
    Matrix<T> p_matrix;
    Point<T> soltn_without_bc;
  };




}






#include "PoissonSolvers.cpp"
#endif
