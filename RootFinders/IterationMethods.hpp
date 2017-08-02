#ifndef ITERATION_METHODS_HPP
#define ITERATION_METHODS_HPP


namespace NumericalAnalysis {

  template<typename T>
  class NewtonMethod1d {
  public:
    NewtonMethod1d(const DifferentiableRealFunction1d<T>* f) : max_steps(10), value_tol(.00001), f(f) {}
    void solve(const T& t0,T& root) const;
  protected:
    unsigned long long max_steps;
    T value_tol;
    const DifferentiableRealFunction1d<T> * f; 
  };
    

  template<typename T>
  class NewtonMethod {
  public:
    NewtonMethod(dim_type dim,const AbstractVectorField<T>* vf,const AbstractDerivativeOfVectorField<T>* dvf_inverse) : max_steps(10), value_tol(.00001), dim(dim), vf(vf), dvf_inverse(dvf_inverse) {}
    void solve(const Point<T>& x0,Point<T>& root) const;
  protected:
    unsigned long long max_steps;
    T value_tol;
    dim_type dim;
    const AbstractVectorField<T> * const vf; 
    const AbstractDerivativeOfVectorField<T> * const dvf_inverse;
  };
 













}




#include "IterationMethods.cpp"

#endif
