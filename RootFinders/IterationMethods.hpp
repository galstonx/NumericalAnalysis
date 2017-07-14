#ifndef __ITERATION_METHODS_HPP
#define __ITERATION_METHODS_HPP


namespace NumericalAnalysis {

  template<typename T>
  class NewtonMethod1d {
  private:
    unsigned long long MAX_STEPS_DEF;
    double VALUE_TOL_DEF;
  protected:
    unsigned long long max_steps;
    T value_tol;
    RealFunction1d<T> const * rf; 
  public:
    NewtonMethod1d(const RealFunction1d<T>&);
    int solve(T,T&);
  };
    

  template<typename T>
  class NewtonMethod {
  private:
    unsigned long long MAX_STEPS_DEF;
    double VALUE_TOL_DEF;
    unsigned long long dim;
  protected:
    unsigned long long max_steps;
    T value_tol;
    VectorField<T> const * vf; 
  public:
    NewtonMethod(unsigned long long,const VectorField<T>&);
    int solve(Point<T>,Point<T>&);
  };
 













}




#include "IterationMethods.cpp"

#endif
