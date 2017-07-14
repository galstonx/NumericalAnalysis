#ifndef __BOUNDARY_VALUE_METHODS_HPP
#define __BOUNDARY_VALUE_METHODS_HPP

#include "../Containers/Functions.hpp"
#include "../RootFinders/EnclosureMethods.hpp"

namespace NumericalAnalysis {

  template<typename T>
  class RealFunction1dForShooting : public RealFunction1d<T> {
  private:
  protected:
    const VF2dFrom2ndOrder<T> vf; // first order vf determined by F
    T t_start;
    T t_stop;
    T x_start;
    T x_stop;
    unsigned long long steps;
    ODEMethod<T>* ode_method;
  public:
    RealFunction1dForShooting(RealFunctionNonAut<T> const *,unsigned long long,T,T,T,T);
    virtual void eval(T,T&) const; // given initial velocity, find position at t_stop minus x_stop
    ~RealFunction1dForShooting();
  };


  template<typename T>
  class ShootingMethodNonlinear1d {
  private:
  protected:
    RealFunctionNonAut<T> const * F; // RHS of 2nd order ODE
    unsigned long long steps_ode_method;
    T t_start;
    T t_stop;
    T x_start;
    T x_stop;
    const RealFunction1dForShooting<T> rf;
  public:
    ShootingMethodNonlinear1d(const RealFunctionNonAut<T>*,unsigned long long,T,T,T,T);
    // get the value of the function with a certain velocity
    int getFunctionValue(T,T&) const;
    // Dirichlet using Bisection method
    int solveDirichletBisection(T,T,T&) const;
  };


  template<typename T>
  class VF2dForEigenvalue : public VectorFieldNonAut<T> {
  private:
  protected:
    RealFunctionNonAut<T> const * F; // RHS of 2nd order ODE
    T lambda;
  public:
    VF2dForEigenvalue(RealFunctionNonAut<T> const *);
    int setEigenvalue(T);
    virtual void eval(T,const Point<T>&,Point<T>&) const;
  };


  template<typename T>
  class RealFunction1dForEigenvalue : public RealFunction1d<T> {
  private:
  protected:
    VF2dForEigenvalue<T> vf; // first order vf determined by F
    T t_start;
    T t_stop;
    unsigned long long steps;
    ODEMethod<T>* ode_method;
    T lambda;
    T vel_def; // default velocity to use for checking
    Point<T> ic;
  public:
    RealFunction1dForEigenvalue(RealFunctionNonAut<T> const *,unsigned long long,T,T);
    int setEigenvalue(T);
    virtual void eval(T&) const; // given lambda, find position at t_stop
    ~RealFunction1dForEigenvalue();
  };


  template<typename T>
  class EigenvalueFor2ndOrder {
  private:
  protected: 
    RealFunctionNonAut<T> const * F; // RHS of 2nd order ODE
    unsigned long long steps_ode_method;
    T t_start;
    T t_stop;
    RealFunction1dForEigenvalue<T> rf;
    T value_tol;
  public:
    EigenvalueFor2ndOrder(const RealFunctionNonAut<T>*,unsigned long long,T,T);
    int getEigenvalues(unsigned long long,T,T,std::vector<T>&);
    


  };





}


#include "BoundaryValueMethods.cpp"



#endif
