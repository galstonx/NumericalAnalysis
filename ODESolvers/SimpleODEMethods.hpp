#ifndef SIMPLE_ODE_METHODS_HPP
#define SIMPLE_ODE_METHODS_HPP

#include "../Sizes.hpp"
#include "../Containers/Functions.hpp"

namespace NumericalAnalysis {

  template<typename T>
  class ODEMethod {
  public:
    // steps, start time,end time,initial condition,vector for answer
    virtual void solve(unsigned long long,T,T,const Point<T>&,std::vector<Point<T> >&) const=0;
    // steps, start time,end time,initial condition,T& for answer
    virtual void solveFinalVal(unsigned long long,T,T,const Point<T>&,Point<T>&) const=0;
  protected:
    dim_type dim;
  };

  template<typename T>
  class EulerMethod : public ODEMethod<T> {
  public:
    EulerMethod(const AbstractVectorFieldNonAut<T>* vf)  : vf(vf) {ODEMethod<T>::dim=vf->getDim();}
    // steps, start time,end time,initial condition,vector for answer
    virtual void solve(unsigned long long n_steps,T t_start,T t_stop,const Point<T>& x0,std::vector<Point<T> >&) const;
    // steps, start time,end time,initial condition,T& for answer
    virtual void solveFinalVal(unsigned long long n_steps,T t_start,T t_stop,const Point<T>& x0,Point<T>& rv) const;
  protected:
    const AbstractVectorFieldNonAut<T> * const vf;
  };



  template<typename T>
  class RungeKuttaMethod : public ODEMethod<T> {
  public:
    RungeKuttaMethod(const AbstractVectorFieldNonAut<T>* vf) : vf(vf) {ODEMethod<T>::dim=vf->getDim();}
    // steps, start time,end time,initial condition,vector for answer
    virtual void solve(unsigned long long,T,T,const Point<T>&,std::vector<Point<T> >&) const;
    // steps, start time,end time,initial condition,T& for answer
    virtual void solveFinalVal(unsigned long long,T,T,const Point<T>&,Point<T>&) const;
  protected:
    const AbstractVectorFieldNonAut<T> * const vf;
  };
















}


#include "SimpleODEMethods.cpp"

#endif
