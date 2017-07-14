#ifndef __SIMPLE_ODE_METHODS_HPP
#define __SIMPLE_ODE_METHODS_HPP

#include "../Containers/Functions.hpp"

namespace NumericalAnalysis {

  template<typename T>
  class ODEMethod1d {
  private:
  protected:
  public:
    // steps, start time,end time,initial condition,vector for answer
    virtual int solve(unsigned long long,T,T,T,std::vector<T>&)=0;
    // steps, start time,end time,initial condition,T& for answer
    virtual int solveFinalVal(unsigned long long,T,T,T,T&)=0;
  };


  template<typename T>
  class EulerMethod1d : public ODEMethod1d<T> {
  private:
  protected:
    unsigned long long n_steps;
    T t_start;
    T t_stop;
    T dt;
    VectorField1dNonAut<T> const * vf;
  public:
    EulerMethod1d(VectorField1dNonAut<T> const * f);
    // steps, start time,end time,initial condition,vector for answer
    int solve(unsigned long long,T,T,T,std::vector<T>&);
    // steps, start time,end time,initial condition,T& for answer
    int solveFinalVal(unsigned long long,T,T,T,T&);
  };


  template<typename T>
  class ODEMethod {
  private:
  protected:
  public:
    // steps, start time,end time,initial condition,vector for answer
    virtual int solve(unsigned long long,T,T,const Point<T>&,std::vector<Point<T> >&)=0;
    // steps, start time,end time,initial condition,T& for answer
    virtual int solveFinalVal(unsigned long long,T,T,const Point<T>&,Point<T>&)=0;
  };

  template<typename T>
  class EulerMethod : public ODEMethod<T> {
  private:
  protected:
    unsigned long long n_steps;
    T t_start;
    T t_stop;
    T dt;
    VectorFieldNonAut<T> const * vf;
    unsigned long long dim;
  public:
    EulerMethod(unsigned long long, VectorFieldNonAut<T> const * f);
    // steps, start time,end time,initial condition,vector for answer
    int solve(unsigned long long,T,T,const Point<T>&,std::vector<Point<T> >&);
    // steps, start time,end time,initial condition,T& for answer
    int solveFinalVal(unsigned long long,T,T,const Point<T>&,Point<T>&);
  };



 template<typename T>
 class RungeKuttaMethod1d : public ODEMethod1d<T> {
  private:
  protected:
    unsigned long long n_steps;
    T t_start;
    T t_stop;
    T dt;
    VectorField1dNonAut<T> const * vf;
  public:
    RungeKuttaMethod1d(VectorField1dNonAut<T> const * f);
    // steps, start time,end time,initial condition,vector for answer
    int solve(unsigned long long,T,T,T,std::vector<T>&);
    // steps, start time,end time,initial condition,T& for answer
    int solveFinalVal(unsigned long long,T,T,T,T&);
  };


  template<typename T>
  class RungeKuttaMethod : public ODEMethod<T> {
  private:
  protected:
    unsigned long long n_steps;
    T t_start;
    T t_stop;
    T dt;
    VectorFieldNonAut<T> const * vf;
    unsigned long long dim;
  public:
    RungeKuttaMethod(unsigned long long, VectorFieldNonAut<T> const * f);
    // steps, start time,end time,initial condition,vector for answer
    int solve(unsigned long long,T,T,const Point<T>&,std::vector<Point<T> >&);
    // steps, start time,end time,initial condition,T& for answer
    int solveFinalVal(unsigned long long,T,T,const Point<T>&,Point<T>&);
  };
















}


#include "SimpleODEMethods.cpp"

#endif
