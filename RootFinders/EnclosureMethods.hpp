#ifndef __BISECTIONMETHOD_HPP
#define __BISECTIONMETHOD_HPP

#include "../Math/Math.hpp"
#include "../Containers/Point.hpp"
#include "../Containers/Functions.hpp"



namespace NumericalAnalysis {

  /*

    ABSTRACT ROOT FINDER CLASS
    T is number type

  */
  template<typename T>
  class RootFinder {
  private:
    const double VALUE_TOL_DEF;
    const long long MAX_STEPS_DEF;
  protected:
    T value_tol;
    long long max_steps;
    const RealFunction1d<T>* rf;
  public:
    RootFinder(const RealFunction1d<T>*);
    virtual int solve(T,T,T&)=0;
  };

  /*

    BISECTION METHOD
    T is number type

   */
  template<typename T>
  class BisectionMethod : public RootFinder<T> {
  private:
  protected:
  public:
    BisectionMethod(const RealFunction1d<T>*);
    virtual int solve(T,T,T&);
  };


  
  /*

    FALSE POSITION METHOD
    T is number type

   */
  template<typename T>
  class FalsePositionMethod : public RootFinder<T> {
  private:
    T root_tol;
    bool use_root_tol;
    bool use_value_tol;
  protected:
  public:
    FalsePositionMethod(RealFunction1d<T>&);
    FalsePositionMethod(RealFunction1d<T>&,T);
    virtual int solve(T,T,T&);
  };








};







#include "EnclosureMethods.cpp"

#endif
