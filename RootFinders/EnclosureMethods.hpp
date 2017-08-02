#ifndef ENCLOSURE_METHODS_HPP
#define ENCLOSURE_METHODS_HPP

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
  public:
    RootFinder(const AbstractRealFunction<T>* f) : rf(f), value_tol(.00001), max_steps(100) {}
    virtual void solve(T,T,T&)=0;
  protected:
    const AbstractRealFunction<T>* rf;
    T value_tol;
    long long max_steps;
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
    BisectionMethod(const AbstractRealFunction<T>* f) : RootFinder<T>(f) {}
    virtual void solve(T,T,T&);
  };


  
  /*

    FALSE POSITION METHOD
    T is number type

   */
  template<typename T>
  class FalsePositionMethod : public RootFinder<T> {
  public:
    FalsePositionMethod(RealFunction1d<T>* f) : RootFinder<T>(f), use_root_tol(false),  use_value_tol(true) {}
    FalsePositionMethod(RealFunction1d<T>* f,T root_tol) : RootFinder<T>(f), use_root_tol(true),  use_value_tol(false), root_tol(root_tol) {}
    virtual void solve(T,T,T&);
  protected:
    bool use_root_tol;
    bool use_value_tol;
    T root_tol;
  };








};







#include "EnclosureMethods.cpp"

#endif
