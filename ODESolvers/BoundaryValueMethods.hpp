#ifndef BOUNDARY_VALUE_METHODS_HPP
#define BOUNDARY_VALUE_METHODS_HPP

#include "../Containers/Functions.hpp"
#include "../RootFinders/EnclosureMethods.hpp"

namespace NumericalAnalysis {

  /*

    Helper classes

  */

  template<typename T>
  class RealFunctionFrom2ndOrderODE : public AbstractRealFunction<T> {
  public:
    RealFunctionFrom2ndOrderODE(RFNA_FunctionPointer<T> f,unsigned long long steps,T t_start,T t_stop,T x_start,T x_stop) : AbstractRealFunction<T>(1), vf(f), t_start(t_start),t_stop(t_stop),x_start(x_start),x_stop(x_stop),steps(steps),ode_method(&vf) {}
    // given initial velocity, find position at t_stop minus x_stop
    virtual void operator()(const Point<T>& vel,T& rv) const {
      Point<T> x(2);
      x[0]=x_start;
      x[1]=vel[0];
      Point<T> final_x(2);
      ode_method.solveFinalVal(steps,t_start,t_stop,x,final_x);
      rv=final_x[0]-x_stop;
    }
  protected:
    const VFNA2dFrom2ndOrder<T> vf; // first order vf determined by F
    T t_start;
    T t_stop;    
    T x_start;
    T x_stop;
    unsigned long long steps;
    RungeKuttaMethod<T> ode_method;
  };

  template<typename T>
  class RealFunctionForEigenvalue : public AbstractRealFunction<T> {
  public:
    RealFunctionForEigenvalue(RFNA_FunctionPointer<T> f,unsigned long long steps,T t_start,T t_stop) : AbstractRealFunction<T>(1), f(f), t_start(t_start),t_stop(t_stop), steps(steps) {}
    // given lambada and intial vel 1, find position at t_stop
    virtual void operator()(const Point<T>& lambda,T& rv) const {
      //RFNA_FunctionPointer<T> f2(0);
      VFNA_2nd_order_to_1st_order<T> vf(2,f,lambda[0]);
      RungeKuttaMethod<T> rkm(&vf);
      Point<T> x(2);
      x[0]=0;
      x[1]=1;
      Point<T> final_x(2);
      rkm.solveFinalVal(steps,t_start,t_stop,x,final_x);
      rv=final_x[0];
    }
  protected:
    RFNA_FunctionPointer<T> f;
    T t_start;
    T t_stop;
    unsigned long long steps;
    template<typename S>
    class VFNA_2nd_order_to_1st_order : public AbstractVectorFieldNonAut<S> {
    public:
      VFNA_2nd_order_to_1st_order(dim_type dim,RFNA_FunctionPointer<S> f,T lambda) : AbstractVectorFieldNonAut<S>(dim), f(f), lambda(lambda) {}
      void operator() (const S& t,const Point<T>& x,Point<S>& rv) const {rv[0]=x[1];f(t,x,rv[1]);rv[1]+=lambda*x[0];}
    protected:
      RFNA_FunctionPointer<S> f;
      T lambda;
    };
  };
  



  /*

    Shooting Method

  */

  template<typename T>
  class ShootingMethodNonlinear1d {
  public:
    ShootingMethodNonlinear1d(RFNA_FunctionPointer<T> f,unsigned long long steps,T t_start,T t_stop,T x_start,T x_stop) :  F(f), steps_ode_method(steps),f(f,steps,t_start,t_stop,x_start,x_stop) {}
    // get the value of the function with a certain velocity
    void getFunctionValue(T vel,T& rv) const {f(Point<T>(1,vel),rv);}
    // Dirichlet using Bisection method
    // vel0,vel1 are initial velocities which make f bigger than 0 and smaller than 0 (or vice verse), rv is return value
    void solveDirichletBisection(T vel0,T vel1,T& rv) const {
      BisectionMethod<T> bm(&f);
      bm.solve(vel0,vel1,rv);
    }
  protected:
    RFNA_FunctionPointer<T> F; // RHS of 2nd order ODE
    unsigned long long steps_ode_method;
    const RealFunctionFrom2ndOrderODE<T> f; // (solution of ODE at t_stop)-x_stop, input=intial velocity, initial po=x_start
  };

  /*
    
    Eigenvalue Methods

  */

  template<typename T>
  class EigenvalueFor2ndOrder {
  public:
    EigenvalueFor2ndOrder(RFNA_FunctionPointer<T> f,unsigned long long steps_ode_method,T t_start,T t_stop) : f(f,steps_ode_method,t_start,t_stop), value_tol(.0001) {}
    void getEigenvalues(unsigned long long steps,T ev_begin,T ev_end,std::vector<T>& rv) {
      T h=(ev_end-ev_begin)/steps;
      Point<T> lambda(1);
      lambda[0]=ev_begin;
      T ans;
      bool skip=false;
      for(int i=0;i<=steps;++i) {
	f(lambda,ans);
	if(abs(ans)<value_tol) {
	  // skip consecutive entries
	  if(! skip) rv.push_back(lambda[0]);
	  skip=true; 
	}
	else {
	  skip=false;
	}
	lambda[0]=lambda[0]+h;
      }
    }
  protected: 
    RealFunctionForEigenvalue<T> f;
    T value_tol;
  };



}





#endif
