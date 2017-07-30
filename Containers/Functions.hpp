#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "../Sizes.hpp"
#include "Point.hpp"

namespace NumericalAnalysis {

  // typedefs for function pointers
  template<typename T>
  using VF_FunctionPointer=void(*const)(const Point<T>& x,Point<T>& rv);
  
  
  // abstract class for vector fields
  // T is the number type
  template<typename T>
  class AbstractVectorField {
  public:
    virtual void operator() (const Point<T>& x,Point<T>& rv) const=0;
  protected:
    AbstractVectorField(dim_type dim) : dim(dim) {}
    const dim_type dim;
  };

  template<typename T>
  class AbstractVectorFieldNonAut {
  public:
    virtual void operator() (const T& t,const Point<T>& x,Point<T>& rv) const=0;
  protected:
    dim_type dim;
  };

  // autonomous vector field
  // T is the number type
  template<typename T>
  class VectorField : public AbstractVectorField<T> {
  public:
    VectorField(dim_type dim,VF_FunctionPointer<T> vf);
    void operator() (const Point<T>& x,Point<T>& rv) const;
  private:
    VF_FunctionPointer<T> vf;
  };
 
  // nonautonomous vector field
  // T is the number type
  template<typename T>
  class VectorFieldNonAut {
  private:
    // first argument is time, then position, then return value
    void (*vf)(T,const Point<T>&,Point<T>&);
  protected:
    VectorFieldNonAut() {};
  public:
    VectorFieldNonAut(void(*)(T,const Point<T>&,Point<T>&));
    virtual void eval(T,const Point<T>&,Point<T>&) const;
  };

  // nonautonomous 1d vector field
  // T is the number type
  template<typename T>
  class VectorField1dNonAut {
  private:
    // first argument is time, then position, then return value
    void (*vf)(T,T,T&);
  protected:
    VectorField1dNonAut() {};
  public:
    VectorField1dNonAut(void(*)(T,T,T&));
    virtual void eval(T,T,T&) const;
  };



  // real valued function
  // T is the number type
  template<typename T>
  class RealFunction {
  private:
    void (*func)(const Point<T>&,T&);
  protected:
    RealFunction() {};
  public:
    RealFunction(void (*) (const Point<T>&,T&));
    virtual void eval(const Point<T>&,T&);
  };

  // "non-autonomous" real valued function (also takes a time t parameter)
  // T is the number type
  //
  template<typename T>
  class RealFunctionNonAut {
  private:
    void (*func)(T,const Point<T>&,T&);
  protected:
    RealFunctionNonAut() {};
  public:
    RealFunctionNonAut(void (*) (T,const Point<T>&,T&));
    virtual void eval(T,const Point<T>&,T&) const;
  };

  // 1d real valued function
  // T is the number type
  template<typename T>
  class RealFunction1d {
  private:
    void (*func)(T,T&);
    void (*d_func)(T,T&);
    bool is_diff;
  protected:
    RealFunction1d() {};
  public:
    RealFunction1d(void (*)(T,T&));
    RealFunction1d(void (*)(T,T&),void (*)(T,T&));
    virtual void eval(T,T&) const;
    virtual void eval_d(T,T&) const;
  };


  // nonautonomous 2d vector field
  // T is the number type
  // to be used for converting 2nd order ODE to 1st order
  template<typename T>
  class VF2dFrom2ndOrder : public VectorFieldNonAut<T> {
  private:
  protected:
    RealFunctionNonAut<T> const * F; // RHS of 2nd order ODE
  public:
    VF2dFrom2ndOrder(RealFunctionNonAut<T> const *);
    virtual void eval(T,const Point<T>&,Point<T>&) const;
  };




}


#include "Functions.cpp"




#endif
