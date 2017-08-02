#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "../Sizes.hpp"
#include "Point.hpp"

namespace NumericalAnalysis {

  // typedefs for functions
  template<typename T>
  using RF_FunctionPointer=void (*const)(const Point<T>& x,T& rv);
  template<typename T>
  using RFNA_FunctionPointer=void (*const)(const T& t,const Point<T>& x,T& rv);
  // typedefs for vector fields
  template<typename T>
  using VF_FunctionPointer=void(*const)(const Point<T>& x,Point<T>& rv);
  template<typename T>
  using VFNA_FunctionPointer=void(*const)(const T& t,const Point<T>& x,Point<T>& rv);
  template<typename T>
  using DerivativeVF_FunctionPointer=void(*const)(const Point<T>& x,const Point<T>& v,Point<T>& rv);
  

  /*

    VECTOR FIELD CLASSES
  */

  // abstract class for vector fields
  // T is the number type
  template<typename T>
  class AbstractVectorField {
  public:
    virtual void operator() (const Point<T>& x,Point<T>& rv) const=0;
    dim_type getDim() const {return dim;}
  protected:
    AbstractVectorField(dim_type dim) : dim(dim) {}
    const dim_type dim;
  };

  template<typename T>
  class AbstractVectorFieldNonAut {
  public:
    virtual void operator() (const T& t,const Point<T>& x,Point<T>& rv) const=0;
    dim_type getDim() const {return dim;}
  protected:
    AbstractVectorFieldNonAut(dim_type dim) : dim(dim) {}
    const dim_type dim;
  };

  // autonomous vector field
  // T is the number type
  template<typename T>
  class VectorField : public AbstractVectorField<T> {
  public:
    VectorField(dim_type dim,VF_FunctionPointer<T> vf) : AbstractVectorField<T>(dim), vf(vf) {}
    virtual void operator() (const Point<T>& x,Point<T>& rv) const {vf(x,rv);}
  private:
    VF_FunctionPointer<T> vf;
  };

  // nonautonomous vector field
  // T is the number type
  template<typename T>
  class VectorFieldNonAut : public AbstractVectorFieldNonAut<T> {
  public:
    VectorFieldNonAut(dim_type dim,VFNA_FunctionPointer<T> vf) : AbstractVectorFieldNonAut<T>(dim), vf(vf) {}
    virtual void operator() (const T& t,const Point<T>& x,Point<T>& rv) const {vf(t,x,rv);}
  private:
    VFNA_FunctionPointer<T> vf;
  };

  // nonautonomous 1d vector field
  // T is the number type
  template<typename T>
  class VectorField1dNonAut : public VectorFieldNonAut<T> {
  public:
    VectorField1dNonAut(VFNA_FunctionPointer<T> vf) : VectorFieldNonAut<T>(1,vf) {}
  };


  // nonautonomous 2d vector field
  // T is the number type
  // to be used for converting 2nd order ODE to 1st order
  // for 1dim domain only
  template<typename T>
  class VFNA2dFrom2ndOrder : public AbstractVectorFieldNonAut<T> {
  public:
    VFNA2dFrom2ndOrder(RFNA_FunctionPointer<T> f) : AbstractVectorFieldNonAut<T>(2), f(f) {}
    virtual void operator() (const T& t,const Point<T>& x,Point<T>& rv) const {rv[0]=x[1];f(t,x,rv[1]);}
  protected:
    RFNA_FunctionPointer<T> f;
  };
  

  // abstract class for vector fields
  // T is the number type
  template<typename T>
  class AbstractDerivativeOfVectorField {
  public:
    virtual void operator() (const Point<T>& x,const Point<T>& v,Point<T>& rv) const=0;
    dim_type getDim() const {return dim;}
  protected:
    AbstractDerivativeOfVectorField(dim_type dim) : dim(dim) {}
    const dim_type dim;
  };




  // autonomous vector field
  // T is the number type
  template<typename T>
  class DerivativeOfVectorField : public AbstractDerivativeOfVectorField<T> {
  public:
    DerivativeOfVectorField(dim_type dim,DerivativeVF_FunctionPointer<T> dvf) : AbstractDerivativeOfVectorField<T>(dim), dvf(dvf) {}
    virtual void operator() (const Point<T>& x,const Point<T>& v,Point<T>& rv) const {dvf(x,v,rv);}
  private:
    DerivativeVF_FunctionPointer<T> dvf;
  };

  /*

    REAL VALUED FUNCTION CLASSES

  */

 
  // abstract real valued function
  // T is the number type
  template<typename T>
  class AbstractRealFunction {
  public:
    virtual void operator()(const Point<T>& x,T& rv) const=0;
    dim_type getDim() const {return dim;}
  protected:
    AbstractRealFunction(dim_type dim) : dim(dim) {}
    dim_type dim;
  };

  // "non-autonomous" real valued function (also takes a time t parameter)
  // T is the number type
  template<typename T>
  class AbstractRealFunctionNonAut {
  public:
    virtual void operator() (const T& t,const Point<T>& x,T& rv) const=0;
    dim_type getDim() const {return dim;}
  protected:
    dim_type dim;
    AbstractRealFunctionNonAut(dim_type dim) : dim(dim) {}
  };

 template<typename T>
 class RealFunction : public AbstractRealFunction<T> {
 public:
   RealFunction(dim_type dim,RF_FunctionPointer<T> f) : AbstractRealFunction<T>(dim), f(f) {}
   virtual void operator()(const Point<T>& x,T& rv) const {f(x,rv);}
 private:
   RF_FunctionPointer<T> f;
 };
  
  template<typename T>
  class RealFunctionNonAut : public AbstractRealFunctionNonAut<T> {
  public:
    RealFunctionNonAut(dim_type dim,RFNA_FunctionPointer<T> f) : AbstractRealFunctionNonAut<T>(dim), f(f) {}
    virtual void operator()(const T& t,const Point<T>& x,T& rv) const {f(t,x,rv);}
  protected:
    RFNA_FunctionPointer<T> f;
  };

  // 1d real valued function
  // T is the number type
  template<typename T>
  class RealFunction1d : public RealFunction<T> {
  public:
    RealFunction1d(RF_FunctionPointer<T> f) : RealFunction<T>(1,f) {}
  };


  template<typename T>
  class DifferentiableRealFunction1d : public RealFunction1d<T> {
  public:
    DifferentiableRealFunction1d(RF_FunctionPointer<T> f,RF_FunctionPointer<T> df) : RealFunction1d<T>(f), df(df) {}
    virtual void derivative(const Point<T>& x,T& rv) const {df(x,rv);}
  protected:
    RF_FunctionPointer<T> df;
  };





}



#endif
