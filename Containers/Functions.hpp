#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP


#include "Point.hpp"

namespace NumericalAnalysis {

  // typedefs for function pointers
  template<typename T>
  using VectorFieldRaw=void(*)(const Point<T>&,Point<T>&);
  
  

  // autonomous vector field
  // T is the number type
  template<typename T>
  class VectorField {
  private:
    void (*vf)(const Point<T>&,Point<T>&);
    void (*d_vf)(const Point<T>&,const Point<T>&,Point<T>&);
    void (*inv_d_vf)(const Point<T>&,Point<T>&,const Point<T>&);
    bool is_diff;
    unsigned long long dim;
  protected:
    VectorField() {};
  public:
    VectorField(unsigned long long,void(*)(const Point<T>&,Point<T>&));
    void set_d_vf(void (*)(const Point<T>&,const Point<T>&,Point<T>&));
    void set_inv_d_vf(void (*)(const Point<T>&,Point<T>&,const Point<T>&));
    void eval(const Point<T>&,Point<T>&) const;
    void d_eval(const Point<T>&,const Point<T>&,Point<T>&) const;
    int inv_d_eval(const Point<T>&,Point<T>&,const Point<T>&) const;
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
