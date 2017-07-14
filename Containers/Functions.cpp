
namespace NumericalAnalysis {

  /*

    VectorField

  */
  template<typename T>
  void VectorField<T>::eval(const Point<T>& x,Point<T>& rv) const {
    vf(x,rv);
  }

  template<typename T>
  void VectorField<T>::d_eval(const Point<T>& x,const Point<T>& dx,Point<T>& rv) const {
    d_vf(x,dx,rv);
  }

  template<typename T>
  int VectorField<T>::inv_d_eval(const Point<T>& x,Point<T>& dx,const Point<T>& rv) const {
    inv_d_vf(x,dx,rv);
    return 0;
  }
  
  template<typename T>
  VectorField<T>::VectorField(unsigned long long d,void (*f)(const Point<T>&,Point<T>&)) {
    vf=f;
    dim=d;
  }

  template<typename T>
  void VectorField<T>::set_d_vf(void (*df)(const Point<T>&,const Point<T>&,Point<T>&)) {
    d_vf=df;
  }

  template<typename T>
  void VectorField<T>::set_inv_d_vf(void (*df)(const Point<T>&,Point<T>&,const Point<T>&)) {
    inv_d_vf=df;
  }

  /*

    RealFunction

  */
  template<typename T>
  void RealFunction<T>::eval(const Point<T>& x,T& rv) {
    func(x,rv);
  }
  
  template<typename T>
  RealFunction<T>::RealFunction(void (*f)(const Point<T>&,T&)) {
    func=f;
  }
  
  /*

    RealFunctionNonAut

  */
  template<typename T>
  void RealFunctionNonAut<T>::eval(T t, const Point<T>& x,T& rv) const {
    func(t,x,rv);
  }
  
  template<typename T>
  RealFunctionNonAut<T>::RealFunctionNonAut(void (*f)(T,const Point<T>&,T&)) {
    func=f;
  }

  /*

    RealFunction1d

  */
  template<typename T>
  void RealFunction1d<T>::eval(T x,T& rv) const {
    func(x,rv);
  }
   
  template<typename T>
  void RealFunction1d<T>::eval_d(T x,T& rv) const {
    d_func(x,rv);
  }

  template<typename T>
  RealFunction1d<T>::RealFunction1d(void (*f)(T,T&)) {
    func=f;
    is_diff=false;
  }
  
  template<typename T>
  RealFunction1d<T>::RealFunction1d(void (*f)(T,T&),void (*df)(T,T&)) {
    func=f;
    d_func=df;
    is_diff=true;
  }
  

  /*

    VectorField1dNonAut

  */

  template<typename T>
  VectorField1dNonAut<T>::VectorField1dNonAut(void(*f)(T,T,T&)) {
    vf=f;
  }

  template<typename T>
  void VectorField1dNonAut<T>::eval(T t,T x,T& rv) const {
    vf(t,x,rv);
  }

 /*

    VectorFieldNonAut

  */

  template<typename T>
  VectorFieldNonAut<T>::VectorFieldNonAut(void(*f)(T,const Point<T>&,Point<T>&)) {
    vf=f;
  }

  template<typename T>
  void VectorFieldNonAut<T>::eval(T t,const Point<T>& x,Point<T>& rv) const {
    vf(t,x,rv);
  }
  


  /*
    
    VF2dFrom2ndOrder
    
  */
  
  template<typename T>
  VF2dFrom2ndOrder<T>::VF2dFrom2ndOrder(RealFunctionNonAut<T> const* f) {
    F=f;
  }
  
  template<typename T>
  void VF2dFrom2ndOrder<T>::eval(T t,const Point<T>& x,Point<T>& rv) const {
    rv[0]=x[1];
    F->eval(t,x,rv[1]);
  }
}
