
namespace NumericalAnalysis {

  /*

    RealFunction1dForShooting

  */
  template<typename T>
  RealFunction1dForShooting<T>::RealFunction1dForShooting(const RealFunctionNonAut<T>* f,unsigned long long s,T a,T b,T x_a,T x_b) : vf(f) {
    steps=s;
    t_start=a;
    t_stop=b;
    x_start=x_a;
    x_stop=x_b;
    ode_method=new RungeKuttaMethod<T>(2,&vf);
  }

  template<typename T>
  RealFunction1dForShooting<T>::~RealFunction1dForShooting() {
    delete(ode_method);
  }

  template<typename T>
  void RealFunction1dForShooting<T>::eval(T vel,T& rv) const {
    Point<T> x(2);
    x[0]=x_start;
    x[1]=vel;
    Point<T> final_x(2);
    ode_method->solveFinalVal(steps,t_start,t_stop,x,final_x);
    rv=final_x[0]-x_stop;
  }

  /*

    ShootingMethodNonlinear1d

  */

  template<typename T>
  ShootingMethodNonlinear1d<T>::ShootingMethodNonlinear1d(const RealFunctionNonAut<T>* f,unsigned long long s,T a,T b,T x_a,T x_b) : rf(f,s,a,b,x_a,x_b) {
    t_start=a;
    t_stop=b;
    x_start=x_a;
    x_stop=x_b;
    steps_ode_method=s;
    F=f;
  }


  template<typename T>
  int ShootingMethodNonlinear1d<T>::solveDirichletBisection(T a,T b,T& rv) const {
    BisectionMethod<T> bm(&rf);
    bm.solve(a,b,rv);
    return 0;
  }


  template<typename T>
  int ShootingMethodNonlinear1d<T>::getFunctionValue(T vel,T& rv) const {
    rf.eval(vel,rv);
    return 0;
  }

  /*

    RealFunction1dForEigenvalue

  */

  template<typename T>
  RealFunction1dForEigenvalue<T>::RealFunction1dForEigenvalue(const RealFunctionNonAut<T>* f,unsigned long long s,T a,T b) : vf(f), ic(2) {
    steps=s;
    t_start=a;
    t_stop=b;
    ode_method=new RungeKuttaMethod<T>(2,&vf);
    vel_def=1.0;
    ic[0]=t_start;
    ic[1]=vel_def;
  }

  template<typename T>
  RealFunction1dForEigenvalue<T>::~RealFunction1dForEigenvalue() {
    delete(ode_method);
  }

  template<typename T>
  void RealFunction1dForEigenvalue<T>::eval(T& rv) const {
    Point<T> ans(2);
    ode_method->solveFinalVal(steps,t_start,t_stop,ic,ans);
    rv=ans[0];
  }
  
  template<typename T>
  int RealFunction1dForEigenvalue<T>::setEigenvalue(T l) {	
    lambda=l;
    vf.setEigenvalue(l);
    return 0;
  }
  /*
    
    VF2dForEigenvalue
    
  */
  
  template<typename T>
  VF2dForEigenvalue<T>::VF2dForEigenvalue(RealFunctionNonAut<T> const* f) {
    F=f;
    lambda=0;
  }

  template<typename T>
  int VF2dForEigenvalue<T>::setEigenvalue(T new_lambda) {
    lambda=new_lambda;
  }
  
  template<typename T>
  void VF2dForEigenvalue<T>::eval(T t,const Point<T>& x,Point<T>& rv) const {
    rv[0]=x[1];
    F->eval(t,x,rv[1]);
    rv[1]=rv[1]+lambda*x[0];
  }

  /*

    EigenvalueFor2ndOrder

  */

  template<typename T>
  EigenvalueFor2ndOrder<T>::EigenvalueFor2ndOrder(const RealFunctionNonAut<T>* f,unsigned long long s,T a,T b) : rf(f,s,a,b) {
    t_start=a;
    t_stop=b;
    steps_ode_method=s;
    F=f;
    value_tol=.0001;
  }

  template<typename T>
  int EigenvalueFor2ndOrder<T>::getEigenvalues(unsigned long long steps,T ev_begin,T ev_end,std::vector<T>& rv) {
    T h=(ev_end-ev_begin)/steps;
    T lambda=ev_begin;
    T ans;
    bool skip=false;
    for(int i=0;i<=steps;++i) {
      rf.setEigenvalue(lambda);
      rf.eval(ans);
      if(abs(ans)<value_tol) {
	// skip consecutive entries
	if(! skip) rv.push_back(lambda);
	skip=true; 
      }
      else {
	skip=false;
      }
      lambda=lambda+h;
    }
    return 0;
  }





}
