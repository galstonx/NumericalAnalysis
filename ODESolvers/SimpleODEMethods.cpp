namespace NumericalAnalysis {


  template<typename T>
  EulerMethod1d<T>::EulerMethod1d(VectorField1dNonAut<T> const * f) : vf(f) {
  }

  template<typename T>
  int EulerMethod1d<T>::solve(unsigned long long steps,T a,T b,T x0,std::vector<T>& rv) {
    if(steps+1 != rv.size() ) {
      return 1;
    }
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    rv[0]=x0;
    T f_val;
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,rv[i-1],f_val);
      rv[i]=rv[i-1]+dt*f_val;
      t=t+dt;
    }
    return 0;
  }


  template<typename T>
  int EulerMethod1d<T>::solveFinalVal(unsigned long long steps,T a,T b,T x0,T& ans) {
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    ans=x0;
    T f_val;
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,ans,f_val);
      ans=ans+dt*f_val;
      t=t+dt;
    }
    return 0;
  }
  
  template<typename T>
  EulerMethod<T>::EulerMethod(unsigned long long d,VectorFieldNonAut<T> const * f) : dim(d), vf(f) {
  }

  template<typename T>
  int EulerMethod<T>::solve(unsigned long long steps,T a,T b,const Point<T>& x0,std::vector<Point<T> >& rv) {
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    rv[0]=x0;
    Point<T> f_val(dim);
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,rv[i-1],f_val);
      for(unsigned long long j=0;j<dim;++j) {
	rv[i][j]=rv[i-1][j]+dt*f_val[j];
      }
      t=t+dt;
    }
    return 0;
  }


  template<typename T>
  int EulerMethod<T>::solveFinalVal(unsigned long long steps,T a,T b,const Point<T>& x0,Point<T>& ans) {
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    ans=x0;
    Point<T> f_val(dim);
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,ans,f_val);
      for(unsigned long long j=0;j<dim;++j) {
	ans[j]=ans[j]+dt*f_val[j];
      }
      t=t+dt;
    }
    return 0;
  }


  template<typename T>
  RungeKuttaMethod1d<T>::RungeKuttaMethod1d(VectorField1dNonAut<T> const * f) : vf(f) {
  }

  template<typename T>
  int RungeKuttaMethod1d<T>::solve(unsigned long long steps,T a,T b,T x0,std::vector<T>& rv) {
    if(steps+1 != rv.size() ) {
      return 1;
    }
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    rv[0]=x0;
    T k1;
    T k2;
    T k3;
    T k4;
    T x2;
    T x3;
    T x4;
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,rv[i-1],k1);
      x2=rv[i-1]+dt/2*k1;
      vf->eval(t+dt/2,x2,k2);
      x3=rv[i-1]+dt/2*k2;
      vf->eval(t+dt/2,x3,k3);
      x4=rv[i-1]+dt*k3;
      vf->eval(t+dt,x4,k4);
      rv[i]=rv[i-1]+dt/6*(k1+2*k2+2*k3+k4);
      t=t+dt;
    }
    return 0;
  }


  template<typename T>
  int RungeKuttaMethod1d<T>::solveFinalVal(unsigned long long steps,T a,T b,T x0,T& ans) {
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    ans=x0;
    T k1;
    T k2;
    T k3;
    T k4;
    T x2;
    T x3;
    T x4;
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,ans,k1);
      x2=ans+dt/2*k1;
      vf->eval(t+dt/2,x2,k2);
      x3=ans+dt/2*k2;
      vf->eval(t+dt/2,x3,k3);
      x4=ans+dt*k3;
      vf->eval(t+dt,x4,k4);
      ans=ans+dt/6*(k1+2*k2+2*k3+k4);
      t=t+dt;
    }
    return 0;
  }
  
  template<typename T>
  RungeKuttaMethod<T>::RungeKuttaMethod(unsigned long long d,VectorFieldNonAut<T> const * f) : dim(d), vf(f) {
  }

  template<typename T>
  int RungeKuttaMethod<T>::solve(unsigned long long steps,T a,T b,const Point<T>& x0,std::vector<Point<T> >& rv) {
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    rv[0]=x0;
    Point<T> k1(dim);
    Point<T> k2(dim);
    Point<T> k3(dim);
    Point<T> k4(dim);
    Point<T> x2(dim);
    Point<T> x3(dim);
    Point<T> x4(dim);
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,rv[i-1],k1);
      for(unsigned long long j=0;j<dim;++j) {
	x2[j]=rv[i-1][j]+dt/2*k1[j];
      }
      vf->eval(t+dt/2,x2,k2);
      for(unsigned long long j=0;j<dim;++j) {
	x3[j]=rv[i-1][j]+dt/2*k2[j];
      }
      vf->eval(t+dt/2,x3,k3);
      for(unsigned long long j=0;j<dim;++j) {
	x4[j]=rv[i-1][j]+dt*k3[j];
      }
      vf->eval(t+dt,x4,k4);
      for(unsigned long long j=0;j<dim;++j) {
	rv[i][j]=rv[i-1][j]+dt/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
      }
      t=t+dt;
    }
    return 0;
  }


  template<typename T>
  int RungeKuttaMethod<T>::solveFinalVal(unsigned long long steps,T a,T b,const Point<T>& x0,Point<T>& ans) {
    n_steps=steps;
    t_start=a;
    t_stop=b;
    dt=(b-a)/steps;
    ans=x0;
    Point<T> k1(dim);
    Point<T> k2(dim);
    Point<T> k3(dim);
    Point<T> k4(dim);
    Point<T> x2(dim);
    Point<T> x3(dim);
    Point<T> x4(dim);
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      vf->eval(t,ans,k1);
      x2.setLinearComb(1,ans,dt/2,k1);
      vf->eval(t+dt/2,x2,k2);
      x3.setLinearComb(1,ans,dt/2,k2);
      vf->eval(t+dt/2,x3,k3);
      x4.setLinearComb(1,ans,dt,k3);
      vf->eval(t+dt,x4,k4);
      for(unsigned long long j=0;j<dim;++j) {
	ans[j]=ans[j]+dt/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
      }
      t=t+dt;
    }
    return 0;
  }

}
