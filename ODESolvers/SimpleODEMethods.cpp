namespace NumericalAnalysis {


  template<typename T>
  void EulerMethod<T>::solve(unsigned long long n_steps,T t_start,T t_end,const Point<T>& x0,std::vector<Point<T> >& rv) const {
    dim_type dim=ODEMethod<T>::dim;
    T dt=(t_end-t_start)/n_steps;
    rv[0]=x0;
    Point<T> f_val(dim);
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      (*vf)(t,rv[i-1],f_val);
      for(unsigned long long j=0;j<dim;++j) {
	rv[i][j]=rv[i-1][j]+dt*f_val[j];
      }
      t=t+dt;
    }
  }


  template<typename T>
  void EulerMethod<T>::solveFinalVal(unsigned long long n_steps,T t_start,T t_end,const Point<T>& x0,Point<T>& ans) const {
    dim_type dim=ODEMethod<T>::dim;
    T dt=(t_end-t_start)/n_steps;
    ans=x0;
    Point<T> f_val(dim);
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      (*vf)(t,ans,f_val);
      for(unsigned long long j=0;j<dim;++j) {
	ans[j]=ans[j]+dt*f_val[j];
      }
      t=t+dt;
    }
  }




  template<typename T>
  void RungeKuttaMethod<T>::solve(unsigned long long n_steps,T t_start,T t_end,const Point<T>& x0,std::vector<Point<T> >& rv) const {
    dim_type dim=ODEMethod<T>::dim;
    T dt=(t_end-t_start)/n_steps;
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
      (*vf)(t,rv[i-1],k1);
      for(unsigned long long j=0;j<dim;++j) {
	x2[j]=rv[i-1][j]+dt/2*k1[j];
      }
      (*vf)(t+dt/2,x2,k2);
      for(unsigned long long j=0;j<dim;++j) {
	x3[j]=rv[i-1][j]+dt/2*k2[j];
      }
      (*vf)(t+dt/2,x3,k3);
      for(unsigned long long j=0;j<dim;++j) {
	x4[j]=rv[i-1][j]+dt*k3[j];
      }
      (*vf)(t+dt,x4,k4);
      for(unsigned long long j=0;j<dim;++j) {
	rv[i][j]=rv[i-1][j]+dt/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
      }
      t=t+dt;
    }
  }


  template<typename T>
  void RungeKuttaMethod<T>::solveFinalVal(unsigned long long n_steps,T t_start,T t_end,const Point<T>& x0,Point<T>& rv) const {
    dim_type dim=ODEMethod<T>::dim;
    T dt=(t_end-t_start)/n_steps;
    rv=x0;
    Point<T> k1(dim);
    Point<T> k2(dim);
    Point<T> k3(dim);
    Point<T> k4(dim);
    Point<T> x2(dim);
    Point<T> x3(dim);
    Point<T> x4(dim);
    T t=t_start;
    for(unsigned long long i=1;i<n_steps+1;++i) {
      (*vf)(t,rv,k1);
      for(unsigned long long j=0;j<dim;++j) {
	x2[j]=rv[j]+dt/2*k1[j];
      }
      (*vf)(t+dt/2,x2,k2);
      for(unsigned long long j=0;j<dim;++j) {
	x3[j]=rv[j]+dt/2*k2[j];
      }
      (*vf)(t+dt/2,x3,k3);
      for(unsigned long long j=0;j<dim;++j) {
	x4[j]=rv[j]+dt*k3[j];
      }
      (*vf)(t+dt,x4,k4);
      for(unsigned long long j=0;j<dim;++j) {
	rv[j]=rv[j]+dt/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
      }
      t=t+dt;
    }
  }

}
