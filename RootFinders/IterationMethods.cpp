namespace NumericalAnalysis {


  template<typename T>
  NewtonMethod1d<T>::NewtonMethod1d(const RealFunction1d<T>& f) : rf(&f), MAX_STEPS_DEF(10), VALUE_TOL_DEF(.00001) {
    max_steps=MAX_STEPS_DEF;
    value_tol=VALUE_TOL_DEF;
  }

  template<typename T>
  int NewtonMethod1d<T>::solve(T x0,T& root) {
    int rv=0;
    unsigned long long counter=0;
    bool done=false;
    root=x0;
    T f_val;
    T df_val;
    while( ! done ) {
      ++counter;
      rf->eval(root,f_val);
      if( abs(f_val) < value_tol ) {
	done=true;
      }
      else {
	rf->eval_d(root,df_val);
	if(df_val==0) {
	  done=true;
	  rv=1;
	}
	root=root-f_val/df_val;
      }
    }
    
    return rv;
  }

  template<typename T>
  NewtonMethod<T>::NewtonMethod(unsigned long long d,const VectorField<T>& f) : vf(&f), MAX_STEPS_DEF(10), VALUE_TOL_DEF(.00001) {
    max_steps=MAX_STEPS_DEF;
    value_tol=VALUE_TOL_DEF;
    dim=d;
  }

  template<typename T>
  int NewtonMethod<T>::solve(Point<T> x0,Point<T>& root) {
    int rv=0;
    unsigned long long counter=0;
    bool done=false;
    root=x0;
    Point<T> f_val(dim);
    Point<T> dx(dim);
    T value_tol_squared=value_tol*value_tol;
    while( ! done ) {
      ++counter;
      vf->eval(root,f_val);
      if( f_val.norm_squared() < value_tol_squared ) {
	done=true;
	std::cout << "met value_tol " << value_tol<<std::endl;
      }
      else if (counter > max_steps) {
	done=true;
	std::cout << "reached max_steps " << counter<<std::endl;
      }
      else {
	if( vf->inv_d_eval(root,dx,f_val) ) { // error finding dx
	  done=true;
	  std::cout << "error finding dx"<<std::endl;
	  rv=1;
	}
	else {
	  for(unsigned long long i=0;i<dim;++i) {
	    root[i]=root[i]-dx[i];
	  }
	}
      }
    }
    
    return rv;
  }




















}
