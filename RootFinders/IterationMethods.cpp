namespace NumericalAnalysis {


  template<typename T>
  void NewtonMethod1d<T>::solve(const T& x0,T& root) const {
    unsigned long long counter=0;
    bool done=false;
    Point<T> root_pt(1,x0);
    root=x0;
    T f_val;
    T df_val;
    while( ! done ) {
      ++counter;
      (*f)(root_pt,f_val);
      if( abs(f_val) < value_tol ) {
	done=true;
      }
      else {
	(*f).derivative(root_pt,df_val);
	if(df_val==0) {
	  done=true;
	}
	root=root-f_val/df_val;
	root_pt[0]=root;
      }
    }
  }

  template<typename T>
  void NewtonMethod<T>::solve(const Point<T>& x0,Point<T>& root) const {
    unsigned long long counter=0;
    bool done=false;
    root=x0;
    Point<T> f_val(dim);
    Point<T> dx(dim);
    T value_tol_squared=value_tol*value_tol;
    while( ! done ) {
      ++counter;
      (*vf)(root,f_val);
      if( f_val.norm_squared() < value_tol_squared ) {
	done=true;
	std::cout << "met value_tol " << value_tol<<std::endl;
      }
      else if (counter > max_steps) {
	done=true;
	std::cout << "reached max_steps " << counter<<std::endl;
      }
      else {
	(*dvf_inverse)(root,f_val,dx);
	for(unsigned long long i=0;i<dim;++i) {
	  root[i]=root[i]-dx[i];
	}
      }
    }
  }
  



















}
