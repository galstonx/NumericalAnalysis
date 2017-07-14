

namespace NumericalAnalysis {







  template<typename T>
  RootFinder<T>::RootFinder(const RealFunction1d<T>* f) : rf(f), VALUE_TOL_DEF(.000001), MAX_STEPS_DEF(100) {
    max_steps=MAX_STEPS_DEF;
    value_tol=VALUE_TOL_DEF;
  }


  template<typename T>
  BisectionMethod<T>::BisectionMethod(const RealFunction1d<T>* f) : RootFinder<T>(f) {
  }



  template<typename T>
  int BisectionMethod<T>::solve(T left,T right, T& root) {
    // needed for two-phase lookup
    T value_tol=RootFinder<T>::value_tol;
    long long max_steps=RootFinder<T>::max_steps;
    const RealFunction1d<T>* rf=RootFinder<T>::rf;

    int rv=0;
    long long counter=0;
    bool done=false;
    T& mid=root;
    T a=left;
    T b=right;
    T rf_left;
    T rf_mid;
    T rf_right;
    while(! done) {
      counter++;
      mid=(a+b)/2;
      rf->eval(a,rf_left);
      rf->eval(b,rf_right);
      rf->eval(mid,rf_mid);
      if( (abs(rf_mid)<value_tol) || (counter>max_steps) ) {
	done=true;
      }
      else {
	// midpoint is negative
	if(rf_mid<0) {
	  if(rf_left<0) {
	    a=mid;
	  }
	  else {
	    b=mid;
	  }
	}
	// midpoint is positive
	else {
	  if(rf_left>0) {
	    a=mid;
	  }
	  else {
	    b=mid;
	  }
	}

      }
    }
    rv=1;
  }


  template<typename T>
  FalsePositionMethod<T>::FalsePositionMethod(RealFunction1d<T>& f) : RootFinder<T>(f) {
    use_root_tol=false;
    use_value_tol=true;
  }

  template<typename T>
  FalsePositionMethod<T>::FalsePositionMethod(RealFunction1d<T>& f,T rt) : RootFinder<T>(f) {
    use_root_tol=true;
    use_value_tol=false;
    root_tol=rt;
  }



  template<typename T>
  int FalsePositionMethod<T>::solve(T left,T right, T& root) {
    // needed for two-phase lookup
    T value_tol=RootFinder<T>::value_tol;
    long long max_steps=RootFinder<T>::max_steps;
    RealFunction1d<T>* rf=RootFinder<T>::rf;
    int rv=0;
    long long counter=0;
    bool done=false;
    T a=left;
    T b=right;
    T tmp1;
    T tmp2;
    T rf_left;
    T rf_root;
    T rf_right;
    T rf_root1;
    T rf_root2;
    while(! done) {
      counter++;
      // main logic for finding next point is here
      rf->eval(a,rf_left);
      rf->eval(b,rf_right);
      tmp2=tmp1;
      tmp1=root;
      root=b-rf_right*(b-a)/(rf_right-rf_left);
      rf->eval(root,rf_root);
      // update endpoints
      // root is negative
      if(rf_root<0) {
	if(rf_left<0) {
	  a=root;
	}
	else {
	  b=root;
	}
      }
      // root is positive
      else {
	if(rf_left>0) {
	  a=root;
	}
	else {
	  b=root;
	}
      }
      // check if we are done
      // first check if we are using value_tol
      if(use_value_tol) {
	if( (abs(rf_root)<value_tol) || (counter>max_steps) ) {
	  done=true;
	}
      }
      // then check if we are using root_tol
      if(use_root_tol && counter>=3) {
	// should put a check here to make sure we are not dividing by zero
	T lambda=(root-tmp1)/(tmp1-tmp2);
	T error=abs(lambda/(lambda-1)*(root-tmp1));
	if(error<root_tol) {
	  done=true;
	}
      }
    }
    rv=1;
  }




};
