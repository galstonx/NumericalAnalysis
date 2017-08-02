

namespace NumericalAnalysis {










  template<typename T>
  void BisectionMethod<T>::solve(T left,T right, T& root) {
    // needed for two-phase lookup
    T value_tol=RootFinder<T>::value_tol;
    long long max_steps=RootFinder<T>::max_steps;
    long long counter=0;
    bool done=false;
    Point<T> a(1,left);
    Point<T> b(1,right);
    Point<T> mid(1);
    T rf_left;
    T rf_mid;
    T rf_right;
    while(! done) {
      counter++;
      mid[0]=(a[0]+b[0])/2;
      (*this->rf)(a,rf_left);
      (*this->rf)(b,rf_right);
      (*this->rf)(mid,rf_mid);
      if( (abs(rf_mid)<value_tol) || (counter>max_steps) ) {
	done=true;
      }
      else {
	// midpoint is negative
	if(rf_mid<0) {
	  if(rf_left<0) {
	    a[0]=mid[0];
	  }
	  else {
	    b[0]=mid[0];
	  }
	}
	// midpoint is positive
	else {
	  if(rf_left>0) {
	    a[0]=mid[0];
	  }
	  else {
	    b[0]=mid[0];
	  }
	}

      }
    }
    root=mid[0];
  }




  template<typename T>
  void FalsePositionMethod<T>::solve(T left,T right, T& root) {
    // needed for two-phase lookup
    T value_tol=RootFinder<T>::value_tol;
    long long max_steps=RootFinder<T>::max_steps;
    const AbstractRealFunction<T>* rf=RootFinder<T>::rf; // actually is a RealFunction1d
    long long counter=0;
    bool done=false;
    Point<T> a(1,left);
    Point<T> b(1,right);
    Point<T> root_pt(1,root);
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
      (*rf)(a,rf_left);
      (*rf)(b,rf_right);
      tmp2=tmp1;
      tmp1=root;
      root=b[0]-rf_right*(b[0]-a[0])/(rf_right-rf_left);
      (*rf)(root,rf_root);
      // update endpoints
      // root is negative
      if(rf_root<0) {
	if(rf_left<0) {
	  a[0]=root;
	  root_pt[0]=root;
	}
	else {
	  b[0]=root;
	  root_pt[0]=root;
	}
      }
      // root is positive
      else {
	if(rf_left>0) {
	  a[0]=root;
	  root_pt[0]=root;
	}
	else {
	  b[0]=root;
	  root_pt[0]=root;
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
  }




};
