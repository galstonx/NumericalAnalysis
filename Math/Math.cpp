namespace NumericalAnalysis {



  template<typename T>
  T abs(T x) {
    T rv;
    if(x<0) rv=-x;
    else rv=x;
    return rv;
  }


};
