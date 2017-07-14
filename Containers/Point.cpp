#include "Point.hpp"


namespace NumericalAnalysis {



  template<typename T>
  Point<T>::Point() {}
  
  template<typename T>
  Point<T>::Point(int size,T t) : p_data(size,t) {}

  template<typename T>
  Point<T>::Point(int size) : p_data(size) {}

  template<typename T>
  T& Point<T>::operator[] (int i) {
    return p_data[i];
  }

  template<typename T>
  const T& Point<T>::operator[] (int i) const {
    return p_data[i];
  }


  template<typename T>
  int midpoint(const Point<T>& a,const Point<T>& b,Point<T>& mid) {
    int rv=1;
    int dim=a.p_data.size();
    if(dim==b.p_data.size() && dim==mid.p_data.size()) {
      for(int i=0;i<dim;++i) {
	mid[i]=(a[i]+b[i])/2;
      }
      rv=0;
    }
    return rv;
  }

  template<typename T>
  T Point<T>::norm_squared() const {
    T rv;
    for(int i=0;i<p_data.size();++i) {
      rv+=p_data[i]*p_data[i];
    }
    return rv;
  }

  template<typename T>
  void Point<T>::setLinearComb(T c1,const Point<T>& p1,T c2,const Point<T>& p2) {
    for(unsigned long long i=0;i<p_data.size();++i) {
      p_data[i]=c1*p1[i]+c2*p2[i];
    }
  }
};
