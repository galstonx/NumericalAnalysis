#ifndef POINT_HPP
#define POINT_HPP

#include <vector>
#include <iostream>

#include "../Sizes.hpp"

namespace NumericalAnalysis {



  // simple class to store a point in Euclidean space
  // the template parameter T is the type of a number
  template<typename T>
  class Point {
  public:
    Point(dim_type dim) : data(dim) {}
    Point(dim_type dim,T t) : data(dim,t) {}
    T norm_squared() const {
      T rv(0);
      for(auto t:data) rv+=t*t;
      return rv;
    }
    T& operator[](dim_type i) {return data[i];}
    const T& operator[](dim_type i) const {return data[i];}
    dim_type size() const {return data.size();}
    //template<typename S> friend std::ostream& operator<< (std::ostream& os,const Point<S>& p);
  protected:
  private:
    std::vector<T> data;
  };
    

  template<typename T> std::ostream& operator<< (std::ostream& os,const Point<T>& p) {
    dim_type dim=p.size();
    for(dim_type i=0;i<dim-1;++i) {
      os << p[i] << ',';
    }
    if(dim>0) os << p[dim-1];
    return os;
  }

}


#endif
