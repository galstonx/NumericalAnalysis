#ifndef __POINT_HPP
#define __POINT_HPP

#include <vector>
#include <iostream>


namespace NumericalAnalysis {



  // simple class to store a point in Euclidean space
  // the template parameter T is the type of a number
  // e.g. double, or other user defined number
  // needs to have arithmetic operations defined
  template<typename T>
  class Point {
  private:
    std::vector<T> p_data;
  protected:
    Point();
  public:
    Point(int);
    Point(int,T);
    T norm_squared() const;
    T& operator[](int);
    const T& operator[](int) const;
    void setLinearComb(T,const Point<T>&,T,const Point<T>&);
    template<typename S> friend int midpoint(const Point<S>&,const Point<S>&,Point<S>&);
  };
    



}



#include "Point.cpp"

#endif
