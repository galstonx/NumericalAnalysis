#include <iostream>
#include "NumericalAnalysis.hpp"
using namespace std;
using namespace NumericalAnalysis;


void exp_vf(const Point<double>& x,Point<double>& rv) {
  rv=x;
}

void real_func(const Point<double>& x,double& rv) {
  rv=x[0]*x[0]-2;
}
void real_func_non_aut(double t,const Point<double>& x,double& rv) {
  rv=t+x[0]+x[1];
}
void F(double t,const Point<double>& x,double& rv) {
  rv=t+x[0]+x[1];
}

int main() {
  // test Point
  Point<double> v(5);
  v[2]=3;
  cout << v[2] << endl;

  // test VectorField
  VectorField<double> vf(2,exp_vf);
  Point<double> p(1);
  p[0]=-4;
  Point<double> rv(1);
  vf.eval(p,rv);
  cout<< rv[0]<<endl;

  //test RealFunction
  cout<<endl<<"Test RealFunction"<<endl;
  RealFunction<double> rf(real_func);
  Point<double> x(2);
  x[0]=5;
  x[1]=6;
  double rv2;
  rf.eval(x,rv2);
  cout<<rv2<<endl;

  //test RealFunctionNonAut
  cout<<endl<<"Test RealFunctionNonAut (expected output=7.1)"<<endl;
  RealFunctionNonAut<double> rf2(real_func_non_aut);
  Point<double> x2(2);
  x2[0]=1;
  x2[1]=2;
  double rv3;
  rf2.eval(4.1,x2,rv3);
  cout << rv3 << endl;

  //test VF2dFrom2ndOrder
  cout << endl << "Test VF2dFrom2ndOrder (expected output=-2 -3.1)"<<endl;
  VF2dFrom2ndOrder<double> vf2(&rf2);
  Point<double> x3(2);
  x3[0]=-1;
  x3[1]=-2;
  Point<double> rv4(2);
  vf2.eval(-.1,x3,rv4);
  cout << rv4[0]<<" "<<rv4[1]<<endl;
  

}
  
