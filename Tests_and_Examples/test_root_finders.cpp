#include <iostream>
#include "NumericalAnalysis.hpp"
using namespace std;
using namespace NumericalAnalysis;


void exp_vf(const Point<double>& x,Point<double>& rv) {
  rv=x;
}

void real_func(double x,double& rv) {
  rv=x*x-2;
}

void d_real_func(double x,double& rv) {
  rv=2*x;
}


void real_func2(double x,double& rv) {
  rv=x*x*x+2*x*x-3*x-1;
}

void vf_func(const Point<double>& x,Point<double>& rv) {
  rv[0]=x[0]*x[0]+x[1];
  rv[1]=x[1]+2;
}

void inv_d_vf_func(const Point<double>& x,Point<double>& dx,const Point<double>& dvf) {
  dx[0]=(dvf[0]-dvf[1])/(2*x[0]);
  dx[1]=dvf[1];
}


int main() {
  Point<double> v(5);
  v[2]=3;
  cout << v[2] << endl;
  VectorField<double> vf(1,exp_vf);
  Point<double> p(1);
  p[0]=-4;
  Point<double> rv(1);
  vf.eval(p,rv);
  cout<< rv[0]<<endl;

  RealFunction1d<double> rf(real_func);
  BisectionMethod<double> bm(rf);
  double a=1;
  double b=2;
  double root=0;
  bm.solve(a,b,root);
  cout << root << endl;

  RealFunction1d<double> rf2(real_func2);
  FalsePositionMethod<double> fpm(rf2,.0001);
  double a2=1;
  double b2=2;
  double root2=0;
  fpm.solve(a2,b2,root2);
  cout << root2 << endl;
  
  //1d Newton Method
  cout << endl;
  RealFunction1d<double> rf3(real_func,d_real_func);
  NewtonMethod1d<double> nm(rf3);
  double root3=0;
  nm.solve(4,root3);
  cout << root3 << endl;

  //multi-variable Newton Method
  VectorField<double> vf2(2,vf_func);
  vf2.set_inv_d_vf(inv_d_vf_func);
  NewtonMethod<double> nm2(2,vf2);
  Point<double> x0(2);
  Point<double> root4(2);
  x0[0]=1;
  x0[1]=-1;
  root4[0]=0;
  root4[1]=0;
  nm2.solve(x0,root4);
  cout << endl<< root4[0] << " "<<root4[1] <<endl;

}
  
