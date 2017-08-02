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

void d_real_func(const Point<double>& x,double& rv) {
  rv=2*x[0];
}


void real_func2(const Point<double>& x,double& rv) {
  rv=x[0]*x[0]*x[0]+2*x[0]*x[0]-3*x[0]-1;
}

void vf_func(const Point<double>& x,Point<double>& rv) {
  rv[0]=x[0]*x[0]+x[1];
  rv[1]=x[1]+2;
}

void inv_d_vf_func(const Point<double>& x,const Point<double>& dvf,Point<double>& dx) {
  dx[0]=(dvf[0]-dvf[1])/(2*x[0]);
  dx[1]=dvf[1];
}


int main() {

  cout.precision(17);

  cout << "Testing BisectionMethod (expected output=sqrt(2))" << endl;
  RealFunction1d<double> rf(real_func);
  BisectionMethod<double> bm(&rf);
  double a=1;
  double b=2;
  double root=0;
  bm.solve(a,b,root);
  cout << root << endl;

  cout << endl << "Testing FalsePositionMethod (expected output=1.1987)" << endl;
  RealFunction1d<double> rf2(real_func2);
  FalsePositionMethod<double> fpm(&rf2,.0001);
  double a2=1;
  double b2=2;
  double root2=0;
  fpm.solve(a2,b2,root2);
  cout << root2 << endl;
  
  
  cout << endl << "Testing DifferentiableRealFunction1d (expected output=-1 2)" << endl;
  DifferentiableRealFunction1d<double> rf3(real_func,d_real_func);
  Point<double> x1(1,1);
  double rv1;
  rf3(x1,rv1);
  cout << rv1 << endl;
  rf3.derivative(x1,rv1);
  cout << rv1 << endl;


  //1d Newton Method
  cout << endl << "Testing NewtonMethod1d (expected output=sqrt(2)" << endl;
  DifferentiableRealFunction1d<double> rf4(real_func,d_real_func);
  NewtonMethod1d<double> nm(&rf3);
  double root3=0;
  nm.solve(4,root3);
  cout << root3 << endl;


  //multi-variable Newton Method
  cout << endl << "Testing NewtonMethod (expected output=sqrt(2),-2)" << endl;
  VectorField<double> vf2(2,vf_func);
  DerivativeOfVectorField<double> inverse_dvf2(2,inv_d_vf_func);
  NewtonMethod<double> nm2(2,&vf2,&inverse_dvf2);
  Point<double> x0(2);
  Point<double> root4(2);
  x0[0]=1;
  x0[1]=-1;
  root4[0]=0;
  root4[1]=0;
  nm2.solve(x0,root4);
  cout << root4 <<endl;

}
  
