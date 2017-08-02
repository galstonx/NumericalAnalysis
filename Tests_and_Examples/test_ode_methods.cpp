#include <iostream>
#include "NumericalAnalysis.hpp"
#include "Sizes.hpp"

using namespace NumericalAnalysis;
using namespace std;


void exp_vf(const double& t,const Point<double>& x,Point<double>& ans) {
  ans[0]=x[0];
}

void spring_vf(const double& t,const Point<double>& x,Point<double>& ans) {
  ans[0]=x[1];
  ans[1]=-x[0];
}

void vf_non_aut1(const double& t,const Point<double>& x,Point<double>& rv) {
  rv[0]=t;
}

void F_non_aut(const double& t,const Point<double>& x,double& rv) {
  rv=t-x[0]-2*x[1];
}

void F_non_aut2(const double& t,const Point<double>& x,double& rv) {
  rv=t;
}

void F_for_2nd_order(const double& t,const Point<double>& x,double& rv) {
  rv=-x[0];
}


int main() {
  cout.precision(17);

  

  cout<<endl;
  cout << endl << "Testing RungeKuttaMethod::solveFinalVal for autonomous 1d ODE (expected output=e)" << endl;
  VectorFieldNonAut<double> vf1(1,exp_vf);
  RungeKuttaMethod<double> rkm1(&vf1);
  Point<double> ans1(1);
  Point<double> x1(1);
  x1[0]=1;
  rkm1.solveFinalVal(1000,0,1,x1,ans1);
  cout << ans1 << endl;

  cout << endl << "Testing RungeKuttaMethod::solve for autonomous 1d ODE (expected output=e)" << endl;
  vector<Point<double>> ans2(1001,Point<double>(1));
  rkm1.solve(1000,0,1,x1,ans2);
  cout << ans2[1000] << endl;


  cout<<endl;
  cout << endl << "Testing RungeKuttaMethod::solveFinalVal for autonomous 2d ODE (expected output=0,-1)" << endl;
  VectorFieldNonAut<double> vf2(2,spring_vf);
  RungeKuttaMethod<double> rkm2(&vf2);
  Point<double> ans3(2);
  Point<double> x2(2);
  x2[0]=1;
  x2[1]=0;
  rkm2.solveFinalVal(1000,0,1.5708,x2,ans3);
  cout << ans3 << endl;

  cout << endl << "Testing RungeKuttaMethod::solve for autonomous 2d ODE (expected output=.707,-.707)" << endl;
  vector<Point<double>> ans4(1001,Point<double>(2));
  rkm2.solve(1000,0,.78540,x2,ans4);
  cout << ans4[1000] << endl;

  cout<< endl;
  cout << endl << "Testing RungeKuttaMethod::solveFinalVal for non-autonomous 1d ODE (expected output=1.5)" << endl;
  VectorFieldNonAut<double> vf3(1,vf_non_aut1);
  RungeKuttaMethod<double> rkm3(&vf3);
  Point<double> ans5(1);
  Point<double> x3(1);
  x3[0]=1;
  rkm3.solveFinalVal(1000,0,1,x3,ans5);
  cout << ans5 << endl;

  cout << endl << "Testing RungeKuttaMethod::solve for non-autonomous 1d ODE (expected output=1.5)" << endl;
  vector<Point<double>> ans6(1001,Point<double>(1));
  rkm3.solve(1000,0,1,x3,ans6);
  cout << ans6[1000] << endl;


  cout << endl << "Testing EulerMethod::solveFinalVal for autonomous 2d ODE (expected output=.707, -.707)" << endl;
  EulerMethod<double> em(&vf2);
  Point<double> x4(2);
  Point<double> ans7(2);
  x4[0]=1;
  x4[1]=0;
  em.solveFinalVal(1000,0,.7853,x4,ans7);
  cout << ans7 <<endl;

  cout << endl << "Testing EulerMethod::solve for autonomous 2d ODE (expected output=.707,-.707)" << endl;
  vector<Point<double>> ans8(1001,Point<double>(2));
  em.solve(1000,0,.7853,x4,ans8);
  cout << ans8[1000] << endl;



  cout << endl << "Test RealFunctionFrom2ndOrderODE (expected output=1.207.. 1.575...)" << endl;
  RealFunctionFrom2ndOrderODE<double> rf_shooting(F_non_aut,1000,0,1,1,0);
  double ans9,ans10;
  rf_shooting(Point<double>(1,1),ans9);
  rf_shooting(Point<double>(1,2),ans10);
  cout << ans9 << " " << ans10 << endl;


  cout << endl << "Test ShootingMethodNonlinear1d (expected output=1.207 -.99 -2.28)" << endl;
  ShootingMethodNonlinear1d<double> sm(F_non_aut,100,0,1,1,0);
  double ans11,ans12,ans13;
  sm.getFunctionValue(1,ans11);
  sm.getFunctionValue(-5,ans12);
  sm.solveDirichletBisection(1,-5,ans13);
  cout << ans11 << " " << ans12 << " " << ans13 << endl;



  cout << endl << "Test RealFunctionForEigenvalue (expected output=0 3.14 -.68)" << endl;
  RealFunctionForEigenvalue<double> rf(F_for_2nd_order,100,0.0,3.14159265358979);
  double ans14,ans15,ans16;
  Point<double> lambda(1);
  lambda[0]=0;
  rf(lambda,ans14);
  lambda[0]=1;
  rf(lambda,ans15);
  lambda[0]=-1;
  rf(lambda,ans16);
  cout << ans14 << " " << ans15 << " " << ans16 << endl;


  cout << endl << "Test EigenvalueFor2ndOrder (expected output=-8 -3 0)" << endl;
  EigenvalueFor2ndOrder<double> ev(F_for_2nd_order,100,0,3.14159265358979);
  vector<double> eigs;
  ev.getEigenvalues(100000,-10,.1,eigs);
  for(int i=0;i<eigs.size();++i) {
    cout << eigs[i] << " ";
  }
  cout << endl;


  return 0;

}
