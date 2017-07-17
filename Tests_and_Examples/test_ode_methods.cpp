#include <iostream>
#include "NumericalAnalysis.hpp"


using namespace NumericalAnalysis;
using namespace std;


void exp_vf(double t,double x,double& ans) {
  ans=x;
}

void spring(double t,const Point<double>& x,Point<double>& ans) {
  ans[0]=x[1];
  ans[1]=-x[0];
}

void F_non_aut(double t,const Point<double>& x,double& rv) {
  rv=t-x[0]-2*x[1];
}

void F_non_aut2(double t,const Point<double>& x,double& rv) {
  rv=-x[0];
}
int main() {
  cout.precision(17);

  VectorField1dNonAut<double> vf1(exp_vf);
  EulerMethod1d<double> nm1(&vf1);
  vector<double> ans(1000001);
  nm1.solve(1000000,0,1,1,ans); // takes about 20 seconds
  cout << ans[1000000] << endl;
  double ans2;
  nm1.solveFinalVal(100000,0.0,1.0,1.0,ans2); // takes about 10 seconds
  cout << ans2 << endl; 

  cout<<endl;
  VectorFieldNonAut<double> vf2(spring);
  vector<Point<double> > ans3(11,Point<double>(2));

  EulerMethod<double> em(2,&vf2);
  Point<double> ic_spring(2);
  ic_spring[0]=1;
  ic_spring[1]=0;
  em.solve(10,0,.7853,ic_spring,ans3);
  cout << ans3[10][0]<<endl;
  Point<double> ans4(2);
  em.solveFinalVal(1000000,0,.785398163,ic_spring,ans4);
  cout << ans4[0] <<endl;
  cout << ans4[1]<<endl;

  RungeKuttaMethod<double> rkm(2,&vf2);
  rkm.solve(10,0,.785398163,ic_spring,ans3);
  cout << ans3[10][0] << endl;
  cout << ans3[10][1] << endl;
  rkm.solveFinalVal(1000,0,.785398163,ic_spring,ans4);
  cout << endl << "RungeKuttaMethod.solveFinalVal:"<<endl;
  cout << ans4[0] <<endl;
  cout << ans4[1]<<endl;

  cout << endl << "RungeKuttaMethod1d"<<endl;
  RungeKuttaMethod1d<double> rkm2(&vf1);
  vector<double> ans5(1001);
  rkm2.solve(1000,0,1,1,ans5);
  cout << ans5[1000] << endl;
  double ans6;
  rkm2.solveFinalVal(2000,0,1,1,ans6); 
  cout << ans6 << endl; 

  cout << endl << "Test RealFunction1dForShooting (expected output=1.207.. 1.575...)" << endl;
  RealFunctionNonAut<double> F(F_non_aut);
  RealFunction1dForShooting<double> rf_shooting(&F,100,0,1,1,0);
  double ans7,ans8;
  rf_shooting.eval(1,ans7);
  rf_shooting.eval(2,ans8);
  cout << ans7 << " " << ans8 << endl;

  cout << endl << "Test ShootingMethodNonlinear1d (expected output=1.207 -.99 -2.28)" << endl;
  ShootingMethodNonlinear1d<double> sm(&F,100,0,1,1,0);
  double ans9,ans10,ans11;
  sm.getFunctionValue(1,ans9);
  sm.getFunctionValue(-5,ans10);
  sm.solveDirichletBisection(1,-5,ans11);
  cout << ans9 << " " << ans10 << " " << ans11 << endl;


  cout << endl << "Test ShootingMethodNonlinear1d again (expected output=-1.79 .782 5.873)" << endl;
  ShootingMethodNonlinear1d<double> sm2(&F,100,0,1,1,3);
  sm2.getFunctionValue(1,ans9);
  sm2.getFunctionValue(8,ans10);
  sm2.solveDirichletBisection(1,8,ans11);
  cout << ans9 << " " << ans10 << " " << ans11 << endl;

  cout << endl << "Test VF2dForEigenvalue (expected output=-3 3.1)" << endl;
  RealFunctionNonAut<double> F2(F_non_aut2);
  VF2dForEigenvalue<double> vf3(&F2);
  vf3.setEigenvalue(4.1);
  Point<double> x(2);
  x[0]=1;
  x[1]=-3;
  Point<double> ans12(2);
  vf3.eval(0,x,ans12);
  cout << ans12[0] << " " << ans12[1] << endl;

  cout << endl << "Test RealFunction1dForEigenvalue (expected output=-.68 0 3.14 11.54)" << endl;
  RealFunction1dForEigenvalue<double> rf(&F2,100,0.0,3.14159265358979);
  double ans13,ans14,ans15,ans16;
  rf.setEigenvalue(-1);
  rf.eval(ans13);
  rf.setEigenvalue(0);
  rf.eval(ans14);
  rf.setEigenvalue(1);
  rf.eval(ans15);
  rf.setEigenvalue(2);
  rf.eval(ans16);
  cout << ans13 << " " << ans14 << " " << ans15 << " " << ans16 << endl;

  cout << endl << "Test EigenvalueFor2ndOrder (expected output=-8 -3 0)" << endl;
  EigenvalueFor2ndOrder<double> ev(&F2,100,0,3.14159265358979);
  vector<double> eigs;
  ev.getEigenvalues(100000,-10,.1,eigs);
  for(int i=0;i<eigs.size();++i) {
    cout << eigs[i] << " ";
  }
  cout << endl;

  return 0;

}
