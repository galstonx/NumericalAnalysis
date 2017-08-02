#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "NumericalAnalysis.hpp"

using namespace std;
using namespace NumericalAnalysis;

void f_nonhomog(const Point<double>& x,double& rv) {
  rv=4;
}

void g_bc(const Point<double>& x,double& rv) {
  rv=x[0]*x[0]+x[1]*x[1];

}

void f_nonhomog2(const Point<double>& x,double& rv) {
  rv=0;
}

void g_bc2(const Point<double>& x,double& rv) {
  rv=x[0]*x[0]-x[1]*x[1];

}

void f_nonhomog3(const Point<double>& x,double& rv) {
  rv=0-50*sin(10*x[0]);
}

void g_bc3(const Point<double>& x,double& rv) {
  rv=x[0]*x[0]-x[1]*x[1]+.5*sin(10*x[0]);

}


void f_nonhomog4(const Point<double>& x,double& rv) {
  rv=0;
}

void g_bc4(const Point<double>& x,double& rv) {
  rv=sin(10*x[0])-cos(5*x[1]);
}


int main() {

  cout << "Test PoissonBasicMethod2d" << endl;
  RealFunction<double> f(2,f_nonhomog4);
  RealFunction<double> g(2,g_bc4);
  double x_min=-1;
  double x_max=1;
  double x_steps=50;
  double y_min=-1;
  double y_max=1;
  double y_steps=50;
  double z_min=-1.5;
  double z_max=1.5;
  PoissonBasicMethod2d<double> pbm(&f,&g,x_min,x_max,y_min,y_max,x_steps,y_steps);
  Point<double> ans1((x_steps+1)*(y_steps+1));
  pbm.solve(ans1);
  cout << "Writing data to file pde_data1" << endl;
  ofstream pde_data1;
  pde_data1.open("pde_data1");
  int counter=0;
  pde_data1 << x_min << " " << x_max << " " << x_steps << endl;
  pde_data1 << y_min << " " << y_max << " " << y_steps << endl;
  pde_data1 << z_min << " " << z_max << endl;
  for(int i=0;i<=y_steps;++i) {
    for(int j=0;j<=x_steps;++j) {
      pde_data1 << ans1[counter] << " ";
      ++counter;
    }
    pde_data1 << endl;
  }
  pde_data1.close();



  return 0;
}
