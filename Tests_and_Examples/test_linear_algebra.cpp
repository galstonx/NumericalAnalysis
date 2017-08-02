#include <iostream>
#include <limits>
#include "NumericalAnalysis.hpp"

using namespace std;
using namespace NumericalAnalysis;


int main() {

  cout.precision(17);

  Matrix<double> A(3,5);
  for(int i=0;i<3;++i) {
    for(int j=0;j<5;++j) {
      A.set(i,j,(double) i+j);
    }
  }

  for(int i=0;i<3;++i) {
    for(int j=0;j<5;++j) {
      cout << A.get(i,j) << " ";
    }
    cout << endl;
  }

  cout<<endl;

  A.rowReduce();
  for(int i=0;i<3;++i) {
    for(int j=0;j<5;++j) {
      cout << A.get(i,j) << " ";
    }
    cout << endl;
  }
  cout<<endl;

  A.swapRows(0,1);
  for(int i=0;i<3;++i) {
    for(int j=0;j<5;++j) {
      cout << A.get(i,j) << " ";
    }
    cout << endl;
  }
  cout<<endl;

  A.rowReduce();
  for(int i=0;i<3;++i) {
    for(int j=0;j<5;++j) {
      cout << A.get(i,j) << " ";
    }
    cout << endl;
  }
  cout<<endl;

  Matrix<double> B(2,2);
  B.set(0,0,7);
  B.set(0,1,1);
  B.set(1,0,2.5);
  B.set(1,1,3.25);
  for(int i=0;i<2;++i) {
    for(int j=0;j<2;++j) {
      cout << B.get(i,j) << " ";
    }
    cout << endl;
  }
  cout<<endl;

  B.rowReduce();
  for(int i=0;i<2;++i) {
    for(int j=0;j<2;++j) {
      cout << B.get(i,j) << " ";
    }
    cout << endl;
  }
  cout<<endl;

  Matrix<double> C(1,2);
  C.set(0,0,0);
  C.set(0,1,2);
  C.rowReduce();
  cout << C << endl;

  Matrix<double> D(1,2);
  D.set(0,0,0);
  D.set(0,1,0);
  D.rowReduce();
  cout << D << endl;

  Matrix<double> E(1,1);
  E.set(0,0,0);
  E.rowReduce();
  cout << E << endl;

  // try solving some equations
  Matrix<double> G(2,4);
  G.set(0,0,1);
  G.set(0,1,2);
  G.set(1,0,3);
  G.set(1,1,4);
  G.set(0,2,1);
  G.set(1,2,1);
  G.set(0,3,10);
  G.set(1,3,20);
  Point<double> ans2(3);
  ans2[2]=1;
  G.rowReduce();
  cout << G << endl;
  G.GaussBackSolve(ans2);
  cout<< "ans=" << ans2[0] << " " << ans2[1] << " "<<ans2[2]<<endl;


  return 0;

}
