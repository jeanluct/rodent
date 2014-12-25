//
// Copyright (c) 2004-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <cstdlib>
#include <jlt/mathvector.hpp>
#include <jlt/math.hpp>
#include <jlt/finitediff.hpp>

using namespace jlt;

typedef double Real;

int main()
{
  using std::cout;
  using std::endl;

  const int N = 10;
  Real dx = (Real)1./N;
  Real fuzz = .1;

  mathvector<Real> x(N), y(N);
  mathvector<Real> yp1(N), yp2(N), yp4(N), ypexact(N);

  for (int i = 0; i < N; ++i) {
    x[i] = (i + fuzz*(rand()/RAND_MAX - 0.5))*dx;
    y[i] = Sin(2*M_PI*x[i]);
    ypexact[i] = 2*M_PI*Cos(2*M_PI*x[i]);
  }

  finitediff1(x,y,yp1);
  finitediff2(x,y,yp2);
  finitediff4(x,y,yp4);

  for (int i = 0; i < N; ++i) {
    yp1[i] = Log10(Abs(yp1[i] - ypexact[i]));
    yp2[i] = Log10(Abs(yp2[i] - ypexact[i]));
    yp4[i] = Log10(Abs(yp4[i] - ypexact[i]));
    cout << x[i] << "\t" << yp1[i] << "\t";
    cout << yp2[i] << "\t" << yp4[i] << endl;
  }
}
