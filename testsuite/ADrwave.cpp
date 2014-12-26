//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <cmath>
#include <rodent/exptimediff.hpp>
#include <jlt/matrix.hpp>
#include <jlt/stlio.hpp>
#include "ADrwave.hpp"


int main(int argc, const char **argv)
{
  using std::cout;
  using std::cerr;
  using std::endl;

  int N = 20;			// Number of modes.
  double D = .001;		// Diffusivity
  double U = 1;			// Magnitude of velocity field.
  double S = 1;			// Source strength.
  int Niter = 200;		// Number of iterations.

  const double L = 1, T = 1;
  const double Deff0 = S*L*L/(4*M_PI*M_PI);
  const double tiny = 1e-10;

  // Get an idea if enough modes are included:
  cerr << "Stiffness:    " << -4*M_PI*M_PI * N*N * D*T/(L*L) << endl;
  cerr << "Diff cut-off: ";
  if (U > 0) {
    // Use strain unit of time.
    cerr << exp(-4*M_PI*M_PI * N*N * D/(U*L)) << endl;
  } else {
    // Use L/U as unit of time.
    cerr << exp(-4*M_PI*M_PI * N*N * D*T/(L*L)) << endl;
  }
  ADrwave rw(N,D,U,L,T);

  // Linear part of the equations.
  std::vector<double> c(rw.size());
  for (int row = 0; row < rw.size(); ++row) {
    int m, n, ri;
    rw.ipk(row,m,n,ri);
    c[row] = -(4*M_PI*M_PI*D/(L*L))*(m*m + n*n);
  }

  std::vector<double> y(rw.size());
  jlt::matrix<std::complex<double> > Th(2*N+1,2*N+1);

  // Initial condition.
  y[rw.pk(1,0,0)] = M_SQRT1_2;

#if 0
  /* Broken? */
  double acc = 1.0e-5;		// Accuracy (tolerance) of integration.
  rodent::AdaptiveETDRK4<ADrwave> int_rw(rw, c);
  int_rw
    .tolerance(acc)
    .setState(0.0,y);
#else
  double ss = 1.0e-3;		// Stepsize of integration.
  rodent::FixedETDRK4<ADrwave> int_rw(rw, c);
  int_rw
    .stepSize(ss)
    .setState(0.0,y);
#endif

  cout.precision(5);
  cout.setf(std::ios::scientific);

  // Running sum for mean variance.
  double mean_var = 0;
  for (int iter = 1; iter <= Niter; ++iter) {
    double t = iter*0.5;
    // Integrate to the very edge of the interval...
    int_rw(t-tiny,y);
    // ... then jump over it.
    int_rw.setState(t+tiny,y);

    mean_var += rw.variance(y);

    // Effective diffusivity.
    double Deff = Deff0 / std::sqrt(mean_var/iter);
    cout << t << "\t" << rw.variance(y) << "\t";
    cout << Deff << endl;
  }

  return 0;
}
