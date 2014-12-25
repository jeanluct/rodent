//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <rodent/explicitrk.hpp>
#include <rodent/data.hpp>
#include "cexp.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using std::vector;
  using jlt::operator<<;
  using rodent::AdaptiveRK4;

  CExp expi(dcomplex(0.,1.));
  vector<dcomplex> y(expi.size());

  // Initial conditions
  y[0] = dcomplex(0.,1.);

  double t_last = 10.;

  double t = 100.0, dtsav = 1.0;

  cout << "Integrate the complex exponential\n";
  cout << "from t = " << 0. << " to " << t;
  cout << " using adaptive RK4.\n";
  cout << "Save data points every " << dtsav << " timestep.\n\n";

  // Initial condition
  y[0] = dcomplex(0.,1.);

  AdaptiveRK4<CExp,vector<dcomplex> > expi_rk(expi);
  expi_rk
    .tolerance(1.0e-10)
    .setState(0.0,y);

  rodent::DataPoints<AdaptiveRK4<CExp,vector<dcomplex> >,vector<dcomplex> >
    expi_data(expi_rk,0.0,t,dtsav);

  cout.precision(10);

  cout << "Last data points:\n";
  expi_data.PrintOn(cout, t - t_last, t);

  cout << "\nSample 20 more points.\n";
  expi_data.Sample(t + t_last,20);
  expi_data.PrintOn(cout, t, t + t_last);

  cout << "\nExact:\n";
  cout << t + t_last << "\t" << expi.Exact(t + t_last, y) << endl;
}
