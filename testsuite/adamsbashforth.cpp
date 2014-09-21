//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <rodent/adamsbashforth.hpp>
#include <rodent/data.hpp>
#include "simpharmonic.hpp"

using namespace rodent;

int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;

  SimpleHarmonic<double> sho(1.);
  std::vector<double> y(sho.size());

  double t_last = 10.;

  double t = 100.0, dtsav = 1.0;

  cout << "Integrate the simple harmonic oscillator\n";
  cout << "from t = " << 0. << " to " << t;
  cout << " using second-order Adams-Bashforth.\n";
  cout << "Save data points every " << dtsav << " timestep.\n\n";

  // Initial conditions
  y[0] = 0.;
  y[1] = 1.;

  cout.precision(10);

  AdamsBashforth2<SimpleHarmonic<double> > sho_rk(sho);
  sho_rk.setState(0.0,y);

  DataPoints<AdamsBashforth2<SimpleHarmonic<double> > >
    sho_data(sho_rk,0.0,t,dtsav);

  cout << "Last data points:\n";
  sho_data.PrintOn(cout, t - t_last, t);

  cout << "\nExact:\n";
  cout << t << "\t" << sho.Exact(t, y) << endl;
}
