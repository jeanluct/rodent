//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <rodent/explicitrk.hpp>
#include <rodent/data.hpp>
#include "simpharmonic.hpp"


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
  cout << " using adaptive RK4.\n";
  cout << "Save data points every " << dtsav << " timestep.\n\n";

  // Initial conditions
  y[0] = 0.;
  y[1] = 1.;

  cout.precision(10);

  rodent::AdaptiveRK4<SimpleHarmonic<double> >
    sho_rk(sho);
  sho_rk
    .tolerance(1.0e-10)
    .setState(0.0,y);

  rodent::DataPoints<rodent::AdaptiveRK4<SimpleHarmonic<double> > >
    sho_data(sho_rk,0.0,t,dtsav);

  cout << "Last data points:\n";
  sho_data.PrintOn(cout, t - t_last, t);

  cout << "\nSample 20 more points.\n";
  sho_data.Sample(t + t_last,20);
  sho_data.PrintOn(cout, t, t + t_last);

  cout << "\nExact:\n";
  cout << t + t_last << "\t" << sho.Exact(t + t_last, y) << endl;
}
