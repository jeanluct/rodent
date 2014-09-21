//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <rodent/explicitrk.hpp>
#include <jlt/stlio.hpp>
#include "simpharmonic.hpp"

using namespace rodent;

int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;

  SimpleHarmonic<double> sho(1.);
  std::vector<double> y(sho.size());

  double ark4_acc = 1.e-10;
  double euler_acc = 1.e-6;
  double rk4_step = 0.01;
  double euler_step = 0.0001;

  double t = 100.0;

  cout << "Integrate the simple harmonic oscillator\n";
  cout << "from t = " << 0. << " to " << t << endl;
  cout << "using several integration methods.\n\n";

  // Initial conditions
  y[0] = 0.;
  y[1] = 1.;

  cout << "Initial conditions: " << y << endl;

  AdaptiveRK4<SimpleHarmonic<double> > sho_ark4(sho);
  sho_ark4
    .tolerance(ark4_acc)
    .setState(0.0,y);

  FixedRK4<SimpleHarmonic<double> > sho_frk4(sho);
  sho_frk4
    .stepSize(rk4_step)
    .setState(0.0,y);

  AdaptiveEuler<SimpleHarmonic<double> > sho_aeuler(sho);
  sho_aeuler
    .tolerance(euler_acc)
    .setState(0.0,y);

  FixedEuler<SimpleHarmonic<double> > sho_feuler(sho);
  sho_feuler
    .stepSize(euler_step)
    .setState(0.0,y);

  cout.precision(10);

  cout << "\nAdaptive RK4:  accuracy = " << ark4_acc << endl;
  cout << t << "\t" << sho_ark4(t) << endl;

  cout << "\nFixed RK4:  step size = " << rk4_step << endl;
  cout << t << "\t" << sho_frk4(t) << endl;

  cout << "\nAdaptive Euler:  accuracy = " << euler_acc << endl;
  cout << t << "\t" << sho_aeuler(t) << endl;

  cout << "\nFixed Euler:  step size = " << euler_step << endl;
  cout << t << "\t" << sho_feuler(t) << endl;

  cout << "\nExact:\n";
  cout << t << "\t" << sho.Exact(t,y) << endl;
}
