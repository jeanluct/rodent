//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <rodent/explicitrk.hpp>
#include <jlt/stlio.hpp>
#include "cexp.hpp"


int main()
{
  using std::cout;
  using std::endl;
  using std::vector;
  using jlt::operator<<;

  CExp expi(dcomplex(0.,1.));
  vector<dcomplex> y(expi.size());

  double ark4_acc = 1.e-10;
  double euler_acc = 1.e-6;
  double rk4_step = 0.01;
  double euler_step = 0.0001;

  double t = 100.0;

  cout << "Integrate the complex exponential\n";
  cout << "from t = " << 0. << " to " << t << endl;
  cout << "using several integration methods.\n\n";

  // Initial conditions
  y[0] = dcomplex(0.,1.);

  cout << "Initial condition: " << y << endl;

  rodent::AdaptiveRK4<CExp,vector<dcomplex> > expi_ark4(expi);
  expi_ark4
    .tolerance(ark4_acc)
    .setState(0.0,y);

  rodent::FixedRK4<CExp,vector<dcomplex> > expi_frk4(expi);
  expi_frk4
    .stepSize(rk4_step)
    .setState(0.0,y);

  rodent::AdaptiveEuler<CExp,vector<dcomplex> > expi_aeuler(expi);
  expi_aeuler
    .tolerance(euler_acc)
    .setState(0.0,y);

  rodent::FixedEuler<CExp,vector<dcomplex> > expi_feuler(expi);
  expi_feuler
    .stepSize(euler_step)
    .setState(0.0,y);

  cout.precision(10);

  cout << "\nAdaptive RK4:  accuracy = " << ark4_acc << endl;
  cout << t << "\t" << expi_ark4(t) << endl;

  cout << "\nFixed RK4:  step size = " << rk4_step << endl;
  cout << t << "\t" << expi_frk4(t) << endl;

  cout << "\nAdaptive Euler:  accuracy = " << euler_acc << endl;
  cout << t << "\t" << expi_aeuler(t) << endl;

  cout << "\nFixed Euler:  step size = " << euler_step << endl;
  cout << t << "\t" << expi_feuler(t) << endl;

  cout << "\nExact:\n";
  cout << t << "\t" << expi.Exact(t,y) << endl;
}
