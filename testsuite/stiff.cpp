//
// Copyright (c) 2004-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <rodent/explicitrk.hpp>
#include <rodent/implicit.hpp>
#include <rodent/data.hpp>
#include <jlt/math.hpp>
#include "stiff.hpp"

using namespace rodent;

int main()
{
  using std::cout;
  using std::endl;
  using jlt::Exp;

  Stiff stiff(.1,0.,-10000.,0.);
  std::vector<double> y(stiff.size());

  double euler_step = 2.e-4;
  double aeul_acc = 1.e-6;
  double grk_acc = 1.e-6;

  double t = 10.0, dtsav = 1.0;

  cout << "Integrate a stiff set of equations ";
  cout << "from t = " << 0. << " to " << 3*t << endl;
  cout << "using an explicit and an implicit method.\n";

  // Initial conditions
  y[0] = 1.;
  y[1] = 1.;

  cout.precision(10);

  FixedEuler<Stiff> feuler(stiff);
  feuler
    .stepSize(euler_step)
    .setState(0.0,y);

  DataPoints<FixedEuler<Stiff> >
    feuler_data(feuler,0.,3*t,dtsav);

  cout << "\nExact:\n";
  cout << 2*t << "\t" << Exp(.1*2*t) << "\t" << Exp(-10000.*2*t) << endl;
  cout << 3*t << "\t" << Exp(.1*3*t) << "\t" << Exp(-10000.*3*t) << endl;

  cout << "\nFixed Explicit Euler:  step size = " << euler_step << endl;
  feuler_data.PrintOn(cout,2*t,3*t);

  FixedImplicitEuler<Stiff> fimplicit_euler(stiff);
  fimplicit_euler
    .stepSize(euler_step)
    .setState(0.0,y);

  DataPoints<FixedImplicitEuler<Stiff> >
    fimplicit_euler_data(fimplicit_euler,0.0,3*t,dtsav);

  cout << "\nFixed Implicit Euler:  step size = " << euler_step << endl;
  fimplicit_euler_data.PrintOn(cout,2*t,3*t);

  AdaptiveImplicitEuler<Stiff> aimplicit_euler(stiff);
  aimplicit_euler
    .tolerance(aeul_acc)
    .setState(0.0,y);

  DataPoints<AdaptiveImplicitEuler<Stiff> >
    aimplicit_euler_data(aimplicit_euler,0.,3*t,dtsav);

  cout << "\nAdaptive Implicit Euler:  accuracy = " << aeul_acc << endl;
    aimplicit_euler_data.PrintOn(cout,2*t,3*t);

  AdaptiveGRK<Stiff> grk(stiff);
  grk
    .stepSize(t)
    .tolerance(grk_acc)
    .setState(0.0,y);

  DataPoints<AdaptiveGRK<Stiff> >
    grk_data(grk,0.,3*t,dtsav);

  cout << "\nAdaptive Generalized Runge-Kutta:  accuracy = ";
  cout << grk_acc << endl;

  grk_data.PrintOn(cout,2*t,3*t);
}
