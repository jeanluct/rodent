#include <iostream>
#include <iomanip>
#include <vector>
#include <rodent/explicitrk.hpp>
#include <rodent/data.hpp>
#include "simpharmonic.hpp"

using namespace rodent;

int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;

  SimpleHarmonic<double> sho(1.);
  const int n = sho.size();
  std::vector<double> y(n), err_rel(n), err_abs(n);

  double t_last = 10.;

  double t = 100.0, dtsav = 1.0;

  cout << "Integrate the simple harmonic oscillator\n";
  cout << "from t = " << 0. << " to " << t;
  cout << " using adaptive RK Cash-Karp.\n";
  cout << "Save data points every " << dtsav << " timestep.\n\n";

  // Initial conditions
  y[0] = 0.;
  y[1] = 1.;

  // Specify error as vector.
  err_rel[0] = 1.0e-10;
  err_rel[1] = 1.0e-8;
  err_abs[0] = 1.0e-8;
  err_abs[1] = 1.0e-10;

  AdaptiveRKCashKarp<SimpleHarmonic<double> > sho_rk(sho, 0.0, y);
  sho_rk
    .absoluteTolerance(err_abs)
    .relativeTolerance(err_rel)
    ;

  DataPoints<AdaptiveRKCashKarp<SimpleHarmonic<double> > >
    sho_data(sho_rk,0.0,t,dtsav);

  cout.precision(10);

  cout << "Last data points:\n";
  sho_data.PrintOn(cout, t - t_last, t);

  cout << "\nSample 20 more points.\n";
  sho_data.Sample(t + 10.,20);
  sho_data.PrintOn(cout, t, t + 10.);

  cout << "\nExact:\n";
  cout << t + 10. << "\t" << sho.Exact(t + 10., y) << endl;
}
