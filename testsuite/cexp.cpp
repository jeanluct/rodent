#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <rodent/explicitrk.hpp>
#include <jlt/stlio.hpp>
#include "cexp.hpp"

using namespace rodent;

int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;

  CExp expi(dcomplex(0.,1.));
  std::vector<dcomplex> y(expi.size());

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

  AdaptiveRK4<CExp,std::vector<dcomplex> >
    expi_ark4(expi, 0.0, y, 0.01, 0, ark4_acc);

  FixedRK4<CExp,std::vector<dcomplex> >
    expi_frk4(expi, 0.0, y, rk4_step);

  AdaptiveEuler<CExp,std::vector<dcomplex> >
    expi_aeuler(expi, 0.0, y, 0.01, 0, euler_acc);

  FixedEuler<CExp,std::vector<dcomplex> >
    expi_feuler(expi, 0.0, y, euler_step);

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
