#include <iostream>
#include <iomanip>
#include <valarray>
#include <rodent/explicitrk.hpp>
#include <rodent/data.hpp>
#include <jlt/math.hpp>

using namespace rodent;
using namespace std;
using namespace jlt;

class valSimpleHarmonic {
private:
  double omega;
  double omega2;

public:
  enum {q = 0, p, num};

private:
  static const int n = num;

public:
  valSimpleHarmonic(double omega_) : omega(omega_), omega2(omega_*omega_) {}

  void operator()(double, const valarray<double>& y, valarray<double>& y_dot)
    {
      y_dot[q] = y[p];
      y_dot[p] = -omega2*y[q];
    }

  valarray<double> Exact(double t, const valarray<double>& yinit) const
    {
      valarray<double> yexact(n);

      yexact[q] = yinit[p]/omega*Sin(omega*t) + yinit[q]*Cos(omega*t);
      yexact[p] = yinit[p]*Cos(omega*t) - yinit[q]*omega*Sin(omega*t);

      return yexact;
    }

  int size() const { return n; }
};

int main()
{
  valSimpleHarmonic sho(1.);
  valarray<double> y(sho.size());

  double t_last = 10.;

  double t = 100.0, dtsav = 1.0;

  cout << "Integrate the simple harmonic oscillator\n";
  cout << "from t = " << 0. << " to " << t;
  cout << " using adaptive RK4 and the valarray class.\n";
  cout << "Save data points every " << dtsav << " timestep.\n\n";

  // Initial conditions
  y[0] = 0.;
  y[1] = 1.;

  AdaptiveRK4<valSimpleHarmonic,valarray<double> >
    sho_rk(sho, 0.0, y, 0.01, 0, 1.0e-10);

  /*
  cout << t+10 << "\t" << sho_rk(t+10) << endl;
  */

  DataPoints<AdaptiveRK4<valSimpleHarmonic,valarray<double> >,valarray<double> >
    sho_data(sho_rk,0.0,t,dtsav);

  cout.precision(10);

  cout << "Last data points:\n";
  sho_data.PrintOn(cout, t - t_last, t);

  cout << "\nSample 20 more points.\n";
  sho_data.Sample(t + t_last,20);
  sho_data.PrintOn(cout, t, t + t_last);

  cout << "\nExact:\n";
  cout << t + t_last << "\t" << sho.Exact(t + t_last, y) << endl;
}