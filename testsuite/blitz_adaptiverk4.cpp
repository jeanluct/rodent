#include <cmath>
#include <iostream>
#include <iomanip>
#include <rodent/explicitrk.hpp>
#include <rodent/data.hpp>
#include <rodent/blitz_traits.hpp>
#include <jlt/matrix.hpp>
#include <blitz/vector.h>

using namespace rodent;
using namespace blitz;


class BlitzSimpleHarmonic {
public:
  enum {q = 0, p, num};

private:
  double omega;
  double omega2;

  static const int n = num;

public:
  BlitzSimpleHarmonic(double omega_) : omega(omega_), omega2(omega_*omega_) {}

  void operator()(double, const Vector<double>& y, Vector<double>& y_dot)
    {
      y_dot[q] = y[p];
      y_dot[p] = -omega2*y[q];
    }

  Vector<double> Exact(double t, const Vector<double>& yinit) const
    {
      Vector<double> yexact(n);

      yexact[q] = yinit[p]/omega*sin(omega*t) + yinit[q]*cos(omega*t);
      yexact[p] = yinit[p]*cos(omega*t) - yinit[q]*omega*cos(omega*t);

      return yexact;
    }

  // Jacobian matrix of the equation.
  void Jacobian(double, const Vector<double>&, const Vector<double>&,
		const double scale, matrix<double>& Jac)
    {
      Jac(q,q) = 0.;
      Jac(q,p) = scale;
      Jac(p,q) = -scale*omega2;
      Jac(p,p) = 0.;
    }

  int size() const { return n; }
};


typedef
AdaptiveRK4<BlitzSimpleHarmonic, Vector<double> >
Integrator;

typedef
DataPoints<Integrator, Vector<double> >
Data;

int main()
{
  BlitzSimpleHarmonic sho(1.);
  Vector<double> y(sho.size());

  double t_last = 10.;

  double t = 100.0, dtsav = 1.0;

  cout << "Integrate the simple harmonic oscillator\n";
  cout << "from t = " << 0. << " to " << t;
  cout << " using adaptive RK4 and the blitz Vector class.\n";
  cout << "Save data points every " << dtsav << " timestep.\n\n";

  // Initial conditions
  y[0] = 0.;
  y[1] = 1.;

  Integrator sho_rk(sho, 0.0, y, 0.01, 0, 1.0e-10);

  Data sho_data(sho_rk,0.0,t,dtsav);

  cout.precision(8);

  cout << "Last data points:\n";
  sho_data.PrintOn(cout, t - t_last, t);

  cout << "\nSample 20 more points.\n";
  sho_data.Sample(t + t_last,20);
  sho_data.PrintOn(cout, t, t + t_last);

  cout << "\nExact:\n";
  cout << t + t_last << "\t" << sho.Exact(t + t_last, y) << endl;
}