#ifndef SIMPHARMONIC_HPP
#define SIMPHARMONIC_HPP

#include <vector>
#include <jlt/matrix.hpp>
#include <jlt/math.hpp>

using namespace std;
using namespace jlt;

template<class Real>
class SimpleHarmonic {
public:
  enum {q = 0, p, num};

private:
  Real omega;
  Real omega2;

  static const int n = 2;

public:
  SimpleHarmonic(Real omega_) : omega(omega_), omega2(omega_*omega_) {}

  void operator()(Real, const vector<Real>& y, vector<Real>& y_dot)
    {
      y_dot[q] = y[p];
      y_dot[p] = -omega2*y[q];
    }

  vector<Real> Exact(Real t, const vector<Real>& yinit) const
    {
      vector<Real> yexact(n);

      yexact[q] = yinit[p]/omega*Sin(omega*t) + yinit[q]*Cos(omega*t);
      yexact[p] = yinit[p]*Cos(omega*t) - yinit[q]*omega*Sin(omega*t);

      return yexact;
    }

  // Jacobian matrix of the equation.
  void Jacobian(Real, const vector<Real>&, const vector<Real>&,
		const Real scale, matrix<Real>& Jac)
    {
      Jac(q,q) = 0.;
      Jac(q,p) = scale;
      Jac(p,q) = -scale*omega2;
      Jac(p,p) = 0.;
    }

  int size() const { return n; }
};

#endif // SIMPHARMONIC_HPP
