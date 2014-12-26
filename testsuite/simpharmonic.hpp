//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef SIMPHARMONIC_HPP
#define SIMPHARMONIC_HPP

#include <vector>
#include <jlt/matrix.hpp>


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

  void operator()(Real, const std::vector<Real>& y, std::vector<Real>& y_dot)
    {
      y_dot[q] = y[p];
      y_dot[p] = -omega2*y[q];
    }

  std::vector<Real> Exact(Real t, const std::vector<Real>& yinit) const
    {
      using std::sin;
      using std::cos;
      std::vector<Real> yexact(n);

      yexact[q] = yinit[p]/omega*sin(omega*t) + yinit[q]*cos(omega*t);
      yexact[p] = yinit[p]*cos(omega*t) - yinit[q]*omega*sin(omega*t);

      return yexact;
    }

  // Jacobian matrix of the equation.
  void Jacobian(Real, const std::vector<Real>&, const std::vector<Real>&,
		const Real scale, jlt::matrix<Real>& Jac)
    {
      Jac(q,q) = 0.;
      Jac(q,p) = scale;
      Jac(p,q) = -scale*omega2;
      Jac(p,p) = 0.;
    }

  int size() const { return n; }
};

#endif // SIMPHARMONIC_HPP
