#ifndef CEXP_HPP
#define CEXP_HPP

#include <jlt/math.hpp>

typedef complex<double> dcomplex;

// Complex exponential function.
// Exactly like the simple Harmonic oscillator.

class CExp {
private:
  dcomplex sigma;

  static const int n = 1;

public:
  CExp(dcomplex sigma_) : sigma(sigma_) {}

  void operator()(double, const vector<dcomplex>& y, vector<dcomplex>& y_dot)
    {
      y_dot[0] = sigma*y[0];
    }

  vector<dcomplex> Exact(double t, const vector<dcomplex>& yinit) const
    {
      vector<dcomplex> yexact(n);

      yexact[0] = yinit[0]*exp(sigma*t);

      return yexact;
    }

  int size() const { return n; }
};

#endif // CEXP_HPP
