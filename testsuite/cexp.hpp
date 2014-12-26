//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef CEXP_HPP
#define CEXP_HPP


typedef std::complex<double> dcomplex;

// Complex exponential function.
// Exactly like the simple Harmonic oscillator.

class CExp {
private:
  dcomplex sigma;

  static const int n = 1;

public:
  CExp(dcomplex sigma_) : sigma(sigma_) {}

  void operator()(double, const std::vector<dcomplex>& y,
		  std::vector<dcomplex>& y_dot)
    {
      y_dot[0] = sigma*y[0];
    }

  std::vector<dcomplex> Exact(double t,
			      const std::vector<dcomplex>& yinit) const
    {
      std::vector<dcomplex> yexact(n);

      yexact[0] = yinit[0]*exp(sigma*t);

      return yexact;
    }

  int size() const { return n; }
};

#endif // CEXP_HPP
