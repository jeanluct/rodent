//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <cmath>
#include <jlt/matrix.hpp>
#include <jlt/stlio.hpp>
#include <jlt/math.hpp>


class ADrwave
{
private:
  const int N, NN, N2, NT;
  const double D, U, L, T, T2, ULc;
  static const int r = 0, i = 1;

public:
  int pk(const int m, const int n, const int ri) const;
  void ipk(const int row, int& m, int& n, int& ri) const;

  ADrwave(int N_, double D_, double U_ = 1, double L_ = 1, double T_ = 1) :
    N(N_), NN(2*N+1), N2(2*N*(N+1)), NT(2*N2),
    D(4*M_PI*M_PI*D_/(L_*L_)), U(U_), L(L_), T(T_), T2(0.5*T),
    ULc(M_SQRT2*M_PI*(U/L)) {}

  void operator()(double t, const std::vector<double>& y,
		  std::vector<double>& y_dot)
  {
    static double X1, X2;
    static int it0 = -1;

    double dt = jlt::Mod(t,T);
    int it = (int)(t/T);

    // Generate random angles at first, or if we reach a new interval.
    // This is dangerous if integrating backwards in time.
    if (it != it0 || it0 == -1) {
      X1 = 2*M_PI * ((double)random() / RAND_MAX);
      X2 = 2*M_PI * ((double)random() / RAND_MAX);
      it0 = it;
    }

    // Diffusion term.
    for (int row = 0; row < size(); ++row) {
      int m, n, ri;
      ipk(row,m,n,ri);
      y_dot[row] = -D*(m*m + n*n) * y[row];
    }
    for (int row = 0; row < size(); ++row) y_dot[row] = 0;

    if (std::abs(dt) < T2) {
      // y-dependent wave in x direction.
      double cX1 = std::cos(X1), sX1 = std::sin(X1);
      for (int m = 1; m <= N; ++m) {
	for (int n = -N; n <= N; ++n) {
	  // Real coefficients.
	  {
	    double term1 = 0, term2 = 0;
	    if (n+1 <= N)
	      term1 = cX1*y[pk(m,n+1,r)] - sX1*y[pk(m,n+1,i)];
	    if (n-1 >= -N)
	      term2 = cX1*y[pk(m,n-1,r)] + sX1*y[pk(m,n-1,i)];
	    y_dot[pk(m,n,r)] += -ULc * m * (term1 - term2);
	  }

	  // Imaginary coefficients.
	  {
	    double term1 = 0, term2 = 0;
	    if (n+1 <= N)
	      term1 = cX1*y[pk(m,n+1,i)] + sX1*y[pk(m,n+1,r)];
	    if (n-1 >= -N)
	      term2 = cX1*y[pk(m,n-1,i)] - sX1*y[pk(m,n-1,r)];
	    y_dot[pk(m,n,i)] += -ULc * m * (term1 - term2);
	  }
	}
      }
    } else {
      // x-dependent wave in y direction.
      double cX2 = std::cos(X2), sX2 = std::sin(X2);
      for (int m = 1; m <= N; ++m) {
	for (int n = -N; n <= N; ++n) {
	  if (n != 0) {
	    // Real coefficients.
	    {
	      double term1 = 0, term2 = 0;
	      if (m+1 <= N)
		term1 = cX2*y[pk(m+1,n,r)] - sX2*y[pk(m+1,n,i)];
	      if (m-1 == 0 && n < 0)
		term2 = cX2*y[pk(0,-n,r)] - sX2*y[pk(0,-n,i)];
	      else
		term2 = cX2*y[pk(m-1,n,r)] + sX2*y[pk(m-1,n,i)];
	      y_dot[pk(m,n,r)] += -ULc * n * (term1 - term2);
	    }

	    // Imaginary coefficients.
	    {
	      double term1 = 0, term2 = 0;
	      if (m+1 <= N)
		term1 = cX2*y[pk(m+1,n,i)] + sX2*y[pk(m+1,n,r)];
	      if (m-1 == 0 && n < 0)
		term2 = -cX2*y[pk(0,-n,i)] - sX2*y[pk(0,-n,r)];
	      else
		term2 = cX2*y[pk(m-1,n,i)] - sX2*y[pk(m-1,n,r)];
	      y_dot[pk(m,n,i)] += -ULc * n * (term1 - term2);
	    }
	  }
	}
      }
      // The m=0 coefficients.
      for (int n = 1; n <= N; ++n) {
	// Real coefficients.
	{
	  double term1 = cX2*y[pk(1,n,r)] - sX2*y[pk(1,n,i)];
	  double term2 = cX2*y[pk(1,-n,r)] - sX2*y[pk(1,-n,i)];
	  y_dot[pk(0,n,r)] += -ULc * n * (term1 - term2);
	}

	// Imaginary coefficients.
	{
	  double term1 = cX2*y[pk(1,n,i)] + sX2*y[pk(1,n,r)];
	  double term2 = -cX2*y[pk(1,-n,i)] - sX2*y[pk(1,-n,r)];
	  y_dot[pk(0,n,i)] += -ULc * n * (term1 - term2);
	}
      }
    }
    // Source: s(x) = S sqrt(2) sin(2pi x/L).
    // S should be passed to the class, as should optionally a forcing vector.
    double S = 1;
    y_dot[pk(1,0,i)] += S * (-M_SQRT1_2);
  }

  void Jacobian(double t, const std::vector<double>& y,
		const std::vector<double>& y_dot,
		const double scale, jlt::matrix<double>& Jac)
  {
    // The Jacobian matrix is sparse.  Inefficient for implicit methods.
  }

  void toMode(const std::vector<double> y,
	      jlt::matrix<std::complex<double> >& Th)
  {
    Th(0,0) = 0;

    for (int n = 1; n <= N; ++n) {
      Th(N+0,N+n) = std::complex<double>(y[pk(0,n,r)],y[pk(0,n,i)]);
    }
    for (int m = 1; m <= N; ++m) {
      for (int n = -N; n <= N; ++n) {
	Th(N+m,N+n) = std::complex<double>(y[pk(m,n,r)],y[pk(m,n,i)]);
      }
    }

    for (int n = -N; n <= -1; ++n) {
      Th(N+0,N+n) = std::complex<double>(y[pk(0,-n,r)],-y[pk(0,-n,i)]);
    }
    for (int m = -N; m <= -1; ++m) {
      for (int n = -N; n <= N; ++n) {
	Th(N+m,N+n) = std::complex<double>(y[pk(-m,-n,r)],-y[pk(-m,-n,i)]);
      }
    }
  }

  std::ostream& printModeList(const std::vector<double>& y,
			      std::ostream& strm) const
  {
    for (int n = 1; n <= N; ++n) {
      strm << 0 << "\t" << n << "\t";
      strm << y[pk(0,n,r)] << "\t" << y[pk(0,n,i)] << std::endl;
    }
    for (int m = 1; m <= N; ++m) {
      for (int n = -N; n <= N; ++n) {
	strm << m << "\t" << n << "\t";
	strm << y[pk(m,n,r)] << "\t" << y[pk(m,n,i)] << std::endl;
      }
    }
    return strm;
  }

  double variance(const std::vector<double>& y)
  {
    double var = 0;

    for (int n = 1; n <= N; ++n) {
      var += pow(y[pk(0,n,r)],2) + pow(y[pk(0,n,i)],2);
    }
    for (int m = 1; m <= N; ++m) {
      for (int n = -N; n <= N; ++n) {
	var += pow(y[pk(m,n,r)],2) + pow(y[pk(m,n,i)],2);
      }
    }
    var *= 2;

    return var;
  }

  // Largest mode number.
  int msize() const { return N; }

  // Total size of system.
  int size() const { return NT; }
};


inline int ADrwave::pk(const int m, const int n, const int ri) const
{
  assert(ri == 0 || ri == 1);
  assert(m <= N && n <= N);
  assert((m > 0 && n >= -N) || (m == 0 && n > 0));

  int b = ri*N2;

  if (m == 0) {
    return (b + n-1);
  }
  b += N;
  return (b + NN*(m-1) + N + n);
}

inline void ADrwave::ipk(const int row, int& m, int& n, int& ri) const
{
  assert(row >=0 && row < 2*N2);

  ri = row / N2;
  int r = row % N2;

  if (r < N) {
    m = 0;
    n = r+1;
    return;
  }
  r -= N;

  n = (r % NN) - N;
  m = (r / NN) + 1;
}
