#ifndef STIFF_HPP
#define STIFF_HPP

#include <vector>
#include <jlt/matrix.hpp>
#include <jlt/math.hpp>


class Stiff {
public:
  enum {u = 0, v, num};

private:
  double uu, uv;
  double vv, vu;
  std::vector<double> y1peps;		// For finite differencing

  static const int n = num;

public:
  Stiff(double uu_, double uv_, double vv_, double vu_)
    : uu(uu_), uv(uv_), vv(vv_), vu(vu_), y1peps(num) {}

  // Good stiff values (NRC 2nd ed., p. 734).
  Stiff() : uu(998.), uv(1998.), vv(-1999.), vu(-999.), y1peps(num) {}

  void operator()(double, const std::vector<double>& y,
		  std::vector<double>& y_dot)
    {
      y_dot[u] = uu*y[u] + uv*y[v];
      y_dot[v] = vu*y[u] + vv*y[v];
    }

  // Solution corresponding to the default values
  // and initial conditions u = 1, v = 0.
  std::vector<double> Exact(double t, const std::vector<double>& yinit) const
    {
      using jlt::Exp;
      std::vector<double> yexact(n);

      yexact[u] = 2.*Exp(-t) - Exp(-1000.*t);
      yexact[v] = -Exp(-t) + Exp(-1000.*t);

      return yexact;
    }

  // Jacobian matrix of the equation.

  // Exact
  void Jacobian(double, const std::vector<double>&, const std::vector<double>&,
		const double scale, jlt::matrix<double>& Jac)
    {
      Jac(u,u) = scale*uu;
      Jac(u,v) = scale*uv;
      Jac(v,u) = scale*vu;
      Jac(v,v) = scale*vv;
    }

#if 0==1
  // By finite differencing
  void Jacobian(double x1, std::vector<double>& y1,
		const std::vector<double>& y1p,
		const double scale, jlt::matrix<double>& Jac)
    {
	double dh, temp, eps = 1.e-4;

	for (int j = 0; j < num; ++j)
	{
	  temp = y1[j];
	  dh = eps*scale;
	  y1[j] += dh;
	  this->operator()(x1, y1, y1peps);
	  y1[j] = temp;
	  for (int i = 0; i < num; ++i)
	  {
	    Jac(i,j) = (y1peps[i] - y1p[i])/eps;
	  }
	}
    }
#endif

  int size() const { return n; }
};

#endif // STIFF_HPP
