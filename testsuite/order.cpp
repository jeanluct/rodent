#include <iostream>
#include <iomanip>
#include <vector>
#include <rodent/explicitrk.hpp>
#include <rodent/adamsbashforth.hpp>
#include <rodent/implicit.hpp>
#include <jlt/stlio.hpp>
#include <jlt/math.hpp>

#if (!defined(ORDER_AB) || !defined(ORDER_FRK)) && !defined(ORDER_ARK)
#  define ORDER_ARK
#endif

using namespace rodent;
using namespace std;
using namespace jlt;

typedef long double Real;

ostream& printErrorOn(ostream& strm, Real err, Real h, int order);

long int factorial(int n);

template<class T>
class ExpM {
private:
  T lambda;
  T y0;

  static const int n = 1;
public:
  static const int u = 0;

  ExpM(T lambda_, T y0_) : lambda(lambda_), y0(y0_) {}

  void operator()(T, const vector<T>& y, vector<T>& y_dot)
    {
      y_dot[u] = lambda*(y[u] + y0);
    }

  vector<T> Exact(T t, const vector<T>& yinit) const
    {
      vector<T> yexact(n);

      // Use Expm1 to get a very accurate value for small t.
      // (especially for yinit[u] = 0.)
      yexact[u] = yinit[u]*Exp(lambda*t) + y0*Expm1(lambda*t);

      return yexact;
    }

  // Jacobian matrix of the equation.
  void Jacobian(T, const vector<T>&, const vector<T>&,
		const T scale, matrix<T>& Jac)
    {
      Jac(u,u) = scale*lambda;
    }

  int size() const { return n; }
};

int main()
{
  int ord = 16;
  Real h, y0ex;
  Real lambda = 1., y0 = 1.;
  ExpM<Real> expm(lambda,y0);
  vector<Real> y(expm.size());

  // Initial condition such that y = t + t^2/2 + t^3/6 +...
  y[0] = 0.;

#if defined(ORDER_AB)
  // Adams-Bashforth methods
  AdamsBashforth2<ExpM<Real>,vector<Real> >
    expm_ab2(expm, 0., y, 0.);
  AdamsBashforth3<ExpM<Real>,vector<Real> >
    expm_ab3(expm, 0., y, 0.);
  AdamsBashforth4<ExpM<Real>,vector<Real> >
    expm_ab4(expm, 0., y, 0.);
  AdamsBashforth5<ExpM<Real>,vector<Real> >
    expm_ab5(expm, 0., y, 0.);

#elif defined(ORDER_FRK)
  // Fixed Runge-Kutta methods
  FixedEuler<ExpM<Real>,vector<Real> >
    expm_feul(expm, 0., y, 0.);
  FixedImplicitEuler<ExpM<Real>,vector<Real> >
    expm_fieul(expm, 0., y, 0.);
  FixedMidpoint<ExpM<Real>,vector<Real> >
    expm_fmid(expm, 0., y, 0.);
  FixedRK4<ExpM<Real>,vector<Real> >
    expm_frk4(expm, 0., y, 0.);

#elif defined(ORDER_ARK)
  // Adaptive Runge-Kutta methods
  AdaptiveEuler<ExpM<Real>,vector<Real> >
    expm_aeul(expm, 0., y, 0., 0., 1.);
  AdaptiveImplicitEuler<ExpM<Real>,vector<Real> >
    expm_aieul(expm, 0., y, 0., 0., 1.);
  AdaptiveMidpoint<ExpM<Real>,vector<Real> >
    expm_amid(expm, 0., y, 0., 0., 1.);
  AdaptiveGRK<ExpM<Real>,vector<Real> >
    expm_agrk(expm, 0., y, 0., 0., 1.);
  AdaptiveRK4<ExpM<Real>,vector<Real> >
    expm_ark4(expm, 0., y, 0., 0., 1.);
  AdaptiveRKCashKarp<ExpM<Real>,vector<Real> >
    expm_rkck(expm, 0., y, 0., 0., 1.);
#endif

  cout.setf(ios::scientific);

  // Output the Log10 of the errors vs the Log10 of the stepsize.
  // Normalize to remove the coefficient of the Taylor expansion of
  // the exact solution (*(order+1)!).

  for (int i = 1; i <= ord; ++i) {
    h = Pow(10.,-i/4.);

    y0ex = expm.Exact(h,y)[0];

    cout << Log10(h);

#ifdef ORDER_AB
    expm_ab2.Restart(0.,y,h);
    printErrorOn(cout,expm_ab2(h)[0] - y0ex,h,2);

    expm_ab3.Restart(0.,y,h);
    printErrorOn(cout,expm_ab3(h)[0] - y0ex,h,3);

    expm_ab4.Restart(0.,y,h);
    printErrorOn(cout,expm_ab4(h)[0] - y0ex,h,4);

    expm_ab5.Restart(0.,y,h);
    printErrorOn(cout,expm_ab5(h)[0] - y0ex,h,5);

#elif defined(ORDER_FRK)
    expm_feul.Restart(0.,y,h);
    printErrorOn(cout,expm_feul(h)[0] - y0ex,h,1);

    expm_fieul.Restart(0.,y,h);
    printErrorOn(cout,expm_fieul(h)[0] - y0ex,h,1);

    expm_fmid.Restart(0.,y,h);
    printErrorOn(cout,expm_fmid(h)[0] - y0ex,h,3);

    expm_frk4.Restart(0.,y,h);
    printErrorOn(cout,expm_frk4(h)[0] - y0ex,h,4);

#elif defined(ORDER_ARK)
    expm_aeul.Restart(0.,y,h);
    printErrorOn(cout,expm_aeul(h)[0] - y0ex,h,2);

    expm_aieul.Restart(0.,y,h);
    printErrorOn(cout,expm_aieul(h)[0] - y0ex,h,2);

    expm_amid.Restart(0.,y,h);
    printErrorOn(cout,expm_amid(h)[0] - y0ex,h,4);

    expm_agrk.Restart(0.,y,h);
    printErrorOn(cout,expm_agrk(h)[0] - y0ex,h,4);

    expm_ark4.Restart(0.,y,h);
    printErrorOn(cout,expm_ark4(h)[0] - y0ex,h,5);

    expm_rkck.Restart(0.,y,h);
    printErrorOn(cout,expm_rkck(h)[0] - y0ex,h,5);
#endif

    cout << endl;
  }
}

ostream& printErrorOn(ostream& strm, Real err, Real h, int order)
{
  strm << "\t";
  if (err != 0) {
    // Multiply by (order+1)! to compensate for the terms in the
    // expansion of expm getting smaller.
    strm << Log10(factorial(order+1)*Abs(err));
  } else {
    // If the error is zero to machine precision, just print the exact
    // value of the error, for comparison purposes.
    // Mark this with a "*".
    strm << Log10(Pow(h,order+1)) << "*";
  }

  return strm;
}

long int factorial(int n)
{
  if (n == 1 || n == 0) return 1;

  assert(n > 0);

  return (n*factorial(n-1));
}
