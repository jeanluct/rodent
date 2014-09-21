//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef RODENT_ADAMSBASHFORTH_HPP
#define RODENT_ADAMSBASHFORTH_HPP

#include <algorithm>

#include <rodent/solver.hpp>
#include <rodent/traits.hpp>
#include <rodent/explicitrk.hpp>

namespace rodent {

// Only fixed stepsize implementations are applicable to
// Adams-Bashforth methods.

// The data for the coefficients is taken from Stoer and Bulirsch,
// "Introduction to Numerical analysis," Second edition, p. 456.

//
// Second-order Adams-Bashforth
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdamsBashforth2
  : public FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>,
                       vecT, vecT_traits>
{
  typedef
  FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
  FSolver;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1;				// Derivative at step n-1.

  // A one-step integrator of the correct order, to get us started.
  FixedMidpoint<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth2(T_Func& _f)
    :
    FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                ( _f.size() ),
    ypm1			( _f.size() ),
    int_startup			( _f ),
    func			( _f )
    {
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      for (int i = 0; i < n; ++i)
	{
	  y1[i] = y[i] + 0.5L*h*(3*yp[i] - ypm1[i]);

	  // Cycle the data from previous steps.
	  ypm1[i] = yp[i];
	}
    }

  void reset()
    {
      vecT ym1(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for the previous step, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.stepSize(dx).setState(x, y);

      // Step n-1
      int_startup(x-dx,ym1);
      func(x-dx, ym1, ypm1);
    }

}; // class AdamsBashforth2


//
// Third-order Adams-Bashforth
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdamsBashforth3
  : public FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>,
                       vecT, vecT_traits>
{
  typedef
  FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
  FSolver;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1, ypm2;			// Derivative at steps n-1 and n-2.

  // A one-step integrator of the correct order, to get us started.
  FixedRK4<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth3(T_Func& _f)
    :
      FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                ( _f.size() ),
      ypm1			( _f.size() ),
      ypm2			( _f.size() ),
      int_startup		( _f ),
      func			( _f )
    {
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      for (int i = 0; i < n; ++i)
	{
	  y1[i] = y[i] + (1.L/12.L)*h*(23*yp[i] - 16*ypm1[i] + 5*ypm2[i]);

	  // Cycle the data from previous steps.
	  ypm2[i] = ypm1[i];
	  ypm1[i] = yp[i];
	}
    }

  void reset()
    {
      vecT ym1(n), ym2(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for previous steps, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.stepSize(dx).setState(x, y);

      // Step n-1
      int_startup(x-dx,ym1);
      func(x-dx, ym1, ypm1);

      // Step n-2
      int_startup(x-2*dx,ym2);
      func(x-2*dx, ym2, ypm2);
    }

}; // class AdamsBashforth3


//
// Fourth-order Adams-Bashforth
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdamsBashforth4
  : public FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>,
                       vecT, vecT_traits>
{
  typedef
  FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
  FSolver;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1, ypm2, ypm3;	// Derivative at steps n-1, n-2, and n-3.

  // A one-step integrator of the correct order, to get us started.
  FixedRK4<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth4(T_Func& _f )
    :
      FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                ( _f.size() ),
      ypm1			( _f.size() ),
      ypm2			( _f.size() ),
      ypm3			( _f.size() ),
      int_startup		( _f ),
      func			( _f)
    {
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      for (int i = 0; i < n; ++i)
	{
	  y1[i] = y[i] + (1.L/24.L)*h*(55*yp[i] - 59*ypm1[i] + 37*ypm2[i]
				       - 9*ypm3[i]);

	  // Cycle the data from previous steps.
	  ypm3[i] = ypm2[i];
	  ypm2[i] = ypm1[i];
	  ypm1[i] = yp[i];
	}
    }

  void reset()
    {
      vecT ym1(n), ym2(n), ym3(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for previous steps, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.stepSize(dx).setState(x, y);

      // Step n-1
      int_startup(x-dx,ym1);
      func(x-dx, ym1, ypm1);

      // Step n-2
      int_startup(x-2*dx,ym2);
      func(x-2*dx, ym2, ypm2);

      // Step n-3
      int_startup(x-3*dx,ym3);
      func(x-3*dx, ym3, ypm3);
    }

}; // class AdamsBashforth4


//
// Fifth-order Adams-Bashforth
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdamsBashforth5
  : public FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>,
                       vecT, vecT_traits>
{
  typedef
  FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
  FSolver;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1, ypm2, ypm3, ypm4;	// Derivative at steps n-1 to n-4.

  // A one-step integrator of the correct order, to get us started.
  // RK Cash-Karp is fifth-order, but is adaptive.  Pass a tolerance
  // of 1 to constructor so every step succeeds.
  AdaptiveRKCashKarp<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth5(T_Func& _f)
    :
      FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                ( _f.size() ),
      ypm1			( _f.size() ),
      ypm2			( _f.size() ),
      ypm3			( _f.size() ),
      ypm4			( _f.size() ),
      int_startup		( _f ),
      func			( _f )
    {
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      for (int i = 0; i < n; ++i)
	{
	  y1[i] = y[i] + (1.L/720.L)*h*(  1901 * yp[i]
					- 2774 * ypm1[i]
					+ 2616 * ypm2[i]
					- 1274 * ypm3[i]
					+  251 * ypm4[i]);

	  // Cycle the data from previous steps.
	  ypm4[i] = ypm3[i];
	  ypm3[i] = ypm2[i];
	  ypm2[i] = ypm1[i];
	  ypm1[i] = yp[i];
	}
    }

  void reset()
    {
      vecT ym1(n), ym2(n), ym3(n), ym4(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for previous steps, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.stepSize(dx).setState(x, y);

      // Step n-1
      int_startup(x-dx,ym1);
      func(x-dx, ym1, ypm1);

      // Step n-2
      int_startup(x-2*dx,ym2);
      func(x-2*dx, ym2, ypm2);

      // Step n-3
      int_startup(x-3*dx,ym3);
      func(x-3*dx, ym3, ypm3);

      // Step n-4
      int_startup(x-4*dx,ym4);
      func(x-4*dx, ym4, ypm4);
    }

}; // class AdamsBashforth5


} // namespace rodent

#endif // RODENT_ADAMSBASHFORTH_HPP
