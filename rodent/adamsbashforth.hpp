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
protected:
  using FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::x;
  using FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::dx;
  using FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::y;
  using FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::yp;
  using FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1;				// Derivative at step n-1.

  // A one-step integrator of the correct order, to get us started.
  FixedMidpoint<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth2(T_Func& _f, const stepT x0, const vecT& y0,
		  const stepT dx0)
    : FixedSolver<AdamsBashforth2<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                  (_f.size(), x0, y0, dx0), ypm1(_f.size()),
		  int_startup(_f,x0,y0,dx0), func(_f)
    {
      reset();
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
#ifdef RODENT_ITERATOR_LOOPS
      CIt yit = y.begin(), ypit = yp.begin();
      for (It y1it = y1.begin(), ypm1it = ypm1.begin();
	   y1it != y1.end(); ++y1it, ++yit, ++ypit, ++ypm1it)
        {
	  *y1it = *yit + 0.5L*h*(3*(*ypit) - *ypm1it);
	  
	  // Cycle the data from previous steps.
	  *ypm1it = *ypit;
	}
#else
      for (int i = 0; i < n; ++i) {
	y1[i] = y[i] + 0.5L*h*(3*yp[i] - ypm1[i]);

	// Cycle the data from previous steps.
	ypm1[i] = yp[i];
      }
#endif
    }

  void reset()
    {
      vecT ym1(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for the previous step, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.Restart(x, y, dx);

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
protected:
  using FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::x;
  using FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::dx;
  using FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::y;
  using FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::yp;
  using FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1, ypm2;			// Derivative at steps n-1 and n-2.

  // A one-step integrator of the correct order, to get us started.
  FixedRK4<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth3(T_Func& _f, const stepT x0, const vecT& y0,
		  const stepT dx0)
    : FixedSolver<AdamsBashforth3<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                  (_f.size(), x0, y0, dx0), ypm1(_f.size()), ypm2(_f.size()),
		  int_startup(_f,x0,y0,dx0), func(_f)
    {
      reset();
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
#ifdef RODENT_ITERATOR_LOOPS
      CIt yit = y.begin(), ypit = yp.begin();
      for (It y1it = y1.begin(), ypm1it = ypm1.begin(),
	            ypm2it = ypm2.begin();
	   y1it != y1.end(); ++y1it, ++yit, ++ypit, ++ypm1it, ++ypm2it)
        {
	  *y1it = *yit + (1.L/12.L)*h*(23*(*ypit) - 16*(*ypm1it)
				       + 5*(*ypm2it));
	  
	  // Cycle the data from previous steps.
	  *ypm2it = *ypm1it;
	  *ypm1it = *ypit;
	}
#else
      for (int i = 0; i < n; ++i) {
	y1[i] = y[i] + (1.L/12.L)*h*(23*yp[i] - 16*ypm1[i] + 5*ypm2[i]);

	// Cycle the data from previous steps.
	ypm2[i] = ypm1[i];
	ypm1[i] = yp[i];
      }
#endif
    }

  void reset()
    {
      vecT ym1(n), ym2(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for previous steps, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.Restart(x, y, dx);

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
protected:
  using FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::x;
  using FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::dx;
  using FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::y;
  using FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::yp;
  using FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1, ypm2, ypm3;	// Derivative at steps n-1, n-2, and n-3.

  // A one-step integrator of the correct order, to get us started.
  FixedRK4<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth4(T_Func& _f, const stepT x0, const vecT& y0,
		  const stepT dx0)
    : FixedSolver<AdamsBashforth4<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                  (_f.size(), x0, y0, dx0), ypm1(_f.size()), ypm2(_f.size()),
                  ypm3(_f.size()), int_startup(_f,x0,y0,dx0), func(_f)
    {
      reset();
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
#ifdef RODENT_ITERATOR_LOOPS
      CIt yit = y.begin(), ypit = yp.begin();
      for (It y1it = y1.begin(), ypm1it = ypm1.begin(),
	            ypm2it = ypm2.begin(), ypm3it = ypm3.begin();
	   y1it != y1.end();
	   ++y1it, ++yit, ++ypit, ++ypm1it, ++ypm2it, ++ypm3it)
        {
	  *y1it = *yit + (1.L/24.L)*h*(55*(*ypit) - 59*(*ypm1it) + 37*(*ypm2it)
				       - 9*(*ypm3it));
	  
	  // Cycle the data from previous steps.
	  *ypm3it = *ypm2it;
	  *ypm2it = *ypm1it;
	  *ypm1it = *ypit;
	}
#else
      for (int i = 0; i < n; ++i) {
	y1[i] = y[i] + (1.L/24.L)*h*(55*yp[i] - 59*ypm1[i] + 37*ypm2[i]
				     - 9*ypm3[i]);

	// Cycle the data from previous steps.
	ypm3[i] = ypm2[i];
	ypm2[i] = ypm1[i];
	ypm1[i] = yp[i];
      }
#endif
    }

  void reset()
    {
      vecT ym1(n), ym2(n), ym3(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for previous steps, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.Restart(x, y, dx);

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
protected:
  using FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::x;
  using FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::dx;
  using FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::y;
  using FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::yp;
  using FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT ypm1, ypm2, ypm3, ypm4;	// Derivative at steps n-1 to n-4.

  // A one-step integrator of the correct order, to get us started.
  // RK Cash-Karp is fifth-order, but is adaptive.  Pass a tolerance
  // of 1 to constructor so every step succeeds.
  AdaptiveRKCashKarp<T_Func,vecT,vecT_traits> int_startup;

public:
  AdamsBashforth5(T_Func& _f, const stepT x0, const vecT& y0,
		  const stepT dx0)
    : FixedSolver<AdamsBashforth5<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                  (_f.size(), x0, y0, dx0), ypm1(_f.size()), ypm2(_f.size()),
                  ypm3(_f.size()), ypm4(_f.size()),
		  int_startup(_f,x0,y0,dx0,1.), func(_f)
    {
      reset();
    }

  T_Func& func;

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
#ifdef RODENT_ITERATOR_LOOPS
      CIt yit = y.begin(), ypit = yp.begin();
      for (It y1it = y1.begin(),
	            ypm1it = ypm1.begin(), ypm2it = ypm2.begin(),
	            ypm3it = ypm3.begin(), ypm4it = ypm4.begin();
	   y1it != y1.end();
	   ++y1it, ++yit, ++ypit, ++ypm1it, ++ypm2it, ++ypm3it, ++ypm4it)
        {
	  *y1it = *yit + (1.L/720.L)*h*(1901*(*ypit) - 2774*(*ypm1it)
			    + 2616*(*ypm2it) - 1274*(*ypm3it) + 251*(*ypm4it));
	  
	  // Cycle the data from previous steps.
	  *ypm4it = *ypm3it;
	  *ypm3it = *ypm2it;
	  *ypm2it = *ypm1it;
	  *ypm1it = *ypit;
	}
#else
      for (int i = 0; i < n; ++i) {
	y1[i] = y[i] + (1.L/720.L)*h*(1901*yp[i] - 2774*ypm1[i] + 2616*ypm2[i]
				      - 1274*ypm3[i] + 251*ypm4[i]);

	// Cycle the data from previous steps.
	ypm4[i] = ypm3[i];
	ypm3[i] = ypm2[i];
	ypm2[i] = ypm1[i];
	ypm1[i] = yp[i];
      }
#endif
    }

  void reset()
    {
      vecT ym1(n), ym2(n), ym3(n), ym4(n);

      // Evaluate derivative vector yp at x.
      func(x, y, yp);

      // Integrate back to find data for previous steps, using an
      // explicit one-step method.

      // Restart the startup integrator.
      int_startup.Restart(x, y, dx);

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