#ifndef RODENT_SOLVER_HPP
#define RODENT_SOLVER_HPP

#include <algorithm>
#ifndef __PGI
#  include <cassert>
#else
#  include <assert.h>
#endif
#include <rodent/base.hpp>
#include <rodent/traits.hpp>
#include <jlt/math.hpp>

using namespace jlt;

namespace rodent {

//
// Fixed Step-size Solver
//

template<class T_Method, class vecT, class vecT_traits>
class FixedSolver
  : public SolverBase<FixedSolver<T_Method,vecT,vecT_traits>, vecT, 
                      vecT_traits>
{
private:
  T_Method& Method() { return static_cast<T_Method&>(*this); }

public:
  typedef typename vecT_traits::step_type	stepT;

protected:
  using
  SolverBase<FixedSolver<T_Method,vecT,vecT_traits>,vecT,vecT_traits>::x;
  using
  SolverBase<FixedSolver<T_Method,vecT,vecT_traits>,vecT,vecT_traits>::dx;
  using
  SolverBase<FixedSolver<T_Method,vecT,vecT_traits>,vecT,vecT_traits>::dx_min;
  using
  SolverBase<FixedSolver<T_Method,vecT,vecT_traits>,vecT,vecT_traits>::y;
  using
  SolverBase<FixedSolver<T_Method,vecT,vecT_traits>,vecT,vecT_traits>::yp;
  using
  SolverBase<FixedSolver<T_Method,vecT,vecT_traits>,vecT,vecT_traits>::n;

public:
  FixedSolver(const int _n, const stepT x0, const vecT& y0, const stepT dx0)
    : SolverBase<FixedSolver<T_Method,vecT,vecT_traits>, vecT,
                 vecT_traits>(_n, x0, y0, dx0)
    {
    }

  bool Step(vecT& y1)
    {
      Method().OneStep(dx, y1);

      // Increment x here.
      x += dx;

      // Update derivative at x + dx.
      Method().func(x, y1, yp);

      // Fixed steps always succeed.
#ifdef RODENT_DEBUG
	  ++n_good_steps;
#endif
      return true;
    }

  void reset()
    {
      Method().reset();
    }

}; // class FixedSolver


//
// Adaptive Step-size Solver
//

template<class T_Method, class vecT, class vecT_traits>
class AdaptiveSolver
  : public SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>, vecT,
                      vecT_traits>
{
private:
  T_Method& Method() { return static_cast<T_Method&>(*this); }

protected:
  using SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::x;
  using SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::dx;
  using SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::dx_min;
  using SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::y;
  using SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::yp;
  using SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
  typedef typename vecT_traits::vec_mag_type	vecmagT;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::const_iterator		CIt;
  typedef typename vecT_traits::mag_iterator		Itmag;
  typedef typename vecT_traits::const_mag_iterator	CItmag;
#endif

private:
  vecmagT err_rel;			// Desired relative accuracy.
  vecmagT err_abs;			// Desired absolute accuracy.

  const magT tiny;

  magT expand(magT err)
    {
      return safety_factor*Exp(expand_factor*Log(err));
    }

  magT shrink(magT err)
    {
      return safety_factor*Exp(shrink_factor*Log(err));
    }

  const int order;			// The order of the method used.
                                        // Error scale as (order+1).
  const magT expand_factor;		// How much to expand a good step.
  const magT shrink_factor;		// How much to shrink a bad step.
  const magT safety_factor;		// Safety factor for step adjustment.

  vecmagT yscal;			// Scale factor vector.

protected:
  vecT yerr;				// Error vector.

public:
  AdaptiveSolver(const int _n, const stepT x0, const vecT& y0,
		 const stepT dx0, const magT _dx_min,
		 const vecmagT _err_rel, const vecmagT _err_abs,
		 const int _order)
    : SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>, vecT, vecT_traits>
                (_n, x0, y0, dx0, _dx_min),
		  err_rel(_err_rel), err_abs(_err_abs), tiny(1.e-30),
		  order(_order),
		  expand_factor(-1./(magT)(_order+1)),
		  shrink_factor(-1./(magT)(_order)),
		  safety_factor(0.9),
		  yscal(_n), yerr(_n)
    {
      assert(order > 0);

      // Check if step size is already too small.
      if (vecT_traits::absval(dx0) < dx_min) {
#ifdef __EXCEPTIONS
	_THROW(stepsize_too_small<magT>
	       ("Initial stepsize too small in rodent::AdaptiveSolver.",dx0));
#else
	cerr << "Initial stepsize too small in rodent::AdaptiveSolver.\n";
	exit(1);
#endif
      }
    }

  AdaptiveSolver(const int _n, const stepT x0, const vecT& y0,
		 const stepT dx0, const magT _dx_min,
		 const magT _err_rel, const magT _err_abs,
		 const int _order)
    : SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>, vecT, vecT_traits>
                (_n, x0, y0, dx0, _dx_min),
		  err_rel(_n,_err_rel), err_abs(_n,_err_abs), tiny(1.e-30),
		  order(_order),
		  expand_factor(-1./(magT)(_order+1)),
		  shrink_factor(-1./(magT)(_order)),
		  safety_factor(0.9),
		  yscal(_n), yerr(_n)
    {
      assert(order > 0);

      // Check if step size is already too small.
      if (vecT_traits::absval(dx0) < dx_min) {
#ifdef __EXCEPTIONS
	_THROW(stepsize_too_small<magT>
	       ("Initial stepsize too small in rodent::AdaptiveSolver.",dx0));
#else
	cerr << "Initial stepsize too small in rodent::AdaptiveSolver.\n";
	exit(1);
#endif
      }
    }

  bool Step(vecT& y1)
    {
      magT errmax;
      stepT h = dx;
      bool result = true;

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; i++) {
	yscal[i] = err_abs[i] + err_rel[i]*(vecT_traits::mag(y[i])
					  + vecT_traits::mag(yp[i]*dx)) + tiny;
      }
#else
      {
	CIt yit = y.begin(), ypit = yp.begin();
	CItmag relit = err_rel.begin(), absit = err_abs.begin();
	for (Itmag yscalit = yscal.begin(); yscalit != yscal.end();
	     ++yscalit, ++yit, ++ypit, ++relit, ++absit)
	  {
	    *yscalit = *absit + *relit*(vecT_traits::mag(*yit)
					+ vecT_traits::mag(*ypit*dx)) + tiny;
	  }
      }
#endif

      for (;;) {
	Method().OneStep(h, y1);
	errmax = 0.;
	// This is the L_infinity norm max_i |v_i|.  Could replace
	// this by an arbitrary norm, provided by traits.
#ifndef RODENT_ITERATOR_LOOPS
	for (int i = 0; i < n; i++) {
	  errmax = max(errmax, vecT_traits::mag(yerr[i]/yscal[i]));
	}
#else
	CIt yerrit = yerr.begin();
	for (CItmag yscalit = yscal.begin(); yerrit != yerr.end();
	     ++yerrit, ++yscalit)
	  {
	    errmax = max(errmax, vecT_traits::mag(*yerrit/(*yscalit)));
	  }
#endif
	if (errmax <= 1) {
	  x += h;
#ifdef RODENT_DEBUG
	  ++n_good_steps;
#endif
	  // Update derivative at x + dx.
	  Method().func(x, y1, yp);

	  dx = h*expand(errmax);
	  break;
	}
#ifdef RODENT_DEBUG
	++n_bad_steps;
#endif
	result = false;
	h *= shrink(errmax);

	// Check if new step size is too small.
	if (vecT_traits::absval(h) < dx_min) {
#ifdef __EXCEPTIONS
	  cerr << "h = " << h << endl;
	  cerr << "dx_min = " << dx_min << endl;
	  _THROW(stepsize_too_small<magT>
		 ("Stepsize too small in rodent::AdaptiveSolver::Step.",h));
#else
	  cerr << "Stepsize too small in rodent::AdaptiveSolver::Step.\n";
	  exit(1);
#endif
	}
      }
      return result;
    }

  void reset()
    {
      Method().reset();
    }

}; // AdaptiveSolver


//
// Implicit Solvers
//
//  Some methods, in particular implicit methods, compute the derivative
//  vector at the new dependent variable value as part of their
//  answer.  For those methods, there is no need to call the
//  right-hand side function at the new time-step.  The new derivative
//  is returned by OneStep.
//
//  There exist explicit methods that have this feature, such as some
//  Runge-Kutta-Fehlberg (embedded) methods, notably Fehlberg's
//  original method (Bulirsch and Stoer, p. 452) and the method of
//  Dormand and Prince (Bulirsch and Stoer, p. 454).
//

//
// Fixed Step-size Implicit Solver
//

template<class T_Method, class vecT, class vecT_traits>
class FixedImplicitSolver
  : public SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>, vecT,
                      vecT_traits>
{
private:
  T_Method& Method() { return static_cast<T_Method&>(*this); }

protected:
  using SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::x;
  using SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::dx;
  using SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::dx_min;
  using SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::y;
  using SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::yp;
  using SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT y1p;				// Derivative at y1.

public:
  FixedImplicitSolver(const int _n, const stepT x0, const vecT& y0,
		      const stepT dx0)
    : SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>, vecT,
                 vecT_traits>(_n, x0, y0, dx0), y1p(_n)
    {
    }

  bool Step(vecT& y1)
    {
      Method().OneStep(dx, y1, y1p);

      // Increment x here.
      x += dx;

      // Update derivative at x + dx.
      for (int i = 0; i < n; ++i) yp[i] = y1p[i];

      // Fixed steps always succeed.
#ifdef RODENT_DEBUG
	  ++n_good_steps;
#endif
      return true;
    }

  void reset()
    {
      Method().reset();
    }

}; // class FixedImplicitSolver


//
// Adaptive Step-size Implicit Solver
//

template<class T_Method, class vecT, class vecT_traits>
class AdaptiveImplicitSolver
  : public SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>, vecT,
                      vecT_traits>
{
private:
  T_Method& Method() { return static_cast<T_Method&>(*this); }

protected:
  using SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::x;
  using SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::dx;
  using SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::dx_min;
  using SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::y;
  using SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::yp;
  using SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>,
		   vecT,vecT_traits>::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
  typedef typename vecT_traits::vec_mag_type	vecmagT;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator		It;
  typedef typename vecT_traits::const_iterator		CIt;
#endif

private:
  magT eps;				// Desired accuracy.

  magT expand(magT err)
    {
      return safety_factor*Exp(expand_factor*Log(err));
    }

  magT shrink(magT err)
    {
      return safety_factor*Exp(shrink_factor*Log(err));
    }

  const int order;			// The order of the method used.
                                        // Error scale as (order+1).
  const magT expand_factor;		// How much to expand a good step.
  const magT shrink_factor;		// How much to shrink a bad step.
  const magT safety_factor;		// Safety factor for step adjustment.

  vecmagT yscal;			// Scale factor vector.
  vecT y1p;				// Derivative at y1.

protected:
  vecT yerr;				// Error vector.

public:
  AdaptiveImplicitSolver(const int _n, const stepT x0, const vecT& y0,
			 const stepT dx0, const magT _dx_min,
			 const magT err, const int _order)
    : SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>, vecT,
                 vecT_traits>
                (_n, x0, y0, dx0, _dx_min),
		  eps(err),
		  order(_order),
		  expand_factor(-1./(magT)(_order+1)),
		  shrink_factor(-1./(magT)(_order)),
		  safety_factor(0.9),
		  yscal(_n), y1p(_n), yerr(_n)
    {
      assert(order > 0);

      // Check if step size is already too small.
      if (vecT_traits::absval(dx0) < dx_min) {
#ifdef __EXCEPTIONS
	_THROW(stepsize_too_small<magT>
	("Initial stepsize too small in rodent::ImplicitAdaptiveSolver.",dx0));
#else
	cerr <<
	  "Initial stepsize too small in rodent::ImplicitAdaptiveSolver.\n";
	exit(1);
#endif
      }
    }

  bool Step(vecT& y1)
    {
      magT errmax;
      stepT h = dx;
      bool result = true;

#ifndef RODENT_ITERATOR_LOOPS
      // for (int i = 0; i < n; i++)
      // yscal[i] = vecT_traits::mag(y[i]) + vecT_traits::mag(yp[i]*dx) + eps;
      for (int i = 0; i < n; i++)
	yscal[i] = vecT_traits::mag(y[i]) + vecT_traits::mag(yp[i]*dx) + 1.e-30;
#else
      {
	CIt yit = y.begin(), ypit = yp.begin();
	for (It yscalit = yscal.begin(); yscalit != yscal.end();
	     ++yscalit, ++yit, ++ypit)
	  {
	    //*yscalit = vecT_traits::mag(*yit)
	    //  + vecT_traits::mag(*ypit*dx) + eps;
	    *yscalit = vecT_traits::mag(*yit)
	      + vecT_traits::mag(*ypit*dx) + 1.e-30;
	  }
      }
#endif

      for (;;) {
	Method().OneStep(h, y1, y1p);
	errmax = 0.;
#ifndef RODENT_ITERATOR_LOOPS
	for (int i = 0; i < n; i++) {
	  errmax = max(errmax, vecT_traits::mag(yerr[i]/yscal[i]));
	}
#else
	for (CIt yerrit = yerr.begin(), yscalit = yscal.begin();
	     yerrit != yerr.end(); ++yerrit, ++yscalit)
	  {
	    errmax = max(errmax, vecT_traits::mag(*yerrit/(*yscalit)));
	  }
#endif
	errmax /= eps;
	if (errmax <= 1) {
	  x += h;
#ifdef RODENT_DEBUG
	  ++n_good_steps;
#endif
	  // Update derivative at x + dx.
	  Method().func(x, y1, yp);

	  dx = h*expand(errmax);
	  break;
	}
#ifdef RODENT_DEBUG
	++n_bad_steps;
#endif
	result = false;
	h *= shrink(errmax);

	// Check if new step size is too small.
	if (vecT_traits::absval(h) < dx_min) {
#ifdef __EXCEPTIONS
	  _THROW(stepsize_too_small<magT>
	    ("Stepsize too small in rodent::ImplicitAdaptiveSolver::Step.",h));
#else
	  cerr <<
	    "Stepsize too small in rodent::ImplicitAdaptiveSolver::Step.\n";
	  exit(1);
#endif
	}
      }
      return result;
    }

  void reset()
    {
      Method().reset();
    }

}; // class AdaptiveImplicitSolver

} // namespace rodent

#endif // RODENT_SOLVER_HPP
