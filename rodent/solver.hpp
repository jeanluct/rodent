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

  typedef
  SolverBase<FixedSolver<T_Method,vecT,vecT_traits>, vecT, vecT_traits>
  SBase;

public:
  typedef typename vecT_traits::step_type	stepT;

protected:
  using SBase::x;
  using SBase::dx;
  using SBase::dx_min;
  using SBase::dx_max;
  using SBase::y;
  using SBase::yp;
  using SBase::n;

public:
  FixedSolver(const int _n)
    : SolverBase<FixedSolver<T_Method,vecT,vecT_traits>, vecT,
                 vecT_traits>(_n)
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

  typedef
  SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>, vecT, vecT_traits>
  SBase;

protected:
  using SBase::x;
  using SBase::dx;
  using SBase::dx_min;
  using SBase::dx_max;
  using SBase::y;
  using SBase::yp;
  using SBase::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
  typedef typename vecT_traits::vec_mag_type	vecmagT;

private:
  vecmagT err_rel;			// Desired relative accuracy.
  vecmagT err_abs;			// Desired absolute accuracy.

  static const magT tiny = 1.0e-30;

  magT expand(magT err)
    {
      return safety_factor*jlt::Exp(expand_factor*jlt::Log(err));
    }

  magT shrink(magT err)
    {
      return safety_factor*jlt::Exp(shrink_factor*jlt::Log(err));
    }

  const int order;			// The order of the method used.
                                        // Error scale as (order+1).
  magT expand_factor;			// How much to expand a good step.
  magT shrink_factor;			// How much to shrink a bad step.
  magT safety_factor;			// Safety factor for step adjustment.

  vecmagT yscal;			// Scale factor vector.

protected:
  vecT yerr;				// Error vector.

public:
  AdaptiveSolver(const int _n, const int _order)
    :
      SolverBase<AdaptiveSolver<T_Method,vecT,vecT_traits>, vecT, vecT_traits>
                                ( _n ),
      err_rel			( _n, 1.0e-6 ),
      err_abs			( _n, 1.0e-6 ),
      order			( _order ),
      expand_factor		( -1./(magT)(_order+1) ),
      shrink_factor		( -1./(magT)(_order) ),
      safety_factor		( 0.9 ),
      yscal			( _n ),
      yerr			( _n )
    {
      assert(order > 0);
    }


  //
  // Methods for Setting Parameters
  //

  T_Method& tolerance(vecmagT _err)
    {
      err_abs = _err;
      err_rel = _err;
      return Method();
    }

  T_Method& tolerance(magT _err)
    {
      for (int i = 0; i < n; ++i)
	{
	  err_abs[i] = err_rel[i] = _err;
	}
      return Method();
    }

  T_Method& absoluteTolerance(vecmagT _err_abs)
    {
      err_abs = _err_abs;
      return Method();
    }

  T_Method& relativeTolerance(vecmagT _err_rel)
    {
      err_rel = _err_rel;
      return Method();
    }

  T_Method& absoluteTolerance(magT _err_abs)
    {
      for (int i = 0; i < n; ++i)
	{
	  err_abs[i] = _err_abs;
	}
      return Method();
    }

  T_Method& relativeTolerance(magT _err_rel)
    {
      for (int i = 0; i < n; ++i)
	{
	  err_rel[i] = _err_rel;
	}
      return Method();
    }

  T_Method& expandFactor(magT _expand_factor)
    {
      expand_factor = _expand_factor;
      return Method();
    }

  T_Method& shrinkFactor(magT _shrink_factor)
    {
      shrink_factor = _shrink_factor;
      return Method();
    }

  T_Method& safetyFactor(magT _safety_factor)
    {
      safety_factor = _safety_factor;
      return Method();
    }

  //
  // Query Parameters
  //


  //
  // Take a Step
  //

  bool Step(vecT& y1)
    {
      magT errmax;
      stepT h = dx;
      bool result = true;

      for (int i = 0; i < n; i++) {
	yscal[i] = err_abs[i] + err_rel[i]*(vecT_traits::mag(y[i])
					  + vecT_traits::mag(yp[i]*dx)) + tiny;
      }

      for (;;) {
	Method().OneStep(h, y1);
	errmax = 0.;
	// This is the L_infinity norm max_i |v_i|.  Could replace
	// this by an arbitrary norm, provided by traits.
	for (int i = 0; i < n; i++) {
	  errmax = std::max(errmax, vecT_traits::mag(yerr[i]/yscal[i]));
	}

	if (errmax <= 1) {
	  x += h;
#ifdef RODENT_DEBUG
	  ++n_good_steps;
#endif
	  // Update derivative at x + dx.
	  Method().func(x, y1, yp);

	  dx = h*expand(errmax);

	  // Check if new step size is too large.
	  if (dx_max != 0) {
	    if (vecT_traits::absval(dx) > dx_max) {
	      dx = (dx >= 0 ? dx_max : -dx_max);
	    }
	  }

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
	  std::cerr << "h = " << h << std::endl;
	  std::cerr << "dx_min = " << dx_min << std::endl;
	  _THROW(jlt::stepsize_too_small<magT>
		 ("Stepsize too small in rodent::AdaptiveSolver::Step.",h));
#else
	  std::cerr << "Stepsize too small in rodent::AdaptiveSolver::Step.\n";
	  std::exit(1);
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

  typedef
  SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>, vecT, vecT_traits>
  SBase;

protected:
  using SBase::x;
  using SBase::dx;
  using SBase::dx_min;
  using SBase::dx_max;
  using SBase::y;
  using SBase::yp;
  using SBase::n;

public:
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT y1p;				// Derivative at y1.

public:
  FixedImplicitSolver(const int _n)
    :
      SolverBase<FixedImplicitSolver<T_Method,vecT,vecT_traits>, vecT,
                 vecT_traits>   ( _n ),
      y1p			( _n )
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

  typedef
  SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>,
	     vecT, vecT_traits>
  SBase;

protected:
  using SBase::x;
  using SBase::dx;
  using SBase::dx_min;
  using SBase::dx_max;
  using SBase::y;
  using SBase::yp;
  using SBase::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
  typedef typename vecT_traits::vec_mag_type	vecmagT;

private:
  magT err;				// Desired accuracy.

  static const magT tiny = 1.0e-30;

  magT expand(magT err)
    {
      return safety_factor*jlt::Exp(expand_factor*jlt::Log(err));
    }

  magT shrink(magT err)
    {
      return safety_factor*jlt::Exp(shrink_factor*jlt::Log(err));
    }

  const int order;			// The order of the method used.
                                        // Error scale as (order+1).
  magT expand_factor;			// How much to expand a good step.
  magT shrink_factor;			// How much to shrink a bad step.
  magT safety_factor;			// Safety factor for step adjustment.

  vecmagT yscal;			// Scale factor vector.
  vecT y1p;				// Derivative at y1.

protected:
  vecT yerr;				// Error vector.

public:
  AdaptiveImplicitSolver(const int _n, const int _order)
    :
      SolverBase<AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>, vecT,
                 vecT_traits>
                                ( _n ),
      err			( 1.0e-6 ),
      order			( _order ),
      expand_factor		( -1./(magT)(_order+1) ),
      shrink_factor		( -1./(magT)(_order) ),
      safety_factor		( 0.9 ),
      yscal			( _n ),
      y1p			( _n ),
      yerr			( _n )
    {
      assert(order > 0);
    }

  //
  // Methods for Setting Parameters
  //

  T_Method& tolerance(magT _err)
    {
      err = _err;
      return Method();
    }

  T_Method& expandFactor(magT _expand_factor)
    {
      expand_factor = _expand_factor;
      return Method();
    }

  T_Method& shrinkFactor(magT _shrink_factor)
    {
      shrink_factor = _shrink_factor;
      return Method();
    }

  T_Method& safetyFactor(magT _safety_factor)
    {
      safety_factor = _safety_factor;
      return Method();
    }

  //
  // Query Parameters
  //


  //
  // Take a Step
  //

  bool Step(vecT& y1)
    {
      magT errmax;
      stepT h = dx;
      bool result = true;

      /*
	Think about the error computation more carefully. Compare to CVODE.
      */
      // for (int i = 0; i < n; i++)
      // yscal[i] = vecT_traits::mag(y[i]) + vecT_traits::mag(yp[i]*dx) + err;
      for (int i = 0; i < n; i++)
	yscal[i] = vecT_traits::mag(y[i]) + vecT_traits::mag(yp[i]*dx) + tiny;

      for (;;) {
	Method().OneStep(h, y1, y1p);
	errmax = 0.;
	for (int i = 0; i < n; i++) {
	  errmax = std::max(errmax, vecT_traits::mag(yerr[i]/yscal[i]));
	}
	errmax /= err;
	if (errmax <= 1) {
	  x += h;
#ifdef RODENT_DEBUG
	  ++n_good_steps;
#endif
	  // Update derivative at x + dx.
	  Method().func(x, y1, yp);

	  dx = h*expand(errmax);

	  // Check if new step size is too large.
	  if (dx_max != 0) {
	    if (vecT_traits::absval(h) > dx_max) {
	      h = (h >= 0 ? dx_max : -dx_max);
	    }
	  }

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
	  _THROW(jlt::stepsize_too_small<magT>
	    ("Stepsize too small in rodent::AdaptiveImplicitSolver::Step.",h));
#else
	  std::cerr <<
	    "Stepsize too small in rodent::AdaptiveImplicitSolver::Step.\n";
	  std::exit(1);
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
