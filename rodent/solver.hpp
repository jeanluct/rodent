//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

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

  static const magT tiny;
  static const magT expand0;

  magT expand(magT err)
    {
      return (err > 0 ?
	      safety_factor*jlt::Exp(expand_factor*jlt::Log(err)) :
	      expand0);
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
      err_rel			( _n ),
      err_abs			( _n ),
      order			( _order ),
      expand_factor		( -1./(magT)(_order+1) ),
      shrink_factor		( -1./(magT)(_order) ),
      safety_factor		( 0.9 ),
      yscal			( _n ),
      yerr			( _n )
    {
      assert(order > 0);
      tolerance(1.0e-6);
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
      absoluteTolerance(_err);
      relativeTolerance(_err);
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
      std::fill(err_abs.begin(),err_abs.end(),_err_abs);
      return Method();
    }

  T_Method& relativeTolerance(magT _err_rel)
    {
      std::fill(err_rel.begin(),err_rel.end(),_err_rel);
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

      //
      // The error weight vector
      //
      // This is an important part of the integration.  It contains
      // an absolute part and a relative part.  In case the absolute
      // part is set to 0, we add a small number "tiny" to the error
      // scale.
      //
      // We also add |yp[i]*dx| to the relative error estimate.  To
      // quote from Numerical Recipes in C, Second Edition, p. 718:

      //   "Here is a more technical point. We have to consider one
      //    additional possibility for yscal. The error criteria
      //    mentioned thus far are “local,” in that they bound the
      //    error of each step individually. In some applications
      //    you may be unusually sensitive about a “global”
      //    accumulation of errors, from beginning to end of the
      //    integration and in the worst possible case where the
      //    errors all are presumed to add with the same sign. Then,
      //    the smaller the stepsize h, the smaller the value that
      //    you will need to impose. Why?  Because there will be
      //    more steps between your starting and ending values of
      //    x. In such cases you will want to set yscal proportional
      //    to h, typically to something like
      //
      //       Delta0 = eps h dydx[i]
      //
      //    This enforces fractional accuracy not on the values of y
      //    but (much more stringently) the increments to those
      //    values at each step. But now look back at (16.2.7) [new
      //    stepsize estimate, based on exapnsion/contraction of
      //    0.2]. If an implicit scaling with h, then the exponent :
      //    is no longer correct: When the stepsize is reduced from
      //    a too-large value, the new predicted value h1 will fail
      //    to meet the desired accuracy when yscal is also altered
      //    to this new h1 value. Instead of 0.20 = 1/5, we must
      //    scale by the exponent 0.25=1/4 for things to work out."
      //
      // Hence, in rodent we set shrink_factor = -1/(order+1) [1/5
      // in the above discussion, since they are talking about RK4],
      // and expand_factor to -1/(order) [1/4 above].
      //

      /* This is a bit of a pain to write in terms of iterators. */
      for (int i = 0; i < n; i++) {
	yscal(i) = err_abs(i) + err_rel(i)*(vecT_traits::mag(y(i))
					  + vecT_traits::mag(yp(i)*dx)) + tiny;
      }

      for (;;) {
	JLT_TRY
	{
	  Method().OneStep(h, y1);
	  errmax = 0.;
	  // This is the L_infinity norm max_i |v_i|.  Could replace
	  // this by an arbitrary norm, provided by traits.
	  for (int i = 0; i < n; i++) {
	    errmax = std::max(errmax, vecT_traits::mag(yerr(i)/yscal(i)));
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
	}
#ifdef __EXCEPTIONS
	catch(std::range_error& oor) {
	  // We've gone outside the domain.
	  // Try to decrease the stepsize.

	  // The solver needs to think the integrator failed
	  // the step, so make errmax larger than 1.
	  errmax = 2;
#ifdef RODENT_DEBUG
	  std::cerr << oor.what();
	  std::cerr << " Reducing stepsize." << std::endl;
#endif
	}
#endif
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
	  JLT_THROW(jlt::stepsize_too_small<magT>
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

// Must initialise (define) tiny and expand0 outside of the class.
// Required by the standard for non-integral types, for some reason.
template<class T_Method, class vecT, class vecT_traits>
const typename AdaptiveSolver<T_Method,vecT,vecT_traits>::magT
AdaptiveSolver<T_Method,vecT,vecT_traits>::tiny = 1.0e-30;

template<class T_Method, class vecT, class vecT_traits>
const typename AdaptiveSolver<T_Method,vecT,vecT_traits>::magT
AdaptiveSolver<T_Method,vecT,vecT_traits>::expand0 = 100;


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
      for (int i = 0; i < n; ++i) yp(i) = y1p(i);

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

  static const magT tiny;
  static const magT expand0;

  magT expand(magT err)
    {
      return (err > 0 ?
	      safety_factor*jlt::Exp(expand_factor*jlt::Log(err)) :
	      expand0);
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
	Think about the error computation more carefully. Compare to
	CVODE.  Maybe update to be like AdaptiveSolver.  Can these two
	classes be merged?  They share a lot of code...

	For instance, implement absolute/relative errors.
      */
      // See also comments in AdaptiveSolver's version of Step.
      // for (int i = 0; i < n; i++)
      // yscal(i) = vecT_traits::mag(y(i)) + vecT_traits::mag(yp(i)*dx) + err;
      for (int i = 0; i < n; i++)
	yscal(i) = vecT_traits::mag(y(i)) + vecT_traits::mag(yp(i)*dx) + tiny;

      for (;;) {
	JLT_TRY
	{
	  Method().OneStep(h, y1, y1p);
	  errmax = 0.;
	  for (int i = 0; i < n; i++) {
	    errmax = std::max(errmax, vecT_traits::mag(yerr(i)/yscal(i)));
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
	}
#ifdef __EXCEPTIONS
	catch(std::range_error& oor) {
	  // We've gone outside the domain.
	  // Try to decrease the stepsize.

	  // The solver needs to think the integrator failed
	  // the step, so make errmax larger than 1.
	  errmax = 2;
#ifdef RODENT_DEBUG
	  std::cerr << oor.what();
	  std::cerr << " Reducing stepsize." << std::endl;
#endif
	}
#endif
#ifdef RODENT_DEBUG
	++n_bad_steps;
#endif
	result = false;
	h *= shrink(errmax);

	// Check if new step size is too small.
	if (vecT_traits::absval(h) < dx_min) {
#ifdef __EXCEPTIONS
	  JLT_THROW(jlt::stepsize_too_small<magT>
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

// Must initialise (define) tiny and expand0 outside of the class.
// Required by the standard for non-integral types, for some reason.
template<class T_Method, class vecT, class vecT_traits>
const typename AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>::magT
AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>::tiny = 1.0e-30;

template<class T_Method, class vecT, class vecT_traits>
const typename AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>::magT
AdaptiveImplicitSolver<T_Method,vecT,vecT_traits>::expand0 = 100;

} // namespace rodent

#endif // RODENT_SOLVER_HPP
