#ifndef RODENT_IMPLICIT_HPP
#define RODENT_IMPLICIT_HPP

#include <iostream>
#include <algorithm>

#include <rodent/solver.hpp>
#include <rodent/traits.hpp>
#include <jlt/exceptions.hpp>
#include <jlt/matrixutil.hpp>


namespace rodent {

//
// Methods
//

template<class T_Func, class vecT, class vecT_traits>
class ImplicitEuler
{
public:
  typedef typename vecT_traits::value_type	value_type;
  typedef typename vecT_traits::step_type	step_type;
  typedef typename vecT_traits::mag_type	mag_type;
  typedef typename vecT_traits::matrix_type	matrix_type;

  ImplicitEuler(T_Func& _f)
    : func(_f), tol(1.e-10), Jac(_f.size(),_f.size()), dy(_f.size())
    {
      idx = new int[func.size()];
    }

  ~ImplicitEuler() { delete[] idx; }

  T_Func& func;

private:
  static const int max_iter = 100;	// Number of Newton iterations
  const mag_type tol;

  matrix_type		Jac;		// -Jacobian
  vecT			dy;		// Solution to J.dy = F
  int			*idx;		// Row permuations for LUdecomp

protected:
  void ieuler_step(const step_type x0, const vecT& y0, const vecT& yp0,
		  const step_type h, vecT& y1, vecT& y1p)
    {
      step_type		x1 = x0 + h;
      mag_type		F;
      int		perm;

      // Initial guess.  Use explicit Euler  (yp0[x] already computed).
      for (int i = 0; i < func.size(); ++i) y1[i] = y0[i] + h*yp0[i];

      // Derivative at current iterate (y1).
      func(x1, y1, y1p);

      for (int ns = 1; ns <= max_iter; ++ns)
      {
	// Jacobian of right-hand side function (scaled by h).
	func.Jacobian(x1, y1, y1p, h, Jac);

	// We need (h Jac - I).
	for (int i = 0; i < func.size(); ++i) Jac(i,i) -= 1.;

	// Compute F(y0), store in dy.
	for (int i = 0; i < func.size(); ++i) dy[i] = y1[i] - y0[i] - h*y1p[i];
	// LU-decompose the Jacobian matrix.
	jlt::LUdecomp<value_type,matrix_type>(Jac, idx, &perm);
	// Solve the equation Jac.dy = -F(y0).
	jlt::LUbacksub<value_type,matrix_type>(Jac, idx, &dy[0]);

	// Compute y1 (new) = y1 (guess) + dy.
	for (int i = 0; i < func.size(); ++i) y1[i] += dy[i];

	// Yet another function call, ready for next step.
	func(x1, y1, y1p);
	// Compute magnitude of F.
	F = 0.;
	for (int i = 0; i < func.size(); ++i)
	  F += vecT_traits::mag(y1[i] - y0[i] - h*y1p[i]);
	// Check to see if convergence was achieved.
	if (F < tol) return;
	// If successful, we have already calculated y1p.  An adaptive
	// solver shouldn't need to compute it again.  Yet we can't
	// overwrite it, since the calling routine might not accept
	// the step.
      }
#ifdef __EXCEPTIONS
      _THROW(jlt::too_many_steps
	("Newton iteration failed to converge in ImplicitEuler ", max_iter));
#else
      std::cerr << "rodent::ImplicitEuler::ieuler_step: ";
      std::cerr << "Failed to converge after " << max_iter << " iterations.\n";
      std::exit(1);
#endif
    }
}; // class ImplicitEuler


//
// Fixed Step Size Implementations
//

template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedImplicitEuler
  : public FixedImplicitSolver<
                       FixedImplicitEuler<T_Func,vecT,vecT_traits>,
                       vecT, vecT_traits>,
    public ImplicitEuler<T_Func,vecT,vecT_traits>
{
protected:
  using FixedImplicitSolver<FixedImplicitEuler<T_Func,vecT,vecT_traits>,
			    vecT, vecT_traits>::x;
  using FixedImplicitSolver<FixedImplicitEuler<T_Func,vecT,vecT_traits>,
			    vecT, vecT_traits>::dx;
  using FixedImplicitSolver<FixedImplicitEuler<T_Func,vecT,vecT_traits>,
			    vecT, vecT_traits>::y;
  using FixedImplicitSolver<FixedImplicitEuler<T_Func,vecT,vecT_traits>,
			    vecT, vecT_traits>::yp;
  using FixedImplicitSolver<FixedImplicitEuler<T_Func,vecT,vecT_traits>,
			    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

  FixedImplicitEuler(T_Func& _f)
    :
      FixedImplicitSolver<FixedImplicitEuler<T_Func,vecT,vecT_traits>,
                          vecT, vecT_traits>	( _f.size() ),
      ImplicitEuler<T_Func,vecT,vecT_traits>	( _f )
    {
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for FixedSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1, vecT& y1p)
    {
      _TRY
      {
	ieuler_step(x, y, yp, h, y1, y1p);
      }
#ifdef __EXCEPTIONS
      catch(jlt::too_many_steps& ex)
      {
	std::cerr << ex.what() << "after " << ex.how_many() << " steps.\n";
	throw;
      }
#endif
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // class FixedImplicitEuler


//
// Adaptive Step Size Implementations
//

//
// First-order Implicit Euler
//
//  Really second-order because doing two half-steps buys an extra
//  order.  However, we only know the error on the 1st order result,
//  so we tell AdaptiveImplicitSolver we are first order (last argument of
//  AdaptiveImplicitSolver constructor is 1).
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveImplicitEuler
  : public AdaptiveImplicitSolver<
                          AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
                          vecT, vecT_traits>,
    public ImplicitEuler<T_Func,vecT,vecT_traits>
{
protected:
  using AdaptiveImplicitSolver<AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
			       vecT, vecT_traits>::x;
  using AdaptiveImplicitSolver<AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
			       vecT, vecT_traits>::dx;
  using AdaptiveImplicitSolver<AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
			       vecT, vecT_traits>::y;
  using AdaptiveImplicitSolver<AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
			       vecT, vecT_traits>::yp;
  using AdaptiveImplicitSolver<AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
			       vecT, vecT_traits>::yerr;
  using AdaptiveImplicitSolver<AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
			       vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;

private:
  vecT y_mid, yp_mid, yh;		// Working variables for OneStep.

  static const int order = 1;

public:
  AdaptiveImplicitEuler(T_Func& _f)
    :
      AdaptiveImplicitSolver<AdaptiveImplicitEuler<T_Func,vecT,vecT_traits>,
                             vecT,vecT_traits>  ( _f.size(), order ),
      ImplicitEuler<T_Func,vecT,vecT_traits>	( _f ),
      y_mid					( n ),
      yp_mid					( n ),
      yh					( n )
    {
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for AdaptiveImplicitSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1, vecT& y1p)
    {
      stepT h_mid = 0.5*h, x_mid = x + h_mid;

      _TRY
      {
	// Two small steps.
	ieuler_step(x, y, yp, h_mid, y_mid, yp_mid);
	// We already have the derivative yp_mid at the midpoint.
	ieuler_step(x_mid, y_mid, yp_mid, h_mid, y1, y1p);	// y1p = y1p_b

	// One large step. (Store derivative in yp_mid.)
	ieuler_step(x, y, yp, h, yh, yp_mid);	// yp_mid = y1p_a

	// Compute error and a 2nd order correction.
	for (int i = 0; i < n; i++) {
	  yerr[i] = y1[i] - yh[i];
	  y1[i] += yerr[i];	// y1_c = 2 y1_b - y1_a

	  // If we take that correction we have to recalculate y1p.
	  // Happily, to the order of the method (2) we can deduce it
	  // without the need for another function call.
	  y1p[i] = 2.*y1p[i] - yp_mid[i];	// y1p_c = 2 y1p_b - y1p_a
	}
      }
#ifdef __EXCEPTIONS
      catch(jlt::too_many_steps& ex)
      {
	std::cerr << ex.what() << "after " << ex.how_many() << " steps.\n";
	throw;
      }
#endif
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // AdaptiveImplicitEuler

} // namespace rodent

#endif // RODENT_IMPLICIT_HPP
