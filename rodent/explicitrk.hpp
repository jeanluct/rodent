#ifndef RODENT_EXPLICITRK_HPP
#define RODENT_EXPLICITRK_HPP

#include <algorithm>

#include <rodent/solver.hpp>
#include <rodent/traits.hpp>
#include <jlt/matrixutil.hpp>


namespace rodent {

//
// Methods
//

//
// Euler step (First-order Runge-Kutta)
//
template<class T_Func, class vecT, class vecT_traits>
class Euler
{
public:
  typedef typename vecT_traits::step_type	step_type;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator	iterator;
  typedef typename vecT_traits::const_iterator	const_iterator;
#endif

  Euler(T_Func& _f) : func(_f) {}

  T_Func& func;

protected:
  void euler_step(const step_type x0, const vecT& y0, const vecT& yp0,
		  const step_type h, vecT& y1) const
    {
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < func.size(); i++) {
	y1[i] = y0[i] + h*yp0[i];
      }
#else
    {
      const_iterator yit = y0.begin(), ypit = yp0.begin();
      for (iterator y1it = y1.begin();
	   y1it != y1.end(); ++y1it, ++yit, ++ypit)
	{
	  *y1it = *yit + *ypit*h;
	}
    }
#endif
    }
}; // class Euler


//
// Midpoint step (Second-order Runge-Kutta)
//
template<class T_Func, class vecT, class vecT_traits>
class Midpoint
{
public:
  typedef typename vecT_traits::step_type	step_type;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator	iterator;
  typedef typename vecT_traits::const_iterator	const_iterator;
#endif

  Midpoint(T_Func& _f) : func(_f), y_2(_f.size()), yp_2(_f.size()) {}

  T_Func& func;

private:
  vecT y_2, yp_2;		// Working variables.

protected:
  void midpoint_step(const step_type x0, const vecT& y0, const vecT& yp0,
		     const step_type h, vecT& y1)
    {
      step_type h_2 = 0.5L*h, x_2 = x0 + h_2;

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < func.size(); i++) { y_2[i] = y0[i] + h_2*yp0[i]; }
#else
      {
	const_iterator y0it = y0.begin(), yp0it = yp0.begin();
	for (iterator y_2it = y_2.begin();
	     y_2it != y_2.end(); ++y_2it, ++y0it, ++yp0it)
	{
	  *y_2it = *y0it + *yp0it*h_2;
	}
      }
#endif
      func(x_2, y_2, yp_2);

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < func.size(); i++) { y1[i] = y0[i] + h*yp_2[i]; }
#else
      {
	const_iterator y0it = y0.begin(), yp_2it = yp_2.begin();
	for (iterator y1it = y1.begin();
	     y1it != y1.end(); ++y1it, ++y0it, ++yp_2it)
	{
	  *y1it = *y0it + *yp_2it*h;
	}
      }
#endif
    }
}; // class Midpoint


//
// Fourth-order Runge-Kutta step
//
template<class T_Func, class vecT, class vecT_traits>
class RK4
{
public:
  typedef typename vecT_traits::step_type	step_type;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator	iterator;
  typedef typename vecT_traits::const_iterator	const_iterator;
#endif

  RK4(T_Func& _f) :
    func(_f), y_2(_f.size()), yp_2(_f.size()), yp1(_f.size()) {}

  T_Func& func;

private:
  vecT y_2, yp_2, yp1;		// Working variables.

protected:
  void rk4_step(const step_type x0, const vecT& y0, const vecT& yp0,
		const step_type h, vecT& y1)
    {
      step_type h_2 = 0.5L*h, h_6 = h/6.L, x_2 = x0 + h_2;

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < func.size(); i++) { y_2[i] = y0[i] + h_2*yp0[i]; }
#else
      {
	const_iterator y0it = y0.begin(), yp0it = yp0.begin();
	for (iterator y_2it = y_2.begin();
	     y_2it != y_2.end(); ++y_2it, ++y0it, ++yp0it)
	{
	  *y_2it = *y0it + *yp0it*h_2;
	}
      }
#endif
      func(x_2, y_2, yp_2);

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < func.size(); i++) { y_2[i] = y0[i] + h_2*yp_2[i]; }
#else
      {
	const_iterator y0it = y0.begin(), yp_2it = yp_2.begin();
	for (iterator y_2it = y_2.begin();
	     y_2it != y_2.end(); ++y_2it, ++y0it, ++yp_2it)
	{
	  *y_2it = *y0it + *yp_2it*h_2;
	}
      }
#endif
      func(x_2, y_2, yp1);

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < func.size(); i++) {
	y_2[i] = y0[i] + h*yp1[i];
	yp1[i] += yp_2[i];
      }
#else
      {
	const_iterator y0it = y0.begin(), yp_2it = yp_2.begin();
	for (iterator y_2it = y_2.begin(), yp1it = yp1.begin();
	     y_2it != y_2.end(); ++y_2it, ++yp1it, ++y0it, ++yp_2it)
	{
	  *y_2it = *y0it + *yp1it*h;
	  *yp1it += *yp_2it;
	}
      }
#endif

      func(x0 + h, y_2, yp_2);
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < func.size(); i++) {
	y1[i] = y0[i] + h_6*(yp0[i] + yp_2[i] + 2.*yp1[i]);
      }
#else
      {
	const_iterator y0it = y0.begin(), yp0it = yp0.begin();
	const_iterator yp_2it = yp_2.begin(), yp1it = yp1.begin();
	for (iterator y1it = y1.begin();
	     y1it != y1.end(); ++y1it, ++y0it, ++yp0it, ++yp_2it, ++yp1it)
	{
	  *y1it = *y0it + h_6*(*yp0it + *yp_2it + 2.*(*yp1it));
	}
      }
#endif
    }
}; // class RK4


//
// Fixed Step Size Implementations
//

// First-order Euler
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedEuler
  : public FixedSolver<FixedEuler<T_Func,vecT,vecT_traits>, vecT, vecT_traits>,
    public Euler<T_Func, vecT, vecT_traits>
{
protected:
  using FixedSolver<FixedEuler<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::x;
  using FixedSolver<FixedEuler<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::dx;
  using FixedSolver<FixedEuler<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::y;
  using FixedSolver<FixedEuler<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::yp;
  using FixedSolver<FixedEuler<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

  FixedEuler(T_Func& _f)
    :
      FixedSolver<FixedEuler<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                                ( _f.size() ),
      Euler<T_Func, vecT, vecT_traits>		( _f )
    {
      reset();
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for AdaptiveSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1) const
    {
      euler_step(x, y, yp, h, y1);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // class FixedEuler


// Second-order Midpoint
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedMidpoint
  : public FixedSolver<FixedMidpoint<T_Func,vecT,vecT_traits>,
                       vecT, vecT_traits>,
    public Midpoint<T_Func, vecT, vecT_traits>
{
protected:
  using FixedSolver<FixedMidpoint<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::x;
  using FixedSolver<FixedMidpoint<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::dx;
  using FixedSolver<FixedMidpoint<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::y;
  using FixedSolver<FixedMidpoint<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::yp;
  using FixedSolver<FixedMidpoint<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

  FixedMidpoint(T_Func& _f)
    :
      FixedSolver<FixedMidpoint<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                                ( _f.size() ),
      Midpoint<T_Func, vecT, vecT_traits>	( _f )
    {
      reset();
    }

  void OneStep(const stepT h, vecT& y1)
    {
      midpoint_step(x, y, yp, h, y1);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // class FixedMidpoint


template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedRK4
  : public FixedSolver<FixedRK4<T_Func,vecT,vecT_traits>, vecT, vecT_traits>,
    public RK4<T_Func, vecT, vecT_traits>
{
protected:
  using FixedSolver<FixedRK4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::x;
  using FixedSolver<FixedRK4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::dx;
  using FixedSolver<FixedRK4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::y;
  using FixedSolver<FixedRK4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::yp;
  using FixedSolver<FixedRK4<T_Func,vecT,vecT_traits>,
		    vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::step_type	stepT;

  FixedRK4(T_Func& _f)
    :
      FixedSolver<FixedRK4<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                        ( _f.size() ),
      RK4<T_Func, vecT, vecT_traits>	( _f )
    {
      reset();
    }

  void OneStep(const stepT h, vecT& y1)
    {
      rk4_step(x, y, yp, h, y1);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // class FixedRK4


//
// Adaptive Step Size Implementations
//

//
// First-order Euler
//
//  Really second-order because doing two half-steps buys an extra
//  order.  However, we only know the error on the 1st order result,
//  so we tell AdaptiveSolver we are first order (last argument of
//  AdaptiveSolver constructor is 1).
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveEuler
  : public AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>, vecT,
                          vecT_traits>,
    public Euler<T_Func, vecT, vecT_traits>
{
protected:
  using AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::x;
  using AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::dx;
  using AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::y;
  using AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yp;
  using AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yerr;
  using AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator	It;
  typedef typename vecT_traits::const_iterator	CIt;
#endif

private:
  // Working variables for OneStep.
  vecT y_mid;		// y at midpoint of interval.
  vecT yp_mid;		// Derivative of y at midpoint of interval.
  vecT yh;		// y after one large step.

  static const int order = 1;

public:
  AdaptiveEuler(T_Func& _f)
    :
      AdaptiveSolver<AdaptiveEuler<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                                ( _f.size(), order ),
      Euler<T_Func, vecT, vecT_traits>		( _f ), 
      y_mid					( n ),
      yp_mid					( n ),
      yh					( n )
    {
      reset();
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for AdaptiveSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      stepT h_mid = 0.5L*h, x_mid = x + h_mid;
      const magT corr = 1;		// Correction is 1/(2^order - 1).

      // Two small steps.
      euler_step(x, y, yp, h_mid, y_mid);
      func(x_mid, y_mid, yp_mid);
      euler_step(x_mid, y_mid, yp_mid, h_mid, y1);

      // One large step.
      euler_step(x, y, yp, h, yh);

      // Compute error and a 2nd order correction.
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; i++) {
	yerr[i] = y1[i] - yh[i];
	y1[i] += corr * yerr[i];
      }
#else
      // The only reason not to use two simple assignments here is to
      // avoid two loops.
      CIt yhit = yh.begin();
      for (It yerrit = yerr.begin(), y1it = y1.begin();
	   yerrit != yerr.end(); ++yerrit, ++y1it, ++yhit)
	{
	  *yerrit = *y1it - *yhit;
	  *y1it += corr * (*yerrit);
	}
#endif
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // class AdaptiveEuler


//
// Second-order Midpoint Method
//
//  Really third-order because doing two half-steps buys an extra
//  order.  However, we only know the error on the 2nd order result,
//  so we tell AdaptiveSolver we are second order (last argument of
//  AdaptiveSolver constructor is 2).
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveMidpoint
  : public AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>, vecT,
                          vecT_traits>,
    public Midpoint<T_Func, vecT, vecT_traits>
{
protected:
  using AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::x;
  using AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::dx;
  using AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::y;
  using AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yp;
  using AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yerr;
  using AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator	It;
  typedef typename vecT_traits::const_iterator	CIt;
#endif

private:
  // Working variables for OneStep.
  vecT y_mid;		// y at midpoint of interval.
  vecT yp_mid;		// Derivative of y at midpoint of interval.
  vecT yh;		// y after one large step.

  static const int order = 2;

public:
  AdaptiveMidpoint(T_Func& _f)
    :
      AdaptiveSolver<AdaptiveMidpoint<T_Func,vecT,vecT_traits>,
                     vecT, vecT_traits>
                                                ( _f.size(), order ),
      Midpoint<T_Func, vecT, vecT_traits>	( _f ),
      y_mid					( n ),
      yp_mid					( n ),
      yh					( n )
    {
      reset();
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for AdaptiveSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      stepT h_mid = 0.5L*h, x_mid = x + h_mid;
      const magT corr = 1./3.L;		// Correction is 1/(2^order - 1).

      // Two small steps.
      midpoint_step(x, y, yp, h_mid, y_mid);
      func(x_mid, y_mid, yp_mid);
      midpoint_step(x_mid, y_mid, yp_mid, h_mid, y1);

      // One large step.
      midpoint_step(x, y, yp, h, yh);

      // Compute error and a 3rd order correction.
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; i++) {
	yerr[i] = y1[i] - yh[i];
	y1[i] += corr * yerr[i];
      }
#else
      {
	CIt yhit = yh.begin();
	for (It yerrit = yerr.begin(), y1it = y1.begin();
	     yerrit != yerr.end(); ++yerrit, ++y1it, ++yhit)
	{
	  *yerrit = *y1it - *yhit;
	  *y1it += corr * (*yerrit);
	}
      }
#endif
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // class AdaptiveMidpoint


//
// Fourth-order Runge-Kutta
//
//  Really fifth-order because doing two half-steps buys an extra
//  order.  However, we only know the error on the 4th order result,
//  so we tell AdaptiveSolver we are fourth order (last argument of
//  AdaptiveSolver constructor is 4).
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveRK4
  : public AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>, vecT,
                          vecT_traits>,
    public RK4<T_Func, vecT, vecT_traits>
{
protected:
  using AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::x;
  using AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::dx;
  using AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::y;
  using AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yp;
  using AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yerr;
  using AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator	It;
  typedef typename vecT_traits::const_iterator	CIt;
#endif

private:
  // Working variables for OneStep.
  vecT y_mid;		// y at midpoint of interval.
  vecT yp_mid;		// Derivative of y at midpoint of interval.
  vecT yh;		// y after one large step.

  static const int order = 4;

public:
  AdaptiveRK4(T_Func& _f)
    :
      AdaptiveSolver<AdaptiveRK4<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                        ( _f.size(), order ),
      RK4<T_Func, vecT, vecT_traits>	( _f ),
      y_mid				( n ),
      yp_mid				( n ),
      yh				( n )
    {
      reset();
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for AdaptiveSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      stepT h_mid = 0.5L*h, x_mid = x + h_mid;
      const magT corr = 1./15.L;	// Correction is 1/(2^order - 1).

      // Two small steps.
      rk4_step(x, y, yp, h_mid, y_mid);
      func(x_mid, y_mid, yp_mid);
      rk4_step(x_mid, y_mid, yp_mid, h_mid, y1);

      // One large step.
      rk4_step(x, y, yp, h, yh);

      // Compute error and a 5th order correction.
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; i++) {
	yerr[i] = y1[i] - yh[i];
	y1[i] += corr * yerr[i];
      }
#else
      {
	CIt yhit = yh.begin();
	for (It yerrit = yerr.begin(), y1it = y1.begin();
	     yerrit != yerr.end(); ++yerrit, ++y1it, ++yhit)
	{
	  *yerrit = *y1it - *yhit;
	  *y1it += corr * (*yerrit);
	}
      }
#endif
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

}; // class AdaptiveRK4


//
// Fifth-order Runge-Kutta-Fehlberg Cash-Karp Method
//
//  This implementation is not derived from a separate Method, because
//  it provides its own error estimate (embedded method).  Therefore
//  it interfaces directly to AdaptiveSolver.  We only know the error
//  on the 4th order result, so we tell AdaptiveSolver we are fourth
//  order (last argument of AdaptiveSolver constructor is 4).
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveRKCashKarp
  : public AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>, vecT,
                          vecT_traits>
{
protected:
  using AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::x;
  using AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::dx;
  using AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::y;
  using AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yp;
  using AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yerr;
  using AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
  typedef typename vecT_traits::vec_mag_type	vecmagT;

private:
  // Working variables for OneStep
  vecT ak1, ak2, ak3, ak4, ak5, ytemp;

  // Parameters for Cash-Karp method, given in NRC 2nd p. 717 (here 0 indexed).
  /** These should really be static const **/
  const stepT
    a1,  a2,  a3,  a4,  a5,
    b10,
    b20, b21,
    b30, b31, b32,
    b40, b41, b42, b43,
    b50, b51, b52, b53, b54,
     c0,       c2,  c3,       c5,
    dc0,      dc2, dc3, dc4, dc5;

  // These values obtained by Cash and Karp (1990) require one less
  // function evaluation than those of Dorman and Prince (1980, see
  // Stoer and Bulirsch p. 454), but the Dorman and Price method gives
  // yp (the derivative at y[x+h]) for the next step.  Similar cost,
  // but the Dorman and Price method would require telling
  // AdaptiveSolver not to re-evaluate yp.

  static const int order = 4;

public:
  // Constructor: absolute and relative errors are equal scalars.
  AdaptiveRKCashKarp(T_Func& _f)
    :
      AdaptiveSolver<AdaptiveRKCashKarp<T_Func,vecT,vecT_traits>, vecT,
                     vecT_traits>
                                ( _f.size(), order ),
      ak1			( n ),
      ak2			( n ),
      ak3			( n ),
      ak4			( n ),
      ak5			( n ),
      ytemp(n),
      a1  (0.2L),		a2  (0.3L),
      a3  (0.6L),		a4  (1.L),
      a5  (0.875L),
      b10 (0.2L),
      b20 (3./40.L),		b21 (9./40.L),
      b30 (0.3L),		b31 (-0.9L),		b32 (1.2L),
      b40 (-11./54.L),		b41 (2.5L),		b42 (-70./27.L),
      b43 (35./27.L),		
      b50 (1631./55296.L),  	b51 (175./512.L),	 b52 (575./13824.L),
      b53 (44275./110592.L),	b54 (253./4096.L),
      c0  (37./378.L),
      c2  (250./621.L),
      c3  (125./594.L),	        c5 (512./1771.L),
      dc0 (c0-2825./27648.L),
      dc2 (c2-18575./48384.L),	dc3 (c3-13525./55296.L),
      dc4 (-277./14336.L),      dc5 (c5-0.25L),
      func(_f)
    {
      reset();
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for AdaptiveSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  void OneStep(const stepT h, vecT& y1)
    {
      for (int i = 0; i < n ; ++i)
        ytemp[i] = y[i] + b10*h*yp[i];

      func(x + a1*h, ytemp, ak1);

      for (int i = 0; i < n ; ++i)
        ytemp[i] = y[i] + h*(b20*yp[i] + b21*ak1[i]);

      func(x + a2*h, ytemp, ak2);

      for (int i = 0; i < n ; ++i)
        ytemp[i] = y[i] + h*(b30*yp[i] + b31*ak1[i] + b32*ak2[i]);

      func(x + a3*h, ytemp, ak3);
      
      for (int i = 0; i < n ; ++i)
        ytemp[i] = y[i] + h*(b40*yp[i] + b41*ak1[i] + b42*ak2[i] + b43*ak3[i]);

      func(x + a4*h, ytemp, ak4);

      for (int i = 0; i < n ; ++i)
        ytemp[i] = y[i] + h*(b50*yp[i] + b51*ak1[i] + b52*ak2[i]
			     + b53*ak3[i] + b54*ak4[i]);

      func(x + a5*h, ytemp, ak5);

      for (int i = 0; i < n ; ++i)
        y1[i] = y[i] + h*(c0*yp[i] + c2*ak2[i] + c3*ak3[i] + c5*ak5[i]);

      for (int i = 0; i < n ; ++i)
        yerr[i] = h*(dc0*yp[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i]
		     + dc5*ak5[i]);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

  T_Func& func;

}; // class AdaptiveRKCashKarp

//
// Fourth-order Generalized Runge-Kutta Methods
//
//  The implementations presented here are due to Kaps and Rentrop.
//  The original method is due to Rosenbrock and Wanner.
//
//  This implementation is not derived from a separate Method, because
//  it provides its own error estimate (embedded method).  Therefore
//  it interfaces directly to AdaptiveSolver.  We only know the error
//  on the 3rd order result, so we tell AdaptiveSolver we are third
//  order (last argument of AdaptiveSolver constructor is 3).
//
//  Even though it is not implicit, the GRK methods have very good
//  stability properties and are ideal for stiff systems.
//
//  References:
//
//  Stoer&Bulirsch 2nd ed, pp. 491-3.
//
//  P. Kaps and P. Rentrop, "Generalized Runge-Kutta Methods of Order
//  Four with Stepsize control for Stiff Ordinary Differential
//  Equations," Numerische Mathematik 33, 55-68 (1979).
//
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveGRK
  : public AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>, vecT,
                          vecT_traits>
{
protected:
  using AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::x;
  using AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::dx;
  using AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::y;
  using AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yp;
  using AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::yerr;
  using AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>::n;

public:
  typedef typename vecT_traits::value_type	T;
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
  typedef typename vecT_traits::matrix_type	matrixT;

private:
  matrixT	Jac, A;			// Jacobian of equations
  vecT		fs1, fs2, fs3, fs4;	// f^* (S&B p. 492)
  vecT		yy1, yy2;
  vecT		yy1p, yy2p;
  int		*idx;			// Row permuations for LUdecomp

  // Parameters for Kaps-Rentrop GRK method.
  const stepT
    g,
    g21,
    g31, g32,
    g41, g42, g43,
    b21, b31, b32,
    c1, c2, c3, c4,
    ch1, ch2, ch3;

  static const int order = 3;

public:
  AdaptiveGRK(T_Func& _f)
    :
      AdaptiveSolver<AdaptiveGRK<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                ( _f.size(), order ),
      Jac			( n,n ),
      A				( n,n ),
      fs1			( n ),
      fs2			( n ),
      fs3			( n ),
      fs4			( n ),
      yy1			( n ),
      yy2			( n ),
      yy1p			( n ),
      yy2p			( n ),
#     if defined(RODENT_GRK4A)

      g   (0.395),
      g21 (-0.767672395484/g),
      g31 (-0.851675323742/g), g32 (0.522967289188/g),
      g41 (0.288463109545/g),  g42 (0.0880214273381/g), g43 (-0.337389840627/g),
      b21 (0.438),
      b31 (0.796920457938),    b32 (0.0730795420615),
       c1 (0.199293275701),     c2 (0.482645235674),     c3 (0.0680614886256),
       c4 (0.25),
      ch1 (0.346325833758),    ch2 (0.285693175712),    ch3 (0.367980990530),

#     elif defined(RODENT_GRKSB)
      // S&B 2nd p. 492: These values have proven pretty lousy: the
      // method seems to be first order and doesn't work well on stiff
      // problems. Maybe the parameters are wrong in S&B.

      g   (0.220428410),
      g21 (0.822867461/g),
      g31 (0.695700194/g), g32 (0./g),
      g41 (3.90481342/g),  g42 (0./g),        g43 (1./g),
      b21 (-0.554591416),
      b31 (0.252787696),   b32 (1.),
       c1 (0.545211088),    c2 (0.301486480),  c3 (0.177064668),
       c4 (-0.0237622363),
      ch1 (-0.162871035),  ch2 (1.18215360),  ch3 (-0.0192825995),

#     else

      g   (0.231),
      g21 (-0.270629667752/g),
      g31 (0.311254483294/g), g32 (0.00852445628482/g),
      g41 (0.282816832044/g), g42 (-0.457959483281/g), g43 (-0.111208333333/g),
      b21 (0.462),
      b31 (-0.0815668168327), b32 (0.961775150166),
       c1 (0.217487371653),    c2 (0.486229037990),     c3 (0.),
       c4 (0.296283590357),
      ch1 (-0.717088504499),  ch2 (1.77617912176),     ch3 (-0.0590906172617),

#     endif
      func(_f)
    {
      idx = new int[n];

      reset();
    }

  // The constructor passes the number of variables according (_f.size())
  // to the constructor for AdaptiveSolver, which passes it to
  // SolverBase.  SolverBase has member n, which is then inherited.

  ~AdaptiveGRK() { delete[] idx; }

  void OneStep(const stepT h, vecT& y1)
    {
      int	perm;

      // Find fs1
      // Jacobian of right-hand side function (scaled by h):
      func.Jacobian(x, y, yp, h, Jac);
      for (int i = 0, j = 0; i < n; ++i, ++j) {
	A(i,j) = -g*Jac(i,j);
      }
      for (int i = 0; i < n; ++i) {
	A(i,i) += 1.;
	fs1[i] = yp[i];
      }

      // LU-decompose the matrix.
      jlt::LUdecomp<T,matrixT>(A, idx, &perm);

      // Solve the equation A.fs1 = yp.
      jlt::LUbacksub<T,matrixT>(A, idx, &fs1[0]);

      for (int i = 0; i < n; ++i) yy1[i] = y[i] + h*b21*fs1[i];

      func(x, yy1, yy1p);

      for (int i = 0; i < n; ++i) fs2[i] = yy1p[i] + g21*fs1[i];

      // Solve the equation A.(fs2 + g21 fs1) = b.
      jlt::LUbacksub<T,matrixT>(A, idx, &fs2[0]);

      for (int i = 0; i < n; ++i) {
	fs2[i] -= g21*fs1[i];
	yy2[i] = y[i] + h*(b31*fs1[i] + b32*fs2[i]);
      }

      func(x, yy2, yy2p);

      for (int i = 0; i < n; ++i) fs3[i] = yy2p[i] + g31*fs1[i] + g32*fs2[i];

      // Solve the equation A.(fs3 + g31 fs1 + g32 fs2) = b.
      jlt::LUbacksub<T,matrixT>(A, idx, &fs3[0]);

      for (int i = 0; i < n; ++i) fs3[i] -= (g31*fs1[i] + g32*fs2[i]);

      for (int i = 0; i < n; ++i) {
	fs4[i] = yy2p[i] + g41*fs1[i] + g42*fs2[i] + g43*fs3[i];
      }

      // Solve the equation A.(fs4 + g41 fs1 + g42 fs2 + g43 fs3) = b.
      jlt::LUbacksub<T,matrixT>(A, idx, &fs4[0]);

      for (int i = 0; i < n; ++i) {
	fs4[i] -= (g41*fs1[i] + g42*fs2[i] + g43*fs3[i]);
      }

      for (int i = 0; i < n; ++i) {
	y1[i] = y[i] + h*(c1*fs1[i] + c2*fs2[i] + c3*fs3[i] + c4*fs4[i]);
	yy1[i] = y[i] + h*(ch1*fs1[i] + ch2*fs2[i] + ch3*fs3[i]);
	yerr[i] = y1[i] - yy1[i];
      }
      return;
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      func(x, y, yp);
    }

  T_Func& func;

}; // AdaptiveGRK

} // namespace rodent

#endif // RODENT_EXPLICITRK_HPP
