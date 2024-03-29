//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef RODENT_EXPTIMEDIFF_HPP
#define RODENT_EXPTIMEDIFF_HPP

#include <algorithm>
#include <rodent/solver.hpp>
#include <rodent/traits.hpp>


namespace rodent {

//
// Methods
//

//
// Explicit Time-Differenced Euler Step (ETD1)
//
template<class T_Func, class vecT, class vecT_traits>
class ETD1
{
public:
  typedef typename vecT_traits::step_type	step_type;
  using T = typename vecT_traits::value_type;

private:
  vecT c;			// Diagonal linear part of equations.

public:
  ETD1(T_Func& _f, vecT& _c) : c(vecT_traits::copy(_c)), func(_f) {}

  T_Func& func;

protected:
  void etd1_step(const step_type x0, const vecT& y0, const vecT& yp0,
		 const step_type h, vecT& y1) const
    {
      for (int i = 0; i < func.size(); ++i)
	{
	  // Should ensure h*c is not too small.
	  T expc = Exp(h*c[i]);
	  y1[i] = y0[i] * expc + yp0[i] * (expc - 1)/c[i];
	}
    }

}; // class ETD1


//
// 2nd Order Explicit Time-Differenced Runge-Kutta Step (ETDRK2)
//
template<class T_Func, class vecT, class vecT_traits>
class ETDRK2
{
public:
  using step_type = typename vecT_traits::step_type;
  using T = typename vecT_traits::value_type;

private:
  vecT c, expc;			// Diagonal linear part of equations.

public:
  ETDRK2(T_Func& _f, vecT& _c) : c(vecT_traits::copy(_c)), expc(_f.size()), 
				 func(_f),
				 a(_f.size()), ap(_f.size()) {}

  T_Func& func;

private:
  vecT a, ap;		// Working variables.

protected:
  void etdrk2_step(const step_type x0, const vecT& y0, const vecT& yp0,
		   const step_type h, vecT& y1)
    {
      double xh = x0 + h;

      // Take an etd1 step.
      for (int i = 0; i < func.size(); ++i)
	{
	  // Should ensure h*c is not too small.
	  expc[i] = Exp(h*c[i]);
	  a[i] = y0[i] * expc[i] + yp0[i] * (expc[i] - 1)/c[i];
	}

      // Evaluate derivative at that point.
      this->func(xh, a, ap);

      for (int i = 0; i < func.size(); ++i)
	{
	  y1[i] = a[i] + (ap[i] - yp0[i])*(expc[i] - 1 - h*c[i])/(h*c[i]*c[i]);
	}
    }

}; // class ETDRK2


//
// 3rd Order Explicit Time-Differenced Runge-Kutta Step (ETDRK3)
//
template<class T_Func, class vecT, class vecT_traits>
class ETDRK3
{
public:
  using step_type = typename vecT_traits::step_type;
  using T = typename vecT_traits::value_type;

private:
  vecT c, expc, expc2;		// Diagonal linear part of equations.

public:
  ETDRK3(T_Func& _f, vecT& _c) : c(vecT_traits::copy(_c)),
				 expc(_f.size()), expc2(_f.size()),
				 func(_f),
				 a(_f.size()), ap(_f.size()),
				 b(_f.size()), bp(_f.size()) {}

  T_Func& func;

private:
  vecT a, ap, b, bp;		// Working variables.

protected:
  void etdrk3_step(const step_type x0, const vecT& y0, const vecT& yp0,
		   const step_type h, vecT& y1)
    {
      step_type h_2 = 0.5L*h, x_2 = x0 + h_2, xh = x0 + h;

      for (int i = 0; i < func.size(); ++i)
	{
	  // Should ensure h_2*c is not too small.
	  expc2[i] = Exp(h_2*c[i]);
	  a[i] = y0[i] * expc2[i] + yp0[i] * (expc2[i] - 1)/c[i];
	}

      this->func(x_2, a, ap);

      for (int i = 0; i < func.size(); ++i)
	{
	  expc[i] = expc2[i]*expc2[i];
	  b[i] = y0[i] * expc[i] + (2*ap[i] - yp0[i]) * (expc[i] - 1)/c[i];
	}

      this->func(xh, b, bp);

      for (int i = 0; i < func.size(); ++i)
	{
	  y1[i] = y0[i] * expc[i] +
	    (yp0[i] * (-4 - h*c[i] + expc[i]*(4 - 3*h*c[i] + h*h*c[i]*c[i]))
	     + 4*ap[i] * (2 + h*c[i] + expc[i]*(-2 + h*c[i]))
	     + bp[i] * (-4 - 3*h*c[i] - h*h*c[i]*c[i] + expc[i]*(4 - h*c[i])))
	    /(h*h*c[i]*c[i]*c[i]);
	}
    }

}; // class ETDRK3


//
// 4th Order Explicit Time-Differenced Runge-Kutta Step (ETDRK4)
//
template<class T_Func, class vecT, class vecT_traits>
class ETDRK4
{
public:
  using step_type = typename vecT_traits::step_type;
  using T = typename vecT_traits::value_type;

private:
  vecT c, expc, expc2;		// Diagonal linear part of equations.

public:
  ETDRK4(T_Func& _f, vecT& _c) : c(vecT_traits::copy(_c)),
				 expc(_f.size()), expc2(_f.size()),
				 func(_f),
				 a(_f.size()), ap(_f.size()),
				 b(_f.size()), bp(_f.size()),
				 d(_f.size()), dp(_f.size()) {}

  T_Func& func;

private:
  vecT a, ap, b, bp, d, dp;	// Working variables.

protected:
  void etdrk4_step(const step_type x0, const vecT& y0, const vecT& yp0,
		   const step_type h, vecT& y1)
    {
      step_type h_2 = 0.5L*h, x_2 = x0 + h_2, xh = x0 + h;

      for (int i = 0; i < func.size(); ++i)
	{
	  // Should ensure h_2*c is not too small.
	  expc2[i] = std::exp(h_2*c[i]);
	  a[i] = y0[i] * expc2[i] + yp0[i] * (expc2[i] - 1)/c[i];
	}

      this->func(x_2, a, ap);

      for (int i = 0; i < func.size(); ++i)
	{
	  b[i] = y0[i] * expc2[i] + ap[i] * (expc2[i] - 1)/c[i];
	}

      this->func(x_2, b, bp);

      for (int i = 0; i < func.size(); ++i)
	{
	  d[i] = a[i] * expc2[i] + (2*bp[i] - yp0[i]) * (expc2[i] - 1)/c[i];
	}

      this->func(xh, d, dp);

      for (int i = 0; i < func.size(); ++i)
	{
	  expc[i] = expc2[i]*expc2[i];
	  y1[i] = y0[i] * expc[i] +
	    (yp0[i] * (-4 - h*c[i] + expc[i]*(4 - 3*h*c[i] + h*h*c[i]*c[i]))
	     + 2*(ap[i] + bp[i])*(2 + h*c[i] + expc[i]*(-2 + h*c[i]))
	     + dp[i]*(-4 - 3*h*c[i] - h*h*c[i]*c[i] + expc[i]*(4 - h*c[i])))
	    /(h*h*c[i]*c[i]*c[i]);
	}
    }

}; // class ETDRK4


//
// Fixed Step Size Implementations
//

// First-order Explicit Time Differencing
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedETD1
  : public FixedSolver<FixedETD1<T_Func,vecT,vecT_traits>, vecT, vecT_traits>,
    public ETD1<T_Func, vecT, vecT_traits>
{
  using FSolver = FixedSolver<FixedETD1<T_Func, vecT, vecT_traits>, vecT, vecT_traits>;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  using stepT = typename vecT_traits::step_type;

  FixedETD1(T_Func& _f, vecT& _c)
    :
      FixedSolver<FixedETD1<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                        ( _f.size() ),
      ETD1<T_Func, vecT, vecT_traits>	( _f,_c )
    {
    }

  void OneStep(const stepT h, vecT& y1) const
    {
      this->etd1_step(x, y, yp, h, y1);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      this->func(x, y, yp);
    }

}; // class FixedETD1


// Second-order Runge-Kutta Explicit Time Differencing
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedETDRK2
  : public FixedSolver<FixedETDRK2<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>,
    public ETDRK2<T_Func, vecT, vecT_traits>
{
  using FSolver = FixedSolver<FixedETDRK2<T_Func, vecT, vecT_traits>, vecT, vecT_traits>;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  using stepT = typename vecT_traits::step_type;

  FixedETDRK2(T_Func& _f, vecT& _c)
    :
      FixedSolver<FixedETDRK2<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                                ( _f.size() ),
      ETDRK2<T_Func, vecT, vecT_traits>		( _f,_c )
    {
    }

  void OneStep(const stepT h, vecT& y1)
    {
      this->etdrk2_step(x, y, yp, h, y1);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      this->func(x, y, yp);
    }

}; // class FixedETDRK2


// Third-order Runge-Kutta Explicit Time Differencing
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedETDRK3
  : public FixedSolver<FixedETDRK3<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>,
    public ETDRK3<T_Func, vecT, vecT_traits>
{
  using FSolver = FixedSolver<FixedETDRK3<T_Func, vecT, vecT_traits>, vecT, vecT_traits>;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  using stepT = typename vecT_traits::step_type;

  FixedETDRK3(T_Func& _f, vecT& _c)
    :
      FixedSolver<FixedETDRK3<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                                ( _f.size() ),
      ETDRK3<T_Func, vecT, vecT_traits>		( _f,_c )
    {
    }

  void OneStep(const stepT h, vecT& y1)
    {
      this->etdrk3_step(x, y, yp, h, y1);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      this->func(x, y, yp);
    }

}; // class FixedETDRK3


// Fourth-order Runge-Kutta Explicit Time Differencing
template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class FixedETDRK4
  : public FixedSolver<FixedETDRK4<T_Func,vecT,vecT_traits>,
		       vecT, vecT_traits>,
    public ETDRK4<T_Func, vecT, vecT_traits>
{
  using FSolver = FixedSolver<FixedETDRK4<T_Func, vecT, vecT_traits>, vecT, vecT_traits>;

protected:
  using FSolver::x;
  using FSolver::dx;
  using FSolver::y;
  using FSolver::yp;
  using FSolver::n;

public:
  using stepT = typename vecT_traits::step_type;

  FixedETDRK4(T_Func& _f, vecT& _c)
    :
      FixedSolver<FixedETDRK4<T_Func,vecT,vecT_traits>, vecT, vecT_traits>
                                                ( _f.size() ),
      ETDRK4<T_Func, vecT, vecT_traits>		( _f,_c )
    {
    }

  void OneStep(const stepT h, vecT& y1)
    {
      this->etdrk4_step(x, y, yp, h, y1);
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      this->func(x, y, yp);
    }

}; // class FixedETDRK4


//
// Adaptive Step Size Implementations
//

template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveETDRK3
  : public AdaptiveSolver<AdaptiveETDRK3<T_Func,vecT,vecT_traits>, vecT,
                          vecT_traits>,
    public ETDRK3<T_Func, vecT, vecT_traits>
{
  using ASolver = AdaptiveSolver<AdaptiveETDRK3<T_Func, vecT, vecT_traits>, vecT, vecT_traits>;

protected:
  using ASolver::x;
  using ASolver::dx;
  using ASolver::y;
  using ASolver::yp;
  using ASolver::yerr;
  using ASolver::n;

public:
  using magT = typename vecT_traits::mag_type;
  using stepT = typename vecT_traits::step_type;

private:
  // Working variables for OneStep.
  vecT y_mid;		// y at midpoint of interval.
  vecT yp_mid;		// Derivative of y at midpoint of interval.
  vecT yh;		// y after one large step.

  static const int order = 4;

public:
  AdaptiveETDRK3(T_Func& _f, vecT& _c)
    : AdaptiveSolver<AdaptiveETDRK3<T_Func,vecT,vecT_traits>, 
		     vecT, vecT_traits>
                                                ( _f.size(), order ),
      ETDRK3<T_Func, vecT, vecT_traits>		( _f,_c ), 
      y_mid					( n ),
      yp_mid					( n ),
      yh					( n )
    {
    }

  void OneStep(const stepT h, vecT& y1)
    {
      stepT h_mid = 0.5L*h, x_mid = x + h_mid;
      const magT corr = 1./15.L;	// Correction is 1/(2^order - 1).

      // Two small steps.
      this->etdrk3_step(x, y, yp, h_mid, y_mid);
      this->func(x_mid, y_mid, yp_mid);
      this->etdrk3_step(x_mid, y_mid, yp_mid, h_mid, y1);

      // One large step.
      this->etdrk3_step(x, y, yp, h, yh);

      // Compute error and a 5th order correction.
      for (int i = 0; i < n; i++)
	{
	  yerr[i] = y1[i] - yh[i];
	  y1[i] += corr * yerr[i];
	}
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      this->func(x, y, yp);
    }

}; // class AdaptiveETDRK3


template
<class T_Func, class vecT = rodent_vec, class vecT_traits = vec_traits<vecT> >
class AdaptiveETDRK4
  : public AdaptiveSolver<AdaptiveETDRK4<T_Func,vecT,vecT_traits>, vecT,
                          vecT_traits>,
    public ETDRK4<T_Func, vecT, vecT_traits>
{
  using ASolver = AdaptiveSolver<AdaptiveETDRK4<T_Func, vecT, vecT_traits>, vecT, vecT_traits>;

protected:
  using ASolver::x;
  using ASolver::dx;
  using ASolver::y;
  using ASolver::yp;
  using ASolver::yerr;
  using ASolver::n;

public:
  using magT = typename vecT_traits::mag_type;
  using stepT = typename vecT_traits::step_type;

private:
  // Working variables for OneStep.
  vecT y_mid;		// y at midpoint of interval.
  vecT yp_mid;		// Derivative of y at midpoint of interval.
  vecT yh;		// y after one large step.

  static const int order = 4;
public:
  AdaptiveETDRK4(T_Func& _f, vecT& _c)
    : AdaptiveSolver<AdaptiveETDRK4<T_Func,vecT,vecT_traits>,
		     vecT, vecT_traits>
                                                ( _f.size(), order ),
      ETDRK4<T_Func, vecT, vecT_traits>		( _f,_c ), 
      y_mid					( n ),
      yp_mid					( n ),
      yh					( n )
    {
    }

  void OneStep(const stepT h, vecT& y1)
    {
      stepT h_mid = 0.5L*h, x_mid = x + h_mid;
      const magT corr = 1./15.L;	// Correction is 1/(2^order - 1).

      // Two small steps.
      this->etdrk4_step(x, y, yp, h_mid, y_mid);
      this->func(x_mid, y_mid, yp_mid);
      this->etdrk4_step(x_mid, y_mid, yp_mid, h_mid, y1);

      // One large step.
      this->etdrk4_step(x, y, yp, h, yh);

      // Compute error and a 5th order correction.
      for (int i = 0; i < n; i++)
	{
	  yerr[i] = y1[i] - yh[i];
	  y1[i] += corr * yerr[i];
	}
    }

  void reset()
    {
      // Evaluate derivative vector yp at x.
      this->func(x, y, yp);
    }

}; // class AdaptiveETDRK4


} // namespace rodent

#endif // RODENT_EXPTIMEDIFF_HPP
