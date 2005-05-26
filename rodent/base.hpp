#ifndef RODENT_BASE_HPP
#define RODENT_BASE_HPP

#include <iostream>
#include <vector>
#include <jlt/exceptions.hpp>
#include <jlt/stlio.hpp>


namespace rodent {


// Default vector type is std::vector class.
typedef std::vector<double> rodent_vec;


template<class T_Control, class vecT, class vecT_traits>
class SolverBase
{
public:

  // From the traits class vecT_traits:
  //
  //   vecT		The vector container to use for T.
  //   T		The type contained in vecT.
  //   stepT		The independent variable type.
  //   magT		The magnitude type.
  //   matrixT		A matrix container (for the Jacobian).
  //
  typedef typename vecT_traits::value_type	T;
  typedef typename vecT_traits::mag_type	magT;
  typedef typename vecT_traits::step_type	stepT;
  typedef typename vecT_traits::vec_mag_type	vecmagT;
  typedef typename vecT_traits::matrix_type	matrixT;

protected:
  const int n;				// Number of integration variables.

  stepT	x;				// Current independent coordinate.
  stepT	dx;				// Current step size.
  magT	dx_min;				// Minimum step size.
  magT	dx_max;				// Maximum step size.
  vecT 	y;				// State vector at x.
  vecT 	yp;				// Derivative vector at x.

#ifdef RODENT_DEBUG
  int n_good_steps, n_bad_steps;	// Number of good and bad steps.
#endif

private:
  T_Control& Control() { return static_cast<T_Control&>(*this); }

  unsigned long int max_steps;	// Maximum number of internal steps.
                                // Largest is 2^32 =~ 4,000,000,000.

public:

  // Constructor: Initialize internal state and save first data point.
  SolverBase(const int _n)
    :
      n				( _n ),
      x				( 0 ),
      dx			( 0.01 ),
      dx_min			( 0 ),
      dx_max			( 0 ),
      y				( _n ),
      yp			( _n ),
      max_steps			( 10000000 )
    {
    }


  //
  // Methods for Setting Parameters
  //

  T_Control& stepSize(const stepT _dx)
    {
      dx = _dx;

      return Control();
    }

  T_Control& minStepSize(const stepT _dx_min)
    {
      dx_min = vecT_traits::absval(_dx_min);
      return Control();
    }

  T_Control& maxStepSize(const stepT _dx_max)
    {
      dx_max = vecT_traits::absval(_dx_max);
      return Control();
    }

  T_Control& maxSteps(unsigned long int _max_steps)
    {
      max_steps = _max_steps;
      return Control();
    }


  // Re-initialize the integrator, recalculating the derivative yp.
  T_Control& setState(const stepT x0, const vecT& y0)
    {
      x = x0;
      for (int i = 0; i < n; ++i) y[i] = y0[i];

      // Update derivative vector yp at x.
      Control().reset();

      return Control();
    }


  // Re-initialize the integrator, explicitly providing the derivative yp.
  T_Control& setState(const stepT x0, const vecT& y0, const vecT& yp0)
    {
      x = x0;
      for (int i = 0; i < n; ++i) y[i] = y0[i];

      // Derivative vector yp at x.
      yp = yp0;

      return Control();
    }


  //
  // Query Parameters
  //

  stepT stepSize() const
    {
      return dx;
    }


  // These are all aliases to obtain the dependent (state) variable y.

  const vecT& dependent() const
    {
      return y;
    }

  const vecT& getState() const
    {
      return dependent();
    }

  const stepT getState(vecT& _y) const
    {
      _y = y;
      return x;
    }

  const stepT getState(vecT& _y, vecT& _yp) const
    {
      _y = y;
      _yp = yp;
      return x;
    }

  // These are all aliases to obtain the independent variable x.

  const stepT independent() const
    {
      return x;
    }

  const stepT time() const
    {
      return independent();
    }

  const stepT position() const
    {
      return independent();
    }

  // Access individual elements of the state vector.
  const T& operator[](int i) const
    {
      return y[i];
    }


  //
  // Integrate!
  //

  inline stepT integrateTo(const stepT x1, vecT& y1);

  stepT integrateTo(const stepT x1)
    {
      vecT y1(n);

      return integrateTo(x1, y1);
    }

  stepT integrateTo(const stepT x1, vecT& y1, vecT& y1p)
    {
      integrateTo(x1, y1);
 
      for (int i = 0; i < n; ++i)
	{
	  y1p[i] = yp[i];
	}

      return x;
    }

  // Integrates to x+dx1.
  stepT integrateAddTo(const stepT dx1)
    {
      return integrateTo(x + dx1);
    }

  // Take one step, return the new x and y. Provide derivative.
  inline stepT takeStep(vecT& y1, vecT& y1p);

  // Take one step, return the new x and y. Compute derivative.
  inline stepT takeStep(vecT& y1);

  // Take one step, return the new x. Compute derivative.
  inline stepT operator++();

  // Aliases for integrateTo.
  const vecT& operator()(const stepT x1)
    {
      integrateTo(x1);

      return y;
    }

  stepT operator()(const stepT x1, vecT& y1)
    {
      return integrateTo(x1,y1);
    }

  stepT operator()(const stepT x1, vecT& y1, vecT& yp1)
    {
      return integrateTo(x1,y1,yp1);
    }


  //
  // Output
  //

  void PrintOn(std::ostream& strm)
    {
      strm << x << "\t" << y << "\t" << yp << std::endl;
    }

}; // class SolverBase


//
// Member Functions Definitions
//

template<class T_Control, class vecT, class vecT_traits>
typename SolverBase<T_Control,vecT,vecT_traits>::stepT
SolverBase<T_Control,vecT,vecT_traits>::integrateTo(const stepT x1, vecT& y1)
{
  if (x == x1)
    {
      // Already there, do nothing, but copy current state to y1.
      for (int i = 0; i < n; ++i) y1[i] = y[i];
      return x;
    }

  // Make sure we are going in the right direction: reverse the
  // sign of dx if needed, but keep the same magnitude.
  if (dx*(x1 - x) < 0) dx = -dx;

  stepT dxold = dx;

#ifdef RODENT_DEBUG
  n_good_steps = n_bad_steps = 0;
#endif

  for (unsigned long int nstp = 1; nstp <= max_steps; ++nstp)
    {
      bool lastStep = false;

      // Check to make sure that we don't overshoot x1.
      if ((dx > 0 && x + dx >= x1) || (dx < 0 && x + dx <= x1))
	{
	  // Save size of last step, in case we are really close to goal.
	  // This makes it easier to restart the integration.
	  dxold = dx;
	  dx = x1 - x;
	  lastStep = true;
	}

      lastStep &= Control().Step(y1);

      // The logical AND ensures that if the last step fails,
      // then it isn't the last one anymore.

      // Next step starts at y1.
      for (int i = 0; i < n; ++i) y[i] = y1[i];

      if (lastStep)
	{
#ifdef RODENT_DEBUG
	  std::cerr << "rodent::integrateTo   Steps--  good = "
		    << n_good_steps;
	  std::cerr << "   bad = " << n_bad_steps;
	  std::cerr << "   x = " << x << endl;
	  std::cerr << "   y = " << y << endl;
	  std::cerr << "  yp = " << yp << endl;
	  std::cerr << "  dx = " << dx << endl;
#endif
	  dx = dxold;	// There might be a bit of an error here...
	  return x;
	}
    }
  _THROW(jlt::too_many_steps
	 ("Too many steps taken in integrateTo ", max_steps));

  return x;
}


template<class T_Control, class vecT, class vecT_traits>
typename SolverBase<T_Control,vecT,vecT_traits>::stepT
SolverBase<T_Control,vecT,vecT_traits>::takeStep(vecT& y1, vecT& y1p)
{
#ifdef RODENT_DEBUG
  n_good_steps = n_bad_steps = 0;
#endif

  Control().Step(y1);

  for (int i = 0; i < n; ++i)
    {
      y[i] = y1[i];
      y1p[i] = yp[i];
    }

#ifdef RODENT_DEBUG
  std::cerr << "rodent::takeStep   Steps--  good = " << n_good_steps;
  std::cerr << "   bad = " << n_bad_steps << endl;
#endif

  return x;
}


template<class T_Control, class vecT, class vecT_traits>
typename SolverBase<T_Control,vecT,vecT_traits>::stepT
SolverBase<T_Control,vecT,vecT_traits>::takeStep(vecT& y1)
{
#ifdef RODENT_DEBUG
  n_good_steps = n_bad_steps = 0;
#endif

  Control().Step(y1);

  for (int i = 0; i < n; ++i)
    {
      y[i] = y1[i];
    }

#ifdef RODENT_DEBUG
  std::cerr << "rodent::takeStep   Steps--  good = " << n_good_steps;
  std::cerr << "   bad = " << n_bad_steps << endl;
#endif

  return x;
}


template<class T_Control, class vecT, class vecT_traits>
typename SolverBase<T_Control,vecT,vecT_traits>::stepT
SolverBase<T_Control,vecT,vecT_traits>::operator++()
{
  vecT y1(n);

#ifdef RODENT_DEBUG
  n_good_steps = n_bad_steps = 0;
#endif

  Control().Step(y1);

  for (int i = 0; i < n; ++i) y[i] = y1[i];

#ifdef RODENT_DEBUG
  std::cerr << "rodent::++   Steps--  good = " << n_good_steps;
  std::cerr << "   bad = " << n_bad_steps << endl;
#endif

  return x;
}


} // namespace rodent

#endif // RODENT_BASE_HPP
