#ifndef RODENT_BASE_HPP
#define RODENT_BASE_HPP

//#ifdef RODENT_DEBUG
//#  define MAT_CHECK_BOUNDS
//#  define VEC_CHECK_BOUNDS
//#endif

// Default vector type is std vector class.
#include <vector>
typedef std::vector<double> rodent_vec;

#include <iostream>
#include <jlt/exceptions.hpp>
#include <jlt/stlio.hpp>


namespace rodent {

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
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename vecT_traits::iterator		It;
  typedef typename vecT_traits::const_iterator		CIt;
  typedef typename vecT_traits::mag_iterator		Itmag;
  typedef typename vecT_traits::const_mag_iterator	CItmag;
#endif

protected:
  const int n;				// Number of integration variables.

  stepT	x;				// Current independent coordinate.
  stepT	dx;				// Current step size.
  const stepT dx_init;			// Initial step size.
  const magT dx_min;			// Minimum step size.
  vecT 	y;				// State vector at x.
  vecT 	yp;				// Derivative vector at x.

#ifdef RODENT_DEBUG
  int n_good_steps, n_bad_steps;	// Number of good and bad steps.
#endif

private:
  T_Control& Control() { return static_cast<T_Control&>(*this); }

  const unsigned long int max_steps;	// Maximum number of internal steps.
                                        // Largest is 2^32 =~ 4,000,000,000.

public:
  // Constructor: Initialize internal state and save first data point.
  SolverBase(const int _n, const stepT x0, const vecT& y0,
	     const stepT dx0, const magT dx_min0 = 0)
    : n(_n), x(x0), dx(dx0), dx_init(dx0), dx_min(dx_min0),
      y(vecT_traits::copy(y0)), yp(_n), max_steps(10000000)
    {
    }

  // Re-initialize the integrator, recalculating the derivative yp.
  void Restart(const stepT x0, const vecT& y0, const stepT dx0)
    {
      // Check if step size is already too small.
      if (vecT_traits::absval(dx0) < dx_min) {
#ifdef __EXCEPTIONS
	_THROW(jlt::stepsize_too_small<magT>
	       ("New stepsize too small in rodent::Restart.",dx0));
#else
	std::cerr << "New stepsize too small in rodent::Restart.\n";
	std::exit(1);
#endif
      }

      x = x0;
      dx = dx0;
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; ++i) y[i] = y0[i];
#else
      y = y0;
#endif
      // Update derivative vector yp at x.
      Control().reset();
    }

  // Re-initialize the integrator, recalculating the derivative yp.
  // Set step size to same as initially.
  void Restart(const stepT x0, const vecT& y0)
    {
      Restart(x0,y0,dx_init);
    }

  // Re-initialize the integrator, explicitly providing the derivative yp.
  void Restart(const stepT x0, const vecT& y0, const stepT dx0,
	       const vecT& yp0)
    {
      if (vecT_traits::absval(dx0) < dx_min) {
#ifdef __EXCEPTIONS
	_THROW(jlt::stepsize_too_small<magT>
	       ("New stepsize too small in rodent::Restart.",dx0));
#else
	std::cerr << "New stepsize too small in rodent::Restart.\n";
	std::exit(1);
#endif
      }

      x = x0;
      dx = dx0;
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; ++i) y[i] = y0[i];
#else
      y = y0;
#endif
      // Derivative vector yp at x.
      yp = yp0;
    }

  // Re-initialize the integrator, explicitly providing the derivative yp.
  // Set step size to same as initially.
  void Restart(const stepT x0, const vecT& y0, const vecT& yp0)
    {
      Restart(x0,y0,dx_init,yp0);
    }

  stepT IntegrateTo(const stepT x1, vecT& y1)
    {
      if (x == x1) {
	// Already there, do nothing, but copy current state to y1.
#ifndef RODENT_ITERATOR_LOOPS
	for (int i = 0; i < n; ++i) y1[i] = y[i];
#else
	y1 = y;
#endif
	return x;
      }

      // Make sure we are going in the right direction: reverse the
      // sign of dx if needed, but keep the same magnitude.
      if (dx*(x1 - x) < 0) dx = -dx;

      stepT dxold = dx;

#ifdef RODENT_DEBUG
      n_good_steps = n_bad_steps = 0;
#endif

      for (unsigned long int nstp = 1; nstp <= max_steps; ++nstp) {
	bool lastStep = false;

	// Check to make sure that we don't overshoot x1.
	if ((dx > 0 && x + dx >= x1) || (dx < 0 && x + dx <= x1)) {
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
#ifndef RODENT_ITERATOR_LOOPS
	for (int i = 0; i < n; ++i) y[i] = y1[i];
#else
	y = y1;
#endif

	if (lastStep) {
#ifdef RODENT_DEBUG
	  std::cerr << "rodent::IntegrateTo   Steps--  good = "
		    << n_good_steps;
	  std::cerr << "   bad = " << n_bad_steps;
	  std::cerr << "   x = " << x << endl;
	  std::cerr << "   y = " << y << endl;
	  std::cerr << "  yp = " << yp << endl;
	  std::cerr << "   dx = " << dx << endl;
#endif
	  dx = dxold;	// There might be a bit of an error here...
	  return x;
	}
      }
      _THROW(jlt::too_many_steps
	      ("Too many steps taken in IntegrateTo ", max_steps));

      return x;
    }

  stepT IntegrateTo(const stepT x1)
    {
      vecT y1(n);

      return IntegrateTo(x1, y1);
    }

  stepT IntegrateTo(const stepT x1, vecT& y1, vecT& y1p)
    {
      IntegrateTo(x1, y1);
 
#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; ++i) {
	y1p[i] = yp[i];
      }
#else
      y1p = yp;
#endif

      return x;
    }

  // Integrates to x+dx1.
  stepT IntegrateAddTo(const stepT dx1)
    {
      return IntegrateTo(x + dx1);
    }

  // Take one step, return the new x.
  stepT TakeStep(vecT& y1, vecT& y1p)
    {

#ifdef RODENT_DEBUG
      n_good_steps = n_bad_steps = 0;
#endif

      Control().Step(y1);

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; ++i) {
	y[i] = y1[i];
	y1p[i] = yp[i];
      }
#else
      y = y1;
      y1p = yp;
#endif

#ifdef RODENT_DEBUG
      std::cerr << "rodent::TakeStep   Steps--  good = " << n_good_steps;
      std::cerr << "   bad = " << n_bad_steps << endl;
#endif

      return x;
    }

  // Take one step, return the new x.
  stepT TakeStep(vecT& y1)
    {

#ifdef RODENT_DEBUG
      n_good_steps = n_bad_steps = 0;
#endif

      Control().Step(y1);

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; ++i) {
	y[i] = y1[i];
      }
#else
      y = y1;
#endif

#ifdef RODENT_DEBUG
      std::cerr << "rodent::TakeStep   Steps--  good = " << n_good_steps;
      std::cerr << "   bad = " << n_bad_steps << endl;
#endif

      return x;
    }

  // Take one step, return the new x.
  stepT operator++()
    {
      vecT y1(n);

#ifdef RODENT_DEBUG
      n_good_steps = n_bad_steps = 0;
#endif

      Control().Step(y1);

#ifndef RODENT_ITERATOR_LOOPS
      for (int i = 0; i < n; ++i) y[i] = y1[i];
#else
      y = y1;
#endif

#ifdef RODENT_DEBUG
      std::cerr << "rodent::++   Steps--  good = " << n_good_steps;
      std::cerr << "   bad = " << n_bad_steps << endl;
#endif

      return x;
    }

  const vecT& operator()(const stepT x1)
    {
      IntegrateTo(x1);

      return y;
    }

  stepT operator()(const stepT x1, vecT& y1)
    {
      return IntegrateTo(x1,y1);
    }

  stepT operator()(const stepT x1, vecT& y1, vecT& yp1)
    {
      return IntegrateTo(x1,y1,yp1);
    }

  stepT stepSize() const
    {
      return dx;
    }

  stepT stepSize(const stepT dx_new)
    {
      dx = dx_new;

      return dx;
    }

  void PrintOn(std::ostream& strm)
    {
      strm << x << "\t" << y << "\t" << yp << std::endl;
    }

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

  const vecT& dependent() const
    {
      return y;
    }

  const vecT& state() const
    {
      return dependent();
    }

  const T& operator[](int i) const
    {
      return y[i];
    }

}; // class SolverBase

} // namespace rodent

#endif // RODENT_BASE_HPP
