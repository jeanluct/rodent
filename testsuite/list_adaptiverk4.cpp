#include <iostream>
#include <iomanip>
#include <list>

// Lists don't support random access [].
// Loops must be done with iterators.
// This limits the methods that can be used.
#ifndef RODENT_ITERATOR_LOOPS
#define RODENT_ITERATOR_LOOPS
#endif

#include <rodent/explicitrk.hpp>
#include <rodent/data.hpp>
#include <jlt/math.hpp>
#include <jlt/stlio.hpp>

using namespace rodent;

//
//  An example of using rodent to integrate a list of elements instead
//  of a vector.  The example is quite contrived, but it illustrates
//  the idea.
//

// Shortcuts
typedef std::list<double>::const_iterator	CIt;
typedef std::list<double>::iterator 		It;

//
//  A list of uncoupled exponentials.
//
class ListExp {
private:
  std::list<double> lambda;

  static const int n = 2;

public:
  ListExp(double lambda1_, double lambda2_)
    {
      lambda.push_back(lambda1_);
      lambda.push_back(lambda2_);
    }

  void operator()(double, const std::list<double>& y, std::list<double>& y_dot)
    {
      CIt yit = y.begin(), lambdait = lambda.begin();
      for (It y_dotit = y_dot.begin();
	   y_dotit != y_dot.end(); ++yit, ++y_dotit, ++lambdait)
	{
	  *y_dotit = *lambdait * *yit;
	}
    }

  const std::list<double> Exact(double t, const std::list<double>& yinit) const
    {
      std::list<double> yexact(n);

      CIt yinitit = yinit.begin(), lambdait = lambda.begin();
      for (It yexactit = yexact.begin();
	   yexactit != yexact.end(); ++yexactit, ++yinitit, ++lambdait)
	{
	  *yexactit = *yinitit * jlt::Exp(*lambdait * t);
	}

      return yexact;
    }

  int size() const { return n; }
};


int main()
{
  using std::cout;
  using std::endl;
  using jlt::operator<<;

  ListExp lexp(1.,-1.);
  std::list<double> y;

  double t = 10.0, dtsav = 1.0;

  cout << "Integrate two real exponentials\n";
  cout << "from t = " << 0. << " to " << t;
  cout << " using adaptive RK4.\n";
  cout << "Save data points every " << dtsav << " timestep.\n\n";
  cout << "The elements are stored as a list, rather than a vector.\n\n";

  // Initial conditions
  y.push_back(1.);
  y.push_back(1.);

  AdaptiveRK4<ListExp,std::list<double> >
    lexp_rk(lexp, 0.0, y, 0.01, 0, 1.0e-10);

  DataPoints<AdaptiveRK4<ListExp,std::list<double> >,std::list<double> >
    lexp_data(lexp_rk,0.0,t,dtsav);

  cout.precision(10);

  lexp_data.PrintOn(cout, 0., t);

  cout << "\nExact:\n";
  cout << t << "\t" << lexp.Exact(t, y) << endl;
}
