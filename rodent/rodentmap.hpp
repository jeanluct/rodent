#ifndef RODENT_RODENTMAP_HPP
#define RODENT_RODENTMAP_HPP

#ifndef __PGI
#  include <cassert>
#else
#  include <assert.h>
#endif
#include <rodent/traits.hpp>
#include <rodent/base.hpp>

using namespace std;

//
// Provide functionality for maps that mimics rodent for ODE's.
//

namespace rodent {

template<class T_Map, class vecT, class vecT_traits = vec_traits<vecT> >
class RodentMap
{
public:
  typedef typename vecT_traits::value_type	T;
  typedef typename vecT_traits::step_type	stepT;

protected:
  int iter;			// The current step.
  vecT y;			// The state vector at iter.

public:
  T_Map& func;

  RodentMap(T_Map& func_, stepT it0, const vecT& y0)
    : iter((int)it0), y(vecT_traits::copy(y0)), func(func_)
    {
      // Verify that the sizes agree.
      assert((int)y.size() == func.size());
    }

  // Reinitialize the map.
  void Restart(const stepT it0, const vecT& y0)
    {
      iter = (int)it0;
      y = y0;
    }

  // Iterate the map to x1.
  stepT IntegrateTo(const stepT x1, vecT& y1)
    {
      if (x1 == (stepT)iter) return x1;	// Already there, do nothing.

      while (TakeStep(y1) < x1);

      return x1;
    }

  stepT IntegrateTo(const stepT x1)
    {
      vecT y1(func.size());

      return IntegrateTo(x1, y1);
    }

  // Iterate the map once.
  stepT TakeStep(vecT& y1)
    {
      // The convention for maps is that they are called with an
      // argument of n=iter and return the value at n+1.  Hence use
      // postincrement.
      func(iter++,y1,y);
      y1 = y;

      return (stepT)iter;
    }
};

} // namespace rodent

#endif // RODENT_RODENTMAP_HPP
