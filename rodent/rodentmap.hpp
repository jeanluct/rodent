//
// Copyright (c) 2004-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef RODENT_RODENTMAP_HPP
#define RODENT_RODENTMAP_HPP

#ifndef __PGI
#  include <cassert>
#else
#  include <assert.h>
#endif
#include <rodent/traits.hpp>
#include <rodent/base.hpp>


//
// Provide functionality for maps that mimics rodent for ODE's.
//

/* Needs updating: many changes to rodent interface are not reflected here. */

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
  void setState(const stepT it0, const vecT& y0)
    {
      iter = (int)it0;
      y = y0;
    }

  // Iterate the map to x1.
  stepT integrateTo(const stepT x1, vecT& y1)
    {
      if (x1 == (stepT)iter) return x1;	// Already there, do nothing.

      while (takeStep(y1) < x1);

      return x1;
    }

  stepT integrateTo(const stepT x1)
    {
      vecT y1(func.size());

      return integrateTo(x1, y1);
    }

  // Iterate the map once.
  stepT takeStep(vecT& y1)
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
