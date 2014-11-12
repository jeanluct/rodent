//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef RODENT_BLITZ_TRAITS_HPP
#define RODENT_BLITZ_TRAITS_HPP

#include <jlt/matrix.hpp>
#include <jlt/math.hpp>
#include <blitz/vector2.h>
#include <rodent/traits.hpp>


namespace rodent {

template<class float_type>
struct vec_traits<blitz::Vector<float_type> >
{
  typedef blitz::Vector<float_type>				vec_type;
  typedef typename blitz::Vector<float_type>::T_numtype 	value_type;
  typedef value_type						step_type;
  typedef value_type						mag_type;
  typedef blitz::Vector<mag_type>				vec_mag_type;
  typedef jlt::matrix<value_type>				matrix_type;
  static inline mag_type	absval(step_type _x) { return jlt::Abs(_x); }
  static inline mag_type	mag(value_type _x) { return jlt::Abs(_x); }

  // For Blitz, the copy constructor has reference semantics, which
  // means that the new vector shares the same memory as the original.
  // The assignment operator (=) copies memory, yielding two disjoint
  // objects.  Hence, use copy() method when duplicating data with
  // copy constructor.
  static inline blitz::Vector<value_type> copy(const blitz::Vector<value_type>& _v)
    {
      return _v.copy();
    }

  /* Probably don't need to do this.  Use iterators everywhere instead. */
#if 0
  static inline value_type at(const blitz::Vector<value_type>& _v, const int i)
    {
      return _v(i);
    }

  static inline value_type& at(blitz::Vector<value_type>& _v, const int i)
    {
      return _v(i);
    }
#endif
};

} // namespace rodent

#endif // RODENT_BLITZ_TRAITS_HPP
