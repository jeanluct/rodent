//
// Copyright (c) 2004-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef RODENT_JLT_TRAITS_HPP
#define RODENT_JLT_TRAITS_HPP

#include <jlt/math.hpp>
#include <jlt/mathvector.hpp>
#include <jlt/mathmatrix.hpp>
#include <rodent/traits.hpp>


namespace rodent {

template<class float_type>
struct vec_traits<jlt::mathvector<float_type> >
{
  typedef jlt::mathvector<float_type>				vec_type;
  typedef typename jlt::mathvector<float_type>::value_type 	value_type;
  typedef value_type						step_type;
  typedef value_type						mag_type;
  typedef jlt::mathvector<mag_type>				vec_mag_type;
  typedef jlt::mathmatrix<value_type>				matrix_type;
  static inline mag_type	absval(step_type _x) { return jlt::Abs(_x); }
  static inline mag_type	mag(value_type _x) { return jlt::Abs(_x); }
  static inline vec_type	copy(const vec_type& _v) { return _v; }
};

} // namespace rodent

#endif // RODENT_JLT_TRAITS_HPP
