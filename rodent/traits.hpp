//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef RODENT_TRAITS_HPP
#define RODENT_TRAITS_HPP

#include <vector>
#include <complex>
#include <list>
#include <cmath>
#include <jlt/matrix.hpp>


namespace rodent {

  //
  // The default template declaration works for all standard library
  // containers of simple types, such as vector<double> and list<float>.
  //
  template<class vecT>
  struct vec_traits
  {
    //
    // vec_type		The trait (vector) type.
    // value_type	The type of the dependent variables (T).
    // step_type	The type of the independent variable ("time step").
    // mag_type		The positive magnitude of vecT.
    // vec_mag_type	A vector container of mag_type.
    // matrix_type	A matrix container (for the Jacobian).
    // mag		A function that returns the positive magnitude of x.
    //
    typedef vecT				vec_type;
    using value_type = typename vecT::value_type;
    using step_type = value_type;
    using mag_type = value_type;
    using vec_mag_type = std::vector<mag_type>;
    using matrix_type = jlt::matrix<value_type>;

    static inline mag_type	absval(step_type _x) { return std::abs(_x); }
    static inline mag_type	mag(value_type _x) { return std::abs(_x); }
    static inline vec_type	copy(const vec_type& _v) { return _v; }
  };

  //
  // Need to specialize traits class when the independent variable
  // is not of value_type, such as for complex types.
  //
  template<class float_type>
  struct vec_traits<std::vector<std::complex<float_type> > >
  {
    using vec_type = std::vector<std::complex<float_type> >;
    using value_type = std::complex<float_type>;
    using step_type = float_type;
    using mag_type = float_type;
    using vec_mag_type = std::vector<mag_type>;
    using matrix_type = jlt::matrix<value_type>;

    static inline mag_type 	absval(step_type _x) { return std::abs(_x); }
    static inline mag_type	mag(value_type _x) { return std::abs(_x); }
    static inline vec_type	copy(const vec_type& _v) { return _v; }
  };

} // namespace rodent

#endif // RODENT_TRAITS_HPP
