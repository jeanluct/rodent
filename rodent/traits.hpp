#ifndef RODENT_TRAITS_HPP
#define RODENT_TRAITS_HPP

#include <vector>
#include <complex>
#include <list>

#include <jlt/matrix.hpp>
#include <jlt/math.hpp>


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
    typedef typename vecT::value_type		value_type;
    typedef value_type				step_type;
    typedef value_type				mag_type;
    typedef std::vector<mag_type>		vec_mag_type;
    typedef jlt::matrix<value_type>		matrix_type;
#ifdef RODENT_ITERATOR_LOOPS
    typedef typename vecT::iterator		iterator;
    typedef typename vecT::const_iterator	const_iterator;
    typedef typename vec_mag_type::iterator	mag_iterator;
    typedef typename vec_mag_type::const_iterator	const_mag_iterator;
#endif

    static inline mag_type	absval(step_type _x) { return jlt::Abs(_x); }
    static inline mag_type	mag(value_type _x) { return jlt::Abs(_x); }
    static inline vec_type	copy(const vec_type& _v) { return _v; }
  };

  //
  // Need to specialize traits class when the independent variable
  // is not of value_type, such as for complex types.
  //
  template<class float_type>
  struct vec_traits<std::vector<std::complex<float_type> > >
  {
    typedef std::vector<std::complex<float_type> >	vec_type;
    typedef std::complex<float_type>			value_type;
    typedef float_type					step_type;
    typedef float_type					mag_type;
    typedef std::vector<mag_type>			vec_mag_type;
    typedef jlt::matrix<value_type>			matrix_type;
#ifdef RODENT_ITERATOR_LOOPS
    typedef typename vec_type::iterator		iterator;
    typedef typename vec_type::const_iterator	const_iterator;
    typedef typename vec_mag_type::iterator	mag_iterator;
    typedef typename vec_mag_type::const_iterator	const_mag_iterator;
#endif

    static inline mag_type 	absval(step_type _x) { return jlt::Abs(_x); }
    static inline mag_type	mag(value_type _x) { return abs(_x); }
    static inline vec_type	copy(const vec_type& _v) { return _v; }
  };

} // namespace rodent

#endif // RODENT_TRAITS_HPP
