#ifndef RODENT_BLITZ_TRAITS_HPP
#define RODENT_BLITZ_TRAITS_HPP

#include <jlt/matrix.hpp>
#include <jlt/math.hpp>
#include <blitz/vector.h>
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
#ifdef RODENT_ITERATOR_LOOPS
  typedef typename blitz::Vector<value_type>::T_iterator	iterator;
  typedef typename blitz::Vector<value_type>::T_iterator	const_iterator;
  typedef typename vec_mag_type::T_iterator		mag_iterator;
  typedef typename vec_mag_type::T_iterator		const_mag_iterator;
#endif
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
};

} // namespace rodent

#endif // RODENT_BLITZ_TRAITS_HPP
