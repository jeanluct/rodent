//
// Copyright (c) 2004-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef JLT_FIXVECTOR_HPP
#define JLT_FIXVECTOR_HPP

//
// fixvector.hpp
//

// Simple vector of fixed length.
//
// Some problems/notes:
//
// -Memory overhead is a bit too large, since neither finish nor
// end_of_storage pointers are needed.
//
// -Should not be able to resize a FixVec, or equate two FixVecs of
// different sizes.
//
// -fixvectors (and fixmatrices) are potentially much faster than
// vectors, especially because I can use buffers to store temporaries.
//
// -Should fixvector really be derived from vector? There are a bunch
// of issues, such as efficiency, and the fact that I don't want it to
// be possible to equate two fixvectors of different lengths, etc. On
// the other hand, it is a lot of rewriting. Also, do I really want a
// fixvector to be a vector? Can I think of applications where I would
// want to substitute a fixvector for a vector?

#include <jlt/vector.hpp>
#include <jlt/mathvector.hpp>

namespace jlt {

template<class T, unsigned int N>
class fixvector : public vector<T>
{
public:

  //
  // Constructors
  //

  // Vector of size N filled with _x.
  explicit fixvector(typename vector<T>::const_reference _x = T())
    : vector<T>(N,_x) {}

  // Copy constructor.
  fixvector(const vector<T>& _v) : vector<T>(_v) {}
}; // class fixvector

template<class T, unsigned int N, class S = T>
class fixmathvector : public mathvector<T,S>
{
public:

  //
  // Constructors
  //

  // Vector of size N filled with _x.
  explicit fixmathvector(typename vector<T>::const_reference _x = T())
    : mathvector<T,S>(N,_x) {}

  // Copy constructor.
  fixmathvector(const mathvector<T,S>& _v) : mathvector<T,S>(_v) {}

#if 0==1
  template<unsigned int M>
  fixmathvector<T,N,S>& operator=(fixmathvector<T,M,S>& v)
    {
      VECTOR_ASSERT(M == N);
    }
#endif
}; // class fixmathvector

} // namespace jlt

#endif // JLT_FIXVECTOR_HPP
