//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#ifndef RODENT_DATA_HPP
#define RODENT_DATA_HPP

#include <iostream>
#include <algorithm>
#include <map>
#include <cassert>
#include <rodent/base.hpp>
#include <rodent/traits.hpp>
#include <jlt/stlio.hpp>


namespace rodent {

template
<
 class T_DataFunc,
 class vecT = rodent_vec,
 class vecT_traits = vec_traits<vecT>
>
class DataPoints
{
public:

  // From the traits class vecT_traits:
  //
  //   vecT		The vector container to use for T.
  //   T		The type contained in vecT.
  //   stepT		The independent variable type.
  //
  typedef typename vecT_traits::value_type	T;
  typedef typename vecT_traits::step_type	stepT;

private:
  // Iterator types
  typedef typename std::map<stepT, vecT>::iterator		mapIter;
  typedef typename std::map<stepT, vecT>::const_iterator	mapCIter;

  std::map<stepT, vecT> y_data;		// Data points.

  stepT dxsav;		// Interval at which points are currently saved.
  stepT x_lastsav;	// Last point saved.

  const unsigned long int max_points;	// Maximum number of saved points.
                                        // Largest is 2^32 =~ 4,000,000,000.

public:
  // static const int n = T_DataFunc::n;	// Number of variables.
  T_DataFunc& data_func;

  // Constructor
  DataPoints(T_DataFunc& f_, const stepT x_min, const stepT x_max,
	     const stepT dxsav_)
    : dxsav(dxsav_), max_points(100000), data_func(f_)
    {
      // With this constructor, the difference between x_min and x_max
      // sets the sign of dxsav.
      if (x_max >= x_min)
	{
	  dxsav = std::abs(dxsav);
	}
      else
	{
	  dxsav = -std::abs(dxsav);
	}

      // Save first point.
      y_data.insert(std::make_pair(x_min,
				   vecT_traits::copy(data_func(x_min))));

      x_lastsav = x_min;

      if (x_max != x_min) Sample(x_max);	// Sample the interval.
    }

  // Constructor
  DataPoints(T_DataFunc& f_, const stepT x_, const stepT dxsav_)
    : dxsav(dxsav_), max_points(100000), data_func(f_)
    {
      // Save first point.
      y_data.insert(make_pair(x_, vecT_traits::copy(data_func(x_))));

      x_lastsav = x_;
    }

  void Sample(int n_sample)
    {
      // Sample the next n_sample data points, at current sampling interval.

      assert(n_sample >= 0);

      stepT xx = x_lastsav;

      for (int i = 1; i <= n_sample; ++i)
	{
	  xx = x_lastsav + i*dxsav;
	  y_data.insert(std::make_pair(xx, vecT_traits::copy(data_func(xx))));
	}
      x_lastsav = xx;
    }

  void Sample(const stepT x2, int n_sample)
    {
      // Sample to x2 with n_sample data points, adjusting dxsav.

      dxsav = (x2 - x_lastsav)/n_sample;

      Sample(n_sample);
    }

  void Sample(const stepT x2, const stepT dxsav_)
    {
      // Sample to x2, using a sampling interval close to dxsav_.

      // Try this interval, but Sample will adjust if necessary.
      dxsav = dxsav_;

      Sample(x2);
    }

  void Sample(const stepT x2)
    {
      // Sample to x2, using a sampling interval close to current dxsav.

      int n_sample = (int)(std::floor(std::abs((x2 - x_lastsav)/dxsav) + 0.5));

      // Adjust the magnitude and sign of dxsav.
      dxsav = (x2 - x_lastsav)/n_sample;

      Sample(n_sample);
    }

  void SampleAndPrintOn(std::ostream& strm, int n_sample)
    {
      // Sample the next n_sample data points, at current sampling interval.
      // Output to strm as we go.

      stepT xx = x_lastsav;

      assert(n_sample >= 0);

      for (int i = 1; i <= n_sample; i++)
	{
	  xx = x_lastsav + i*dxsav;
	  y_data.insert(make_pair(xx, vecT_traits::copy(data_func(xx))));
	  strm << xx << "\t" << y_data[xx] << std::endl;
	}
      x_lastsav = xx;
    }

  void SampleAndPrintOn(std::ostream& strm, const stepT x2, int n_sample)
    {
      // Sample to x2 with n_sample data points, adjusting dxsav.
      // Output to strm as we go.

      dxsav = (x2 - x_lastsav)/n_sample;

      SampleAndPrintOn(strm,n_sample);
    }

  void SampleAndPrintOn(std::ostream& strm, const stepT x2, const stepT dxsav_)
    {
      // Sample to x2, using a sampling interval close to dxsav_.
      // Output to strm as we go.

      dxsav = dxsav_;

      SampleAndPrintOn(strm,x2);
    }

  void SampleAndPrintOn(std::ostream& strm, const stepT x2)
    {
      // Sample to x2, using a sampling interval close to current dxsav.
      // Output to strm as we go.

      int n_sample = (int)(std::floor(std::abs((x2 - x_lastsav)/dxsav) + 0.5));

      // Adjust the magnitude and sign of dxsav.
      dxsav = (x2 - x_lastsav)/n_sample;

      SampleAndPrintOn(strm,n_sample);
    }

  void PrintOn(std::ostream& strm, stepT x0, stepT x1)
    {
      using jlt::operator<<;

      for (mapCIter it = y_data.lower_bound(x0);
	   it != y_data.upper_bound(x1); ++it)
	{
	  strm << (*it).first << "\t" << (*it).second  << std::endl;
	}
    }

  void PrintOn(std::ostream& strm, stepT x0)
    {
      using jlt::operator<<;

      mapCIter it = y_data.lower_bound(x0);

      strm << (*it).first << "\t" << (*it).second  << std::endl;
    }

  void PrintOn(std::ostream& strm)
    {
      using jlt::operator<<;

      strm << y_data;
    }

  //
  // Iterators
  //

  // Only the constant iterator is publicly available.
  typedef mapCIter const_iterator;

  const_iterator begin() const { return y_data.begin(); }
  const_iterator end() const { return y_data.end(); }

}; // class DataPoints

} // namespace rodent

#endif // RODENT_DATA_HPP
