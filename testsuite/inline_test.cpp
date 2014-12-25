//
// Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
//
// See the file LICENSE for copying permission.
//

#include <vector>
#include <rodent/explicitrk.hpp>
#include "simpharmonic.hpp"


//
// Test inlining of RHS function.
//
// Compile this program (to object with -c) with/without the -O flag.
// Then run
//
//   objdump -Cd inline_test.o | grep SimpleHarmonic\<double\>::operator\(\)
//
// and see if the function is called or not.
//

std::vector<double> sho(double t)
{
  typedef SimpleHarmonic<double> Flow;
  typedef rodent::AdaptiveRK4<Flow> Flow_int;

  Flow flow(1.);
  std::vector<double> y(flow.size());

  double acc = 1.e-10;

  // Initial conditions
  y[0] = 0.;
  y[1] = 1.;

  Flow_int flow_int(flow);
  flow_int
    .tolerance(acc)
    .setState(0.0,y);

  return flow_int(t);
}
