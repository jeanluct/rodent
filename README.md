# rodent

*rodent* (the Rapid ODE iNTegrator) is a C++ library for integrating ordinary differential equations (ODEs).  The function to be integrated is passed as a template parameter, rather than a function pointer, to allow the compiler to aggressively inline.

### contributors

*rodent* was written and is maintained by [Jean-Luc Thiffeault][1].

### documentation

The header files in the rodent folder are standalone with no associated .cpp files, so no special linking is required.  Examples of usage are in the [testsuite][2] folder.  When compiling, the root folder of the *rodent* project must be in the include path.  For example, you can compile the program [simpharmonic.cpp][3] in the [testsuite][2] folder with
```
g++ -O -I.. -o simpharmonic simpharmonic.cpp
```
To compile all the examples, run `make` from the [testsuite][2] folder.

To see if inlining is working properly, type `make inline_test` from the [testsuite][2] folder:
```Shell
> make inline_test

Compile without inlining:

g++ -Wall -c -I.. -I../extern/jlt inline_test.cpp

Disassembly shows the right-hand side function is called:

objdump -Cd inline_test.o \
		| grep --color SimpleHarmonic\<double\>::operator\(\)
0000000000000000 <SimpleHarmonic<double>::operator()(double, std::vector<double, std::allocator<double> > const&, std::vector<double, std::allocator<double> >&)>:
  26:	e8 00 00 00 00       	callq  2b <SimpleHarmonic<double>::operator()(double, std::vector<double, std::allocator<double> > const&, std::vector<double, std::allocator<double> >&)+0x2b>
  3a:	e8 00 00 00 00       	callq  3f <SimpleHarmonic<double>::operator()(double, std::vector<double, std::allocator<double> > const&, std::vector<double, std::allocator<double> >&)+0x3f>
  51:	e8 00 00 00 00       	callq  56 <SimpleHarmonic<double>::operator()(double, std::vector<double, std::allocator<double> > const&, std::vector<double, std::allocator<double> >&)+0x56>
  62:	f2 0f 10 05 00 00 00 	movsd  0x0(%rip),%xmm0        # 6a <SimpleHarmonic<double>::operator()(double, std::vector<double, std::allocator<double> > const&, std::vector<double, std::allocator<double> >&)+0x6a>
  7f:	e8 00 00 00 00       	callq  84 <SimpleHarmonic<double>::operator()(double, std::vector<double, std::allocator<double> > const&, std::vector<double, std::allocator<double> >&)+0x84>

Compile with inlining (-O):

g++ -Wall -c -O -I.. -I../extern/jlt inline_test.cpp
Function not found in objdump (good)!

The right-hand side function does not appear in the disassembly.
```
The first compilation, without inlining, shows that the function call to the right-hand side function `SimpleHarmonic<double>::operator()` does indeed show up in the object dump.  In the second, with inlining turned on, the `grep` command finds no call to the function.  This is possible because the right-hand side function is passed to the integrator as a template argument rather than a function pointer, so is available at compile time.  This also means that a new integrator is created for each function to be integrated.

### license

*rodent* is released under the [MIT License][4].

[1]: http://www.math.wisc.edu/~jeanluc/
[2]: https://github.com/jeanluct/rodent/raw/master/testsuite
[3]: https://github.com/jeanluct/rodent/raw/master/testsuite/simpharmonic.cpp
[4]: https://github.com/jeanluct/rodent/raw/master/LICENSE

[![Analytics](https://ga-beacon.appspot.com/UA-58116885-1/rodent/readme)](https://github.com/igrigorik/ga-beacon)
