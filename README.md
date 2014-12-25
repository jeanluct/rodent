# rodent

*rodent* (the Rapid ODE iNTegrator) is a C++ library for integrating ordinary differential equations (ODEs).  The function to be integrated is passed as a template parameter, rather than a function pointer, to allow the compiler to aggressively inline.

### contributors

*rodent* was written and is maintained by [Jean-Luc Thiffeault][1].

### documentation

The header files in the rodent folder are standalone with no associated .cpp files, so no special linking is required.  Examples of usage are in the [testsuite][2] folder.  To compile the examples, run `autoconf; ./configure` from the base folder, then `cd testsuite; make` to compile some examples.

To see if inlining is working properly, type `make inline_test` from the testuite folder.

### license

*rodent* is released under the [MIT License][3].

[1]: http://www.math.wisc.edu/~jeanluc/
[2]: https://github.com/jeanluct/rodent/raw/master/testsuite
[3]: https://github.com/jeanluct/rodent/raw/master/LICENSE
