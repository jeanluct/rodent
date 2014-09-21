# rodent

*rodent* (the Rapid ODE iNTegrator) is a C++ library for integrating ordinary differential equations (ODEs).  The function to be integrated is passed as a template parameter, rather than a function pointer, to allow the compiler to aggressively inline.

### contributors

*rodent* was written and is maintained by [Jean-Luc Thiffeault][1].

### documentation

The header files in the rodent folder are standalone with no associated .cpp files.  To see examples of usage, run "autoconf; ./configure" from the base folder, then "cd testsuite; make" to compile some examples.

### license

*rodent* is released under the [MIT License][2].

[1]: http://www.math.wisc.edu/~jeanluc/
[2]: http://bitbucket.org/jeanluc/jlt/raw/tip/LICENSE
