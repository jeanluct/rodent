#
# Copyright (c) 1999-2014 Jean-Luc Thiffeault <jeanluc@mailaps.org>
#
# See the file LICENSE for copying permission.
#

srcdir = .

RODENTDIR = ..
JLTDIR = ../extern/jlt

CXXFLAGS = -Wall -O -I$(srcdir)/$(RODENTDIR) -I$(JLTDIR)

LIBS = -lm

############################################################################

TESTSUITE = \
	simpharmonic \
	cexp \
	adaptiverk4 \
	cadaptiverk4 \
	valadaptiverk4 \
	arkCashKarp \
	stiff \
	order \
	adamsbashforth \
	ADrwave

.SUFFIXES: .o .cpp

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $(srcdir)/$*.cpp

all:		$(TESTSUITE)
tests:		$(TESTSUITE)
testsuite:	$(TESTSUITE)

$(TESTSUITE):
	$(CXX) $(LDFLAGS) $@.o -o $@ $(LIBS)

###########################################################################

simpharmonic:		simpharmonic.o
cexp:			cexp.o
adaptiverk4:		adaptiverk4.o
cadaptiverk4:		cadaptiverk4.o
valadaptiverk4:		valadaptiverk4.o
arkCashKarp:		arkCashKarp.o
stiff:			stiff.o
order:			order.o
adamsbashforth:		adamsbashforth.o
ADrwave:		ADrwave.o

inline_test:
	@rm -f inline_test.o
	@echo
	@echo 'Compile without inlining:'
	@echo
	$(CXX) -Wall -c -I$(srcdir)/$(RODENTDIR) -I$(JLTDIR) inline_test.cpp
	@echo
	@echo 'Disassembly shows the right-hand side function is called:'
	@echo
	objdump -Cd inline_test.o | grep --color Harmonic\<double\>::operator\(\)
	@echo
	@echo 'Compile with inlining (-O):'
	@echo
	$(CXX) -Wall -c -O -I$(srcdir)/$(RODENTDIR) -I$(JLTDIR) inline_test.cpp
	@-objdump -Cd inline_test.o \
		| grep Harmonic\<double\>::operator\(\) \
		&& ([ $$? -eq 0 ] && echo "Oops... function found.") \
		|| echo "Function not found in objdump (good)!"
	@echo
	@echo 'The right-hand side function does not appear in the disassembly.'

###########################################################################

clean:
	rm -f *.o $(TESTSUITE)
