integrators = euler.o runge-kutta.o adams-bashforth.o
equations1 =   duffing.o lorenz.o linear-oscillator.o 
equations2 =   linear-oscillator.o 
objects1 =   $(integrators) $(equations1) ode_solve.o
objects2 =   $(integrators) $(equations2) test_convergence.o

CXXFLAGS = -g -Wall

all: ode_solve test_convergence 

ode_solve : $(objects1)
	$(CXX) -o $@ $^

test_convergence : $(objects2)
	$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
