Object-Oriented-ODE-Integrator
==============================

Integrator interface to solve the forced Duffing oscillator, linear oscillator and the lorenz system, using 3 different
integration schemes:

1) 4th order Runge Kutta 
2) 2nd order Adams bashforth
3) Explicit Euler

Description of files:

• model.h: an abstract base class Model for evaluating f(x, t)

• integrator.h: an abstract base class Integrator for ODE integrators

• duffing.{h,cc}: an implementation of the Duffing oscillator (a derived
class of Model).

• euler.{h,cc}: an implementation of an explicit Euler integrator (a derived
class of Integrator).

• duffing_solve.cc: a driver program to integrate the Duffing oscillator
for the given parameter values and initial conditions.

• Makefile: builds the program duffing_solve.

The general driver program ode_solve that lets you choose between
the different equations and integration schemes. It reads commandline
arguments in the following format:

ode_solve {equation} {integrator} {timestep} {numsteps}

The program test_convergence does  convergence test of the different
integrators. The parameters should be specified on the command line
as

test_convergence {integrator} {timestep} {numsteps}

The python script "run_tests" runs test_convergence for all three integrators,
up to time t = 25, for timesteps 0.1, 0.01, 0.001, 0.0001, 0.00001.

