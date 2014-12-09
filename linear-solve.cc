#include <stdio.h>
#include <stdlib.h>
#include "linear-oscillator.h"
#include "euler.h"
#include "runge-kutta.h"
#include "adams-bashforth.h"
// Print a line
//    time x[0] x[1] ...
// to standard out
void PrintState(double n, double t, const double *x) {
  printf("%15.8f", t);
  for (int i = 0; i < n; ++i) {
    printf("%15.8f", x[i]);
  }
  printf("\n");
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("USAGE: %s <timestep> <numsteps>\n", argv[0]);
    exit(1);
  }
  const double dt = atof(argv[1]);
  const int nsteps = atoi(argv[2]);

  const double beta = 0.1;
  const double gamma = 1.;
  const double omega = 0.9;
  Model *model = new Linear(beta, gamma, omega);

  // Initial condition
  const int dimen = model->dimen();   // number of states
  double x[dimen];
  for (int i = 0; i < dimen; ++i) {
    x[i] = 0;
  }

  // Set up integration scheme
  Integrator *integrator = new Euler(dt, *model);

  //Integrator *first_step_integrator = new RungeKutta4(dt, *model);
  //Integrator *integrator = new AdamsBashforth(dt, *model, first_step_integrator );

  double t = 0;
  PrintState(dimen, t, x);
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, x);
    t = (i+1) * dt;
    PrintState(dimen, t, x);
  }

  delete integrator;
  delete model;
  return 0;
}
