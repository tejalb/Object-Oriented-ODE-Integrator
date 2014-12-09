#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "lorenz.h"
#include "linear-oscillator.h"
#include "duffing.h"
#include "euler.h"
#include "runge-kutta.h"
#include "adams-bashforth.h"
#include <iostream>
using namespace std;

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
  if(argc!=5){
    printf("USAGE: %s <equation> <integrator> <timestep> <numsteps>\n", argv[0]);
    exit(1);
  }
  
  const double dt = atof(argv[3]);
  const int nsteps = atoi(argv[4]);
  const string equation=(argv[1]);
  const string integrator_type = (argv[2]);
  
  // Initial condition
 
  double *x;
 
  cout<<equation<<" "<<integrator_type<<endl;
  Model *model;

  if(equation=="duffing"){
    const double delta = 0.2;
    const double gamma = 0.3;
    const double omega = 1.0;
    model = new Duffing(delta, gamma, omega);
    x = new double[model->dimen()];
    for (int i = 0; i < model->dimen(); ++i) {
      x[i] = 0;
    }
  }
  else if(equation=="lorenz"){  
    const double sigma = 10.;
    const double rho = 28.;
    const double beta = 8./3.;
    model = new Lorenz(sigma, rho, beta);
    x=new double[model->dimen()];
    x[0] = 0;
    x[1]=0.01;
    x[2]=0;
  
  }
  else if(equation=="linear"){
    const double beta = 0.1;
    const double gamma = 1.;
    const double omega = 0.9;
    model = new Linear(beta, gamma, omega);
    x=new double[model->dimen()];
    for (int i = 0; i < model->dimen(); ++i) {
      x[i] = 0;
    }
  }
  else{cout<<"Equation should be 'linear' or 'lorenz' or 'duffing'"<<endl; return 0;}
 Integrator *integrator;

  if(integrator_type=="euler"){
    integrator=  new Euler(dt, *model);
  }
  else if(integrator_type=="rk4"){
    integrator=  new RungeKutta4(dt, *model);
  }
  else if(integrator_type=="ab2"){
    Integrator *first_step_integrator = new RungeKutta4(dt, *model);
    integrator=  new AdamsBashforth(dt, *model, first_step_integrator); //first_step_integrator's scope is only in this block, but 
    //the pointer is copied to AB and the memory allocated on the heap is persistent
  }
  else{cout<<"Integrator should be 'euler' or 'rk4' or 'ab2'"<<endl; return 0;}
 
  // Initial condition
  const int dimen = model->dimen();   // number of states
  double t = 0;
  PrintState(dimen, t, x);
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, x);
    t = (i+1) * dt;
    PrintState(dimen, t, x);
  }

  delete integrator;
  delete x;
  delete model;
  return 0;
}
