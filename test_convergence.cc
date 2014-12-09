#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "linear-oscillator.h"
#include "runge-kutta.h"
#include "adams-bashforth.h"
#include "euler.h"
#include <math.h>
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

inline double sq(double x){
  return x*x;
}

int main(int argc, char *argv[]) {
  /*     if(argc!=4){
   printf("USAGE: %s  <integrator> <timestep> <numsteps>\n", argv[0]);
    return 0;
    }
  */
  const double dt = atof(argv[2]);
  const int nsteps = atoi(argv[3]);
  const string integrator_type = (argv[1]);
  
  // Initial condition
 
  double *x;
  double *xtrue; 
  Model *model;

    const double beta = 0.1;
    const double gamma = 1.;
    const double omega = 0.9;
    model = new Linear(beta, gamma, omega);
    x=new double[model->dimen()];
    xtrue=new double[model->dimen()];
    for (int i = 0; i < model->dimen(); ++i) {
      x[i] = 0;
      xtrue[i]=0;
    }

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
 
  const double a=((1-sq(omega))*gamma)/(sq(1-sq(omega))+4*sq(beta)*sq(omega));
  const double b=(2*beta*omega*gamma)/(sq(1-sq(omega))+4*sq(beta)*sq(omega));
  const double c1=-a;
  const double omega_damped=sqrt(1-pow(beta,2));
  const double c2=(-beta*a - omega*b)/omega_damped;


  // Initial condition
  //  const int dimen = model->dimen();   // number of states
  double t = 0;
  //  PrintState(dimen, t, x);
  double enorm = 0;

  for (int i = 0; i < nsteps; ++i) {
    
    integrator->Step(t, x);
    t = (i+1) * dt;
    xtrue[0]=exp(-beta*t)*(c1*cos(omega_damped*t)+c2*(sin(omega_damped*t))) + a*cos(omega*t) + b*sin(omega*t);
    xtrue[1]=exp(-beta*t)*(-beta*(c1*cos(omega_damped*t)+c2*(sin(omega_damped*t)))+omega_damped*(-c1*sin(omega_damped*t)+c2*cos(omega_damped*t))) -omega*(a*sin(omega*t) - b*cos(omega*t));
    enorm+=(sq(xtrue[0]-x[0])+sq(xtrue[1]-x[1]))*dt;
    //    PrintState(dimen, t, x);
    //    PrintState(dimen, t, xtrue);
  }
  enorm=sqrt(enorm);
  cout<<enorm<<endl;
  delete integrator;
  delete model;
  delete xtrue;
  delete x;
  return 0;
}
