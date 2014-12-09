#ifndef ADAMSBASHFORTH_H_
#define ADAMSBASHFORTH_H_

#include "integrator.h"
class Model;

class AdamsBashforth : public Integrator {
 public:
  AdamsBashforth(double dt, const Model &model, Integrator *_first_step_integrator);
  ~AdamsBashforth();
  int Step(double t, double *x);
  
  
 private:
  const int dimen_;                     // dimension of state x
  const double dt_;                     // timestep
  const Model &model_;                  // functor to evaluate f(x,t)
  Integrator *first_step_integrator;
 
  double *fx_;                          // temporary space to hold f(x,t)
  double *prevk;                          // temporary space to hold previous f(x,t)
};

#endif  // ADAMSBASHFORTH_H
