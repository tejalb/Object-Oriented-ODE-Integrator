#ifndef RUNGEKUTTA4_H_
#define RUNGEKUTTA4_H_

#include "integrator.h"
class Model;

class RungeKutta4 : public Integrator {
 public:
  RungeKutta4(double dt, const Model &model);
  ~RungeKutta4();
  int Step(double t, double *x);
 private:
  const int dimen_;                     // dimension of state x
  const double dt_;                     // timestep
  const Model &model_;                  // functor to evaluate f(x,t)
  double *fx_;                          // temporary space to hold f(x,t)
};

#endif  // RUNGEKUTTA_H
