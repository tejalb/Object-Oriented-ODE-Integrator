#include "runge-kutta.h"
#include "model.h"

RungeKutta4::RungeKutta4(double dt, const Model &model)
    : dimen_(model.dimen()),
      dt_(dt),
      model_(model) {
  fx_ = new double[dimen_];
}

RungeKutta4::~RungeKutta4() {
  delete [] fx_;
}

int RungeKutta4::Step(double t, double *x) {
  
 
  double k2[dimen_];
  double k3[dimen_];
  double k4[dimen_];
  double h=dt_;
  
  
  model_.rhs(t,x,fx_);

  double x1[dimen_],x2[dimen_],x3[dimen_];
  for (int i=0;i<dimen_;i++) {
    x1[i]=x[i]+(h/2.)*fx_[i];
  }
  
  model_.rhs(t+(h/2),x1,k2);

  for (int i=0;i<dimen_;i++) {
    x2[i]=x[i]+(h/2.)*k2[i];
  }

  model_.rhs(t+(h/2.),x2,k3);

  for (int i=0;i<dimen_;i++) {
    x3[i]=x[i]+h*k3[i];
  }

  model_.rhs(t+h,x3,k4);   // k1,k2,k3,k4 are updated 
 
  for(int i=0;i<dimen_;i++) {
    
    x[i]+=(h/6.0)*fx_[i] + (h/3.0)*(k2[i]+k3[i]) + (h/6.0)*k4[i];
  }
  return 0;
}
