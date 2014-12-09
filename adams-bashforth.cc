#include "adams-bashforth.h"
#include "model.h"

AdamsBashforth::AdamsBashforth(double dt, const Model &model, Integrator *integrator_input)
    : dimen_(model.dimen()),
      dt_(dt),
      model_(model),
      first_step_integrator(integrator_input) {
  fx_ = new double[dimen_];
  prevk= new double[dimen_];
}

AdamsBashforth::~AdamsBashforth() {
  delete [] fx_;
}

int AdamsBashforth::Step(double t, double *x) {


  double h=dt_;
  
   
  static int firststep=1;	
// First step handled separately	     
  if(firststep==1) {
    firststep=0;
    model_.rhs(t,x,fx_);
    first_step_integrator->Step(t,x);
    for(int i=0;i<dimen_;i++){
    prevk[i]=fx_[i];
    }
    delete first_step_integrator;
  } 
  else {
  model_.rhs(t,x,fx_);
 
  for(int i=0;i<dimen_;i++) {
    
    x[i]+=(h*1.5)*fx_[i] - (h*0.5)*prevk[i];
  }
    
  for(int i=0;i<dimen_;i++){
    prevk[i]=fx_[i];
   }
  } 
  return 0;

}
