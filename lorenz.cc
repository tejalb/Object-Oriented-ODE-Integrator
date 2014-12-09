#include "lorenz.h"
#include <math.h>

Lorenz::Lorenz(double sigma, double rho, double beta)
    : sigma_(sigma),
      rho_(rho),
      beta_(beta)
{}

Lorenz::~Lorenz()
{}

int Lorenz::rhs(double t, const double *x, double *fx) const {
  fx[0] = sigma_*(x[1]-x[0]);
  fx[1] = rho_*x[0] -x[1] -x[0]*x[2];
  fx[2] = -beta_*x[2] + x[0]*x[1];
  return 0;
}
