#ifndef LORENZ_H_
#define LORENZ_H_

#include "model.h"

// Lorenz system:
//   \dot x = sigma*(y-x)
//   \dot y = rho*x - y - x*z
//   \dot z = -beta*z + x*y
class Lorenz : public Model {
 public:
  Lorenz(double sigma, double rho, double beta);
  ~Lorenz();
  int rhs(double t, const double *x, double *fx) const;
  int dimen() const { return kDimen; }
 private:
  const double sigma_;                        // damping constant
  const double rho_;                        // amplitude of forcing
  const double beta_;                        // frequency of forcing
  static const int kDimen = 3;
};

#endif  // LORENZ_H_
