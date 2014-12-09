#ifndef LINEAR_H_
#define LINEAR_H_

#include "model.h"

// Forced linear oscillator:
//   \dot x = y
//   \dot y = -x - 2*beta*y -  gamma * cos(omega t)
class Linear : public Model {
 public:
  Linear(double beta, double gamma, double omega);
  ~Linear();
  int rhs(double t, const double *x, double *fx) const;
  int dimen() const { return kDimen; }
 private:
  const double beta_;                        // damping constant
  const double gamma_;                        // amplitude of forcing
  const double omega_;                        // frequency of forcing
  static const int kDimen = 2;
};

#endif  // LINEAR_H_
