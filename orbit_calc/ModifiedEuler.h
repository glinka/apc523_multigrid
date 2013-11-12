#ifndef MODIFIEDEULER_H
#define MODIFIEDEULER_H
#include "Integrator.h"

class ModifiedEuler: public Integrator {
 public:
 ModifiedEuler(fn f, const double dt = 1e-3): Integrator(f, dt) {};
  ~ModifiedEuler() {};
  vec step(const vec x) {
    int i;
    vec k1 = f(x);
    vec x2 = x;
    for(i = 0 ; i < x.size() - 1; i++) {
      x2[i] += dt*k1[i];
    }
    x2[i] += dt;
    vec k2 = f(x2);
    vec new_x;
    for(int i = 0; i < x.size(); i++) {
      new_x.push_back(x[i] + dt*(k1[i] + k2[i])/2);
    }
    return new_x;
  };
};

#endif
