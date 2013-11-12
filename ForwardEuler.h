#ifndef FORWARDEULER_H
#define FORWARDEULER_H
#include "Integrator.h"

class ForwardEuler: public Integrator {
 public:
 ForwardEuler(fn f, const double dt = 1e-3) : Integrator(f, dt) {};
  ~ForwardEuler() {};
  vec step(const vec x) {
    vec k = f(x);
    vec new_x;
    for(int i = 0; i < x.size(); i++) {
      new_x.push_back(x[i] + dt*k[i]);
    }
    return new_x;
  };
};

#endif
