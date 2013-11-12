#ifndef RK4_H
#define RK4_H
#include "Integrator.h"

class RK4: public Integrator {
 public:
 RK4(fn f, const double dt = 1e-3): Integrator(f, dt) {};
    ~RK4() {};
  vec step(const vec x) {
    int i;
    vec k1 = f(x);
    vec x2 = x;
    for(i = 0; i < x.size() - 1; i++) {
      x2[i] += k1[i]/2;
    }
    x2[i] += dt/2;
    vec k2 = f(x2);
    vec x3 = x;
    for(i = 0; i < x.size() - 1; i++) {
      x3[i] += k2[i]/2;
    }
    x3[i] += dt/2;
    vec k3 = f(x3);
    vec x4 = x;
    for(i = 0; i < x.size() - 1; i++) {
      x4[i] += k3[i];
    }
    x4[i] += dt;
    vec k4 = f(x4);
    vec new_x;
    for(i = 0;i < x.size(); i++) {
      new_x.push_back(x[i] + dt*(k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6));
    }
    return new_x;
  };
};

#endif
