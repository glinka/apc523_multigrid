#ifndef INTEGRATORS_H
#define INTEGRATORS_H
#include <vector>

typedef std::vector<double> vec;
typedef vec (*fn)(vec);

//in all cases, the vector input must be organized as x1, x2, ..., xn, t

class Integrator {
 protected:
  fn f;
  const double dt;
 public:
  virtual vec step(const vec x) = 0;
 Integrator(fn f, const double dt = 1e-3): f(f), dt(dt) {};
  ~Integrator() {};
};

#endif
