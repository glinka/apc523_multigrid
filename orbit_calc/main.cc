#include <ctime>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include "ForwardEuler.h"
#include "ModifiedEuler.h"
#include "RK4.h"


vec my_fn(vec x) {
  vec eval;
  eval.push_back(x[2]);
  eval.push_back(x[3]);
  eval.push_back(-x[0]/pow((pow(x[0], 2) + pow(x[1], 2)), 3./2));
  eval.push_back(-x[1]/pow((pow(x[0], 2) + pow(x[1], 2)), 3./2));
  eval.push_back(1);
  return eval;
}

void save_data(std::vector< vec > data, std::string filename) {
  std::ofstream file(filename);
  for(std::vector< vec >::iterator v = data.begin(); v != data.end(); v++) {
    for(vec::iterator val = (*v).begin(); val != (*v).end() - 1; val++) {
      file << *val << ",";
    }
    file << (*v).back() << std::endl;
  }
  file.close();
}

int main(int argc, char *argv[]) {
  int alg_selection = atoi(argv[1]);
  long int nsteps = atoi(argv[2]);
  double dt = atof(argv[3]);
  std::stringstream ss;
  ss << "data/integration_data_" << alg_selection << "_" << nsteps << "_" << dt << ".csv";
  std::string filename = ss.str();
  Integrator *alg;
  switch(alg_selection) {
  case 1:
    alg = new ForwardEuler(my_fn, dt);
    break;
  case 2:
    alg = new ModifiedEuler(my_fn, dt);
     break;
  case 3:
    alg = new RK4(my_fn, dt);
    break;
  }
  std::vector<vec> saved_data(nsteps);
  vec x;
  x.push_back(1.5);
  x.push_back(1.0);
  x.push_back(-1.0);
  x.push_back(2.0);
  x.push_back(0);
  double r0 = pow((pow(x[0], 2) + pow(x[1], 2)), 0.5);
  double theta0 = acos(x[0]/r0);
  const double h = pow(r0, 2)*((x[0]*x[3]-x[1]*x[2])/(pow(x[0], 2) + pow(x[1], 2)));
  const double e = pow(h, 2)/r0 - 1;
  long int clock_ticks = 0;
  const int STEPS_PER_TIMING = 1000;
  std::clock_t t = clock();
  for(long int i = 0; i < nsteps; i++) {
    x = alg->step(x);
    if((i+1) % STEPS_PER_TIMING == 0) {
      t = clock() - t;
      clock_ticks += t;
      t = clock();
    }
    saved_data[i] = x;
    //calculate analytical values
    double r = pow((pow(x[0], 2) + pow(x[1], 2)), 0.5);
    double theta = acos(x[0]/r);
    r = pow(h, 2)/(1+e*cos(theta-theta0));
    //double dtheta = 2*3.14*i/nsteps
    double x_real = cos(theta)*r;
    double y_real = sin(theta)*r;
    //calculate energy and momentum
    double v = pow((pow(x[2], 2) + pow(x[3], 2)), 0.5);
    double energy = pow(v, 2)/2 - 1/r;
    double ang_momentum = v*r;
    saved_data[i].push_back(x_real);
    saved_data[i].push_back(y_real);
    saved_data[i].push_back(energy);
    saved_data[i].push_back(ang_momentum);
    saved_data[i].push_back(clock_ticks);
  }
  save_data(saved_data, filename);
  delete alg;
}
