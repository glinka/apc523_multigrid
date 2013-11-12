#include <cstdlib>
#include <sstream>
#include <fstream>
#include <cmath>
#include "ForwardEuler.h"
#include "ModifiedEuler.h"
#include "RK4.h"


vec my_fn(vec x) {
  vec eval;
  eval.push_back(x[2]);
  eval.push_back(x[3]);
  eval.push_back(x[0]/pow((pow(x[0], 2) + pow(x[1], 2)), 2));
  eval.push_back(x[1]/pow((pow(x[0], 2) + pow(x[1], 2)), 2));
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
  x.push_back(1);
  x.push_back(1);
  x.push_back(0);
  x.push_back(0);
  x.push_back(0);
  for(long int i = 0; i < nsteps; i++) {
    x = alg->step(x);
    saved_data[i] = x;
  }
  save_data(saved_data, filename);
  delete alg;
}
