#include <vector>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>

typedef std::vector< double > vect;
const double PI = 3.14159;
const int STARTING_RES = 8;

void jacobi(double** A, vect &x, const double* b, const int n) {
  int i, j;
  double temp[n];
  for(i = 0; i < n; i++) {
    temp[i] = x[i];
  }
  for(i = 0; i < n; i++) {
    double sum = 0;
    for(j = 0; j < n; j++) {
      sum += A[i][j]*temp[j];
    }
    sum -= A[i][i]*temp[i];
    x[i] = (b[i] - sum)/A[i][i];
  }
}

double* smoothen(vect& x, const int j) {
  int n = pow(2, (j+1)) - 1;
  double* newx = new double[n];
  int i;
  for(i = 0; i < n; i++) {
    if(i != 0 && i != n-1) {
      if((i+1) % 2 == 0) {
	newx[i] = x[(i-1)/2];
      }
      else {
	newx[i] = (x[(i-2) / 2] + x[i / 2])/2;
      }
    }
    else if(i == 0) {
      newx[i] = x[i];
    }
    else if(i == n-1) {
      newx[i] = x[i];
    }
  }
  return newx;
}

double* coarsen(vect &resid, const int j) {
  int n = pow(2, (j-1)) - 1;
  double* coarse_resid = new double[n];
  int i;
  for(i = 0; i < n; i++) {
    coarse_resid[i] = (resid[2*i] + resid[2*i+1] + resid[2*i+2])/3;
  }
  return coarse_resid;
}

void save_data(const std::vector< vect > data, std::ofstream &filehandle) {
  int n = data[0].size();
  int m  = data.size();
  int i, j;
  for(i = 0; i < n; i++) {
    for(j = 0; j < m-1; j++) {
      filehandle << data[j][i] << ",";
    }
    filehandle << data[j][i] << std::endl;
  }
}

void init_eqn(double **A, const int n) {
  int i, j;
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      A[i][j] = 0;
    }
    if(i != 0 && i != n-1) {
      A[i][i] = -2;
      A[i][i-1] = 1;
      A[i][i+1] = 1;
    }
    else if(i == 0) {
      A[i][i] = 1;
    }
    else if(i == n-1) {
      A[i][i] = 1;
    }
  }
}

vect get_resid(double **A, vect x, const double *b, const int n) {
  vect neg_resid(n);
  for(int i = 0; i < n; i++) {
    double sum = 0;
    for(int j = 0; j < n; j++) {
      sum += A[i][j]*x[j];
    }
    neg_resid[i] = b[i] - sum;
  }
  return neg_resid;
}

void multigrid(double **A, std::vector< vect > &xs, const double *b, int n, int depth, const int maxdepth) {
  int i;
  int njacobi_iters = 10;
  int level = STARTING_RES - depth;
  vect x = xs[depth];
  double** newA;
  double *newb;
  int newn;
  for(i = 0; i < njacobi_iters; i++) {
    jacobi(A, x, b, n);
  }
  xs[depth] = x;
  //if not at bottom, coarsen and descend
  if(depth < maxdepth - 1) {
    vect resid = get_resid(A, x, b, n);
    newb = coarsen(resid, level);
    newn = pow(2, level - 1) - 1;
    newA = new double*[newn];
    for(int j = 0; j < newn; j++) {
      newA[j] = new double[newn];
    }
    init_eqn(newA, newn);
    xs.push_back(vect(newn, 0));
    multigrid(newA, xs, newb, newn, depth + 1, maxdepth);
  }
  //if at bottom, smoothen and ascend
  if(depth > 0) {
    x = xs[depth];
    for(i = 0; i < njacobi_iters; i++) {
      jacobi(A, x, b, n);
    }
    double* smoothed_x = smoothen(x, level);
    int nextn = pow(2, level + 1) - 1;
    for(i = 0; i < nextn; i++) {
      xs[depth-1][i] += smoothed_x[i];
    }
  }
  else{
    x = xs[depth];
    for(i = 0; i < njacobi_iters; i++) {
      jacobi(A, x, b, n);
    }
    xs[depth] = x;
  }
}

int main(int argc, char *argv[]) {
  int i, j;
  int n = atoi(argv[1]);
  const int nsteps = atoi(argv[2]);
  double h = 1.0/(n-1);
  //init matrices/vectors
  double *calc_x = new double[n];
  double *analytical_x = new double[n];
  double *b = new double[n];
  double **A = new double*[n];
  for(i = 0; i < n; i++) {
    double x = i*h;
    A[i] = new double[n];
    analytical_x[i] = 1 + 12*x - 10*pow(x, 2) + 0.5*sin(20*PI*pow(x, 3));
    b[i] = h*h*(-20 + 120*0.5*PI*x*cos(20*PI*pow(x, 3)) - 0.5*pow(60*PI*pow(x, 2), 2)*sin(20*PI*pow(x, 3)));
    calc_x[i] = 0;
  }
  b[0] = 1;
  b[n-1] = 3;
  init_eqn(A, n);
  std::stringstream ss;
  ss << "data/jacobi_" << n << ".csv";
  std::string filename = ss.str();
  std::ofstream output(filename);
  std::vector< vect > data;
  //straight dope jacobi

  vect vcalc_x(calc_x, calc_x + n);
  for(i = 0; i < nsteps; i++) {
    jacobi(A, vcalc_x, b, n);
    if(i+1 == 20 || i+1 == 100 || i+1 == 1000) {
      data.push_back(vcalc_x);
    }
  }
  data.push_back(vcalc_x);
  data.push_back(vect(analytical_x, analytical_x+n));

  //straight garbage multigrid
  /**
  const int maxdepth = 7;
  std::vector< vect > xs;
  xs.push_back(vect(calc_x, calc_x + n));
  multigrid(A, xs, b, n, 0, maxdepth);
  data.push_back(xs.front());
  **/


  save_data(data, output);
  output.close();
  delete[] calc_x;
  delete[] analytical_x;
  delete[] b;
  for(i = 0; i < n; i++) {
    delete[] A[i];
  }
  delete[] A;
}
