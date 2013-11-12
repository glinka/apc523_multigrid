#include <vector>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>

const double PI = 3.14159;

void jacobi(double** A, double* x, const double* b, const int n) {
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

void smoothen(double* x, const int j) {
  int n = pow(2, (j+1)) - 1;
  double temp[n];
  int i;
  for(i = 0; i < n; i++) {
    if(i != 0 && i != n-1) {
      if(i % 2 == 0) {
	temp[i] = x[i%2];
      }
      else {
	temp[i] = (x[(i+1) % 2] + x[(i-1) % 2])/2;
      }
    }
    else if(i == 0) {
      temp[i] = x[i];
    }
    else if(i == n-1) {
      temp[i] = x[i];
    }
  }
  delete[] x;
  x = new double[n];
  for(i = 0; i < n; i++) {
    x[i] = temp[i];
  }
}

void coarsen(double* resid, const int j) {
  int n = pow(2, (j-1)) - 1;
  double temp[n];
  int i;
  for(i = 0; i < n; i++) {
    temp[i] = (resid[i-1] + resid[i] + resid[i+1])/3;
  }
  delete[] resid;
  resid = new double[n];
  for(i = 0; i < n; i++) {
    resid[i] = temp[i];
  }
}

void save_data(const std::vector< std::vector<double> > data, std::ofstream &filehandle) {
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

void init_eqn(double **A, double *calc_x, double *b, const int n) {
  double h = 1.0/(n-1);
  int i, j;
  for(i = 0; i < n; i++) {
    double x = i*h;
    b[i] = h*h*(-20 + 120*0.5*PI*x*cos(20*PI*pow(x, 3)) - 0.5*pow(60*PI*pow(x, 2), 2)*sin(20*PI*pow(x, 3)));
    calc_x[i] = 0;
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
      b[i] = 1;
    }
    else if(i == n-1) {
      A[i][i] = 1;
      b[i] = 3;
    }
  }
}

void multigrid(double **A, double *x, const double *b, double* resid, int n, int depth, const int maxdepth) {
  int i;
  int njacobi_iters = 10;
  int level = maxdepth - depth;
  for(i = 0; i < njacobi_iters; i++) {
    jacobi(A, x, b, n);
  }
  resid = get_neg_resid(A, x, b, n);
  if(depth < maxdepth - 1) {
    coarsen(resid, level);
    int newn = pow(2, level - 1) - 1;
    double newA = new double*[newn];
    for(int j = 0; j < newn; j++) {
      newA[j] = new double[newn];
    }
    double newx[newn];
    double newb[newn];
    init_eqn(newA, newx, newb, newn);
    multigrid(newA, newx, resid, newn, depth + 1, maxdepth);
  }
  else if(depth == maxdepth - 1) {
    for(i = 0; i < njacobi_iters; i++) {
      jacobi(A, x, b, n);
    }
    smoothen(x, level);
  }


  for(i = 0; i < n; i++) {
    
  for(i = 0; i < njacobi_iters; i++) {
    jacobi(A, x, b, n);
  }
  smoothen(x, level);
  }
  //pass resid through all levels, add to individual xs
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
    A[i] = new double[n];
    analytical_x[i] = 1 + 12*x - 10*pow(x, 2) + 0.5*sin(20*PI*pow(x, 3));
  }
  init_eqn(A, calc_x, b, n);
  std::stringstream ss;
  ss << "data/jacobi_" << n << ".csv";
  std::string filename = ss.str();
  std::ofstream output(filename);
  std::vector< std::vector<double> > data;
  //straight dope jacobi
  /**
  for(i = 0; i < nsteps; i++) {
    jacobi(A, calc_x, b, n);
    if(i+1 == 20 || i+1 == 100 || i+1 == 1000) {
      data.push_back(std::vector<double>(calc_x, calc_x+n));
    }
  }
  **/  
  const int maxdepth = 7;
  multigrid(A, calc_x, b, n, 0, 7);
  data.push_back(std::vector<double>(calc_x, calc_x+n));
  data.push_back(std::vector<double>(analytical_x, analytical_x+n));
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

  /**
  for(j = 0; j < n; j++) {
    A[0][j] /= A[0][0];
  }
  b[0] /= A[0][0];
  //row-diag
  for(i = 1; i < n; i++) {
    //pivot
    double max = A[i][i];
    int max_row = i;
    for(k = i; k < n; k++) {
      if(A[k][i] > max) {
	max = A[k][i];
	max_row = k;
      }
    }
    for(k = 0; k < n; k++) {
      temp[k] = A[i][k];
      A[i][k] = A[max_row][k];
      A[max_row][k] = temp[k];
    }
    for(j = i-1; j < n; j++) {
      A[i][j] -= A[i][i-1]*A[i-1][j];
      A[i][j] /= A[i][i];
    }
  }
  **/
  
