#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <vector>
#include <ctime>
#include <stdexcept>
#include <fstream>
#include <algorithm>

#include <Eigen/Dense>

extern "C" {
void minimize_linear_square_divergence(
    double rho,
    double *z,
    int z_size,
    double mass,
    double scale,
    double *out);
}

using namespace Eigen;
static const double TOL = 0.000001;

double solve_inner_lambda_eta(double lambda_, const std::vector<double> &y,
			      int z_size, double mass) {
  double cumulative_sum = 0.0;
  double eta;
  double next_eta = 0.0;

  for (int i = 0; i < z_size; ++i) {
    cumulative_sum += y[i];
    eta = next_eta;
    next_eta = (-cumulative_sum - (z_size * mass - i - 1) * lambda_) / (i + 1);
    if (lambda_ < y[i] + next_eta) {
      return eta;
    }
  }
  return next_eta;
}

void minimize_linear_square_divergence(
    double rho, double *z, int z_size, double mass,
    double scale, double *out) {
  double lambda_min = 0.0;
  double lambda_max = 1.0;
  double lambda_ = 0.5;
  double eta = 0;
  Map< VectorXd > z_eigen(z, z_size);

  std::vector<double> y(z,z+z_size);
  int min_ind = std::min_element(y.begin(), y.end()) - y.begin();
  std::sort(y.begin(), y.end());

  Map< VectorXd > x(out, z_size);
  x[min_ind] = 1.0;

  const double compare_value =
    2 * (mass + rho) / scale - z_size / std::pow(scale, 2);

  if (x.squaredNorm() <= compare_value) {
    return;
  }
  while (x.squaredNorm() > compare_value) {
    lambda_max *= 2;
    eta = solve_inner_lambda_eta(lambda_max, y, z_size, mass);
    // Eigen won't let you do a scalar - vector
    x = (((-z_eigen).array() + (lambda_max - eta)) / (lambda_max * z_size)).matrix().cwiseMax(0.0);
  }
  while (abs(lambda_max - lambda_min) > TOL) {
    lambda_ = (lambda_max + lambda_min) / 2;
    eta = solve_inner_lambda_eta(lambda_, y, z_size, mass);
    // Eigen won't let you do a scalar - vector
    x = (((-z_eigen).array() + (lambda_ - eta)) / (lambda_ * z_size)).matrix().cwiseMax(0.0);
    if (x.squaredNorm() > compare_value) {
      lambda_min = lambda_;
    } else {
      lambda_max = lambda_;
    }
  }
}
