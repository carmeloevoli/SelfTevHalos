#include "utilities/misc.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#define index(i, j) ((j) + (i)*Y.size())

namespace utils {

double interpolate(double x, const std::vector<double> &X, const std::vector<double> &Y) {
  std::vector<double>::const_iterator it = std::upper_bound(X.begin(), X.end(), x);
  if (it == X.begin()) return Y.front();
  if (it == X.end()) return Y.back();

  size_t i = it - X.begin() - 1;
  return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
}

double interpolate2d(double x, double y, const std::vector<double> &X, const std::vector<double> &Y,
                     const std::vector<double> &Z) {
  std::vector<double>::const_iterator itx = std::upper_bound(X.begin(), X.end(), x);
  std::vector<double>::const_iterator ity = std::upper_bound(Y.begin(), Y.end(), y);

  if (x > X.back() || x < X.front()) return 0;
  if (y > Y.back() || y < Y.front()) return 0;

  if (itx == X.begin() && ity == Y.begin()) return Z.front();
  if (itx == X.end() && ity == Y.end()) return Z.back();

  size_t i = itx - X.begin() - 1;
  size_t j = ity - Y.begin() - 1;

  double Q11 = Z[index(i, j)];
  double Q12 = Z[index(i, j + 1)];
  double Q21 = Z[index(i + 1, j)];
  double Q22 = Z[index(i + 1, j + 1)];

  double R1 = ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q11 + ((x - X[i]) / (X[i + 1] - X[i])) * Q21;
  double R2 = ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q12 + ((x - X[i]) / (X[i + 1] - X[i])) * Q22;

  return ((Y[j + 1] - y) / (Y[j + 1] - Y[j])) * R1 + ((y - Y[j]) / (Y[j + 1] - Y[j])) * R2;
}

double interpolateEquidistant(double x, double lo, double hi, const std::vector<double> &Y) {
  if (x <= lo) return Y.front();
  if (x >= hi) return Y.back();

  double dx = (hi - lo) / (Y.size() - 1);
  double p = (x - lo) / dx;
  size_t i = floor(p);
  return Y[i] + (p - i) * (Y[i + 1] - Y[i]);
}

size_t closestIndex(double x, const std::vector<double> &X) {
  size_t i1 = std::lower_bound(X.begin(), X.end(), x) - X.begin();
  if (i1 == 0) return i1;
  size_t i0 = i1 - 1;
  if (std::fabs(X[i0] - x) < std::fabs(X[i1] - x))
    return i0;
  else
    return i1;
}

void printVector(std::vector<double> const &input) {
  std::for_each(input.begin(), input.end(), [](const auto &e) { std::cout << e << " "; });
  std::cout << "\n";
}

// size_t find_nearest(const std::vector<double>& x, const double& x0) {
//   double min_distance = 1e100;
//   size_t min_value = 0;
//   for (size_t i = 0; i < x.size(); ++i)
//     if (fabs(x[i] - x0) < min_distance) {
//       min_distance = fabs(x[i] - x0);
//       min_value = i;
//     }
//   return min_value;
// }

}  // namespace utils