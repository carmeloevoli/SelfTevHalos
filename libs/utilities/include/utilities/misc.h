#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>

namespace utils {

// Perform linear interpolation on a set of n tabulated data points X[0 .. n-1] -> Y[0 .. n-1]
// Returns Y[0] if x < X[0] and Y[n-1] if x > X[n-1]
double interpolate(double x, const std::vector<double>& X, const std::vector<double>& Y);

// Perform bilinear interpolation on a set of (n,m) tabulated data points X[0 .. n-1], Y[0 .. m-1] -> Z[0.. n-1*m-1]
// Returns 0 if x < X[0] or x > X[n-1] or y < Y[0] or y > Y[m-1]
double interpolate2d(double x, double y, const std::vector<double>& X, const std::vector<double>& Y,
                     const std::vector<double>& Z);

// Perform linear interpolation on equidistant tabulated data
// Returns Y[0] if x < lo and Y[n-1] if x > hi
double interpolateEquidistant(double x, double lo, double hi, const std::vector<double>& Y);

// Find index of value in a sorted vector X that is closest to x
size_t closestIndex(double x, const std::vector<double>& X);

// Print a vector
void printVector(std::vector<double> const& input);

// pow implementation as template for integer exponents pow_integer<2>(x)
// evaluates to x*x
template <unsigned int exponent>
inline double pow_integer(double base) {
  return pow_integer<(exponent >> 1)>(base * base) * (((exponent & 1) > 0) ? base : 1);
}

template <>
inline double pow_integer<0>(double base) {
  return 1;
}

}  // namespace utils

#endif