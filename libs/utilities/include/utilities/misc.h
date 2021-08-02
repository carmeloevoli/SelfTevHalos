#ifndef UTILITIES_H
#define UTILITIES_H

#include <algorithm>
#include <cmath>
#include <iostream>
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
template <typename T>
void printVector(std::vector<T> const& input) {
  std::for_each(input.begin(), input.end(), [](const auto& e) { std::cout << e << " "; });
  std::cout << "\n";
}

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

template <typename T>
std::vector<T> buildLogAxis(const T& min, const T& max, const size_t& size) {
  using std::exp;
  using std::log;
  const T delta_log = exp(log(max / min) / (size - 1));
  std::vector<T> axis(size);
  int i = 0;
  generate(axis.begin(), axis.end(), [min, delta_log, &i]() { return exp(log(min) + (T)(i++) * log(delta_log)); });
  return axis;
}

template <typename T>
std::vector<T> buildLinAxis(const T& min, const T& max, const size_t& size) {
  const T dx = (max - min) / (T)(size - 1);
  std::vector<T> axis(size);
  int i = 0;
  generate(axis.begin(), axis.end(), [min, dx, &i]() { return min + dx * (T)(i++); });
  return axis;
}

double NIntegrate(std::vector<double> f, double h);

template <typename T>
void infoAxis(std::vector<T> const& input, const std::string& name, const T& units) {
  const auto min = *min_element(input.begin(), input.end()) / units;
  const auto max = *max_element(input.begin(), input.end()) / units;
  const auto size = input.size();
  std::cout << name << " axis in : " << min << " ... " << max << " with " << size << " elements.\n";
}

}  // namespace utils

#endif