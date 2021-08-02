#include <iomanip>
#include <iostream>

#include "utilities/units.h"

template <size_t N>
double integrate_power_law(double slope) {
  const auto x = utils::buildLogAxis(1., 1e3, N);
  std::vector<double> y;
  for (auto xi : x) {
    y.push_back(xi * std::pow(xi, -slope));
  }
  return utils::NIntegrate(y, std::log(x[1] / x[0]));
}

double integral_value(double slope, double x_min, double x_max) {
  double value = 1. / (1. - slope);
  value *= std::pow(x_max, -slope + 1) - std::pow(x_min, -slope + 1);
  return value;
}

int main() {
  const double slope = 2.4;
  std::cout << std::setprecision(10);
  double I = integral_value(slope, 1., 1e3);
  {
    double NI = integrate_power_law<8>(slope);
    std::cout << NI << " " << (NI - I) / I << "\n";
  }
  {
    double NI = integrate_power_law<16>(slope);
    std::cout << NI << " " << (NI - I) / I << "\n";
  }
  {
    double NI = integrate_power_law<32>(slope);
    std::cout << NI << " " << (NI - I) / I << "\n";
  }
  {
    double NI = integrate_power_law<64>(slope);
    std::cout << NI << " " << (NI - I) / I << "\n";
  }
}