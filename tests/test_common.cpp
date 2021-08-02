#include "common.h"
#include "utilities/units.h"

using namespace CRWAVES;

int main() {
  std::cout << larmor_radius(cgs::PV, cgs::muG) / cgs::parsec << "\n";
  std::cout << larmor_radius(1e-3 * cgs::GV, cgs::gauss) / cgs::cm << "\n";
  std::cout << resonant_momentum(1. / cgs::parsec, cgs::muG) / cgs::PV << "\n";
  std::cout << gamma2(1e-5 * cgs::GV) << " " << gamma2(1e5 * cgs::GV) << "\n";
  std::cout << beta(1e-5 * cgs::GV) << " " << beta(1e5 * cgs::GV) << "\n";
  std::cout << alfvenSpeed(1. / cgs::cm3, cgs::muG) / (cgs::km / cgs::second) << "\n";
}