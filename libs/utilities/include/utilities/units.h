// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef UNITS_H_
#define UNITS_H_

#include <cmath>

#include "utilities/misc.h"

// https://en.wikipedia.org/wiki/Gaussian_units
namespace cgs {

#define pow2 utils::pow_integer<2>

// basic units
static const double centimeter = 1;   // LENGTH
static const double gram = 1;         // MASS
static const double second = 1;       // TIME
static const double kelvin = 1;       // TEMPERATURE
static const double steradian = 1;    // SOLID ANGLE
static const double dyne = 1;         // FORCE
static const double erg = 1;          // ENERGY
static const double kayser = 1;       // WAVENUMBER
static const double hertz = 1;        // FREQUENCY
static const double statcoulomb = 1;  // ELECTRIC CHARGE
static const double gauss = 1;        // MAGNETIC FLUX DENSITY

// derived units
static const double parsec = 3.0856775807e18 * centimeter;
static const double electronvolt = 1.60217657e-12 * erg;
static const double year = 3.15569e7 * second;

// physical constants
static const double c_light = 2.99792458e10 * centimeter / second;
static const double electron_mass = 9.10938370e-28 * gram;
static const double proton_mass = 1.6726219e-24 * gram;
static const double k_B = 1.380649e-16 * erg / kelvin;
static const double elementary_charge = 4.80320427e-10 * statcoulomb;
static const double sigma_th = 6.6524e-25 * utils::pow_integer<2>(centimeter);

// derived quantities
static const double me_c = electron_mass * c_light;
static const double me_c2 = electron_mass * utils::pow_integer<2>(c_light);
static const double volt = electronvolt / c_light;
static const double cm2 = utils::pow_integer<2>(centimeter);
static const double cm3 = utils::pow_integer<3>(centimeter);

// prefixes
static const double peta = 1e15;
static const double tera = 1e12;
static const double giga = 1e9;
static const double mega = 1e6;
static const double kilo = 1e3;
static const double milli = 1e-3;
static const double micro = 1e-6;
static const double nano = 1e-9;
static const double pico = 1e-12;
static const double meter = 1e2 * centimeter;
static const double kilometer = kilo * meter;
static const double kiloyear = kilo * year;
static const double kiloparsec = kilo * parsec;
static const double megaparsec = mega * parsec;
static const double microgauss = micro * gauss;
static const double gigavolt = giga * volt;
static const double teravolt = tera * volt;
static const double petavolt = peta * volt;

// abbreviations
static const double cm = centimeter;
static const double km = kilometer;
static const double kpc = kiloparsec;
static const double kyr = kiloyear;
static const double muG = microgauss;
static const double eV = electronvolt;
static const double GV = gigavolt;
static const double TV = teravolt;
static const double PV = petavolt;
static const double K = kelvin;

}  // namespace cgs

// namespace mks {

// // MKS units
// static const double meter = 1;
// static const double kilogram = 1;
// static const double second = 1;
// static const double ampere = 1;
// static const double kelvin = 1;
// static const double mole = 1;
// static const double joule = kilogram * utils::pow_integer<2>(meter / second);
// static const double tesla = kilogram / ampere / utils::pow_integer<2>(second);
// static const double coulomb = ampere * second;

// // derived units
// static const double centimeter = 1e-2 * meter;
// static const double kilometer = 1e3 * meter;
// static const double gram = 1e-3 * kilogram;
// static const double erg = 1e-7 * joule;
// static const double gauss = 1e-4 * tesla;
// static const double microgauss = 1e-6 * gauss;
// static const double nanogauss = 1e-9 * gauss;

// // abbreviations
// static const double s = second;
// static const double K = kelvin;
// static const double kg = kilogram;
// static const double cm = centimeter;
// static const double km = kilometer;
// static const double cm2 = cm * cm;
// static const double cm3 = cm * cm * cm;
// static const double m2 = meter * meter;
// static const double m3 = meter * meter * meter;
// static const double muG = microgauss;
// static const double nG = nanogauss;

// // electron volt
// static const double electronvolt = 1.60217657e-19 * joule;
// static const double kiloelectronvolt = 1e3 * electronvolt;
// static const double megaelectronvolt = 1e6 * electronvolt;
// static const double gigaelectronvolt = 1e9 * electronvolt;
// static const double teraelectronvolt = 1e12 * electronvolt;
// static const double petaelectronvolt = 1e15 * electronvolt;
// static const double exaelectronvolt = 1e18 * electronvolt;
// static const double eV = electronvolt;

// // time
// static const double year = 3.15569e7 * second;
// static const double kiloyear = 1e3 * year;
// static const double megayear = 1e6 * year;
// static const double gigayear = 1e9 * year;
// static const double kyr = kiloyear;
// static const double Myr = megayear;
// static const double Gyr = gigayear;

// // parsec
// static const double parsec = 3.0856775807e18 * centimeter;
// static const double kiloparsec = 1e3 * parsec;
// static const double megaparsec = 1e6 * parsec;
// static const double gigaparsec = 1e9 * parsec;
// static const double pc = parsec;
// static const double kpc = kiloparsec;
// static const double Mpc = megaparsec;
// static const double Gpc = gigaparsec;

// // physical constants
// static const double c_light = 2.99792458e10 * centimeter / second;
// static const double c_squared = c_light * c_light;
// static const double electron_mass = 9.10938356e-31 * kilogram;
// static const double electron_mass_c = electron_mass * c_light;
// static const double electron_mass_c2 = electron_mass * c_squared;
// static const double proton_mass = 1.6726219e-27 * kilogram;
// static const double sun_mass = 1.989e30 * kg;
// static const double sigma_th = 6.6524e-29 * utils::pow_integer<2>(meter);
// static const double h_planck = 6.62607004e-34 * utils::pow_integer<2>(meter) * kg / s;
// static const double k_boltzmann = 1.3806488e-23 * joule / kelvin;
// static const double elementary_charge = 1.60217662e-19 * coulomb;
// static const double vacuum_permeability = 4e-7 * M_PI * tesla * meter / ampere;
// static const double electron_radius = 2.8179403227e-15 * meter;

// // derived units with physical constants
// static const double keV_c = kiloelectronvolt / c_light;
// static const double MeV_c = megaelectronvolt / c_light;
// static const double GeV_c = gigaelectronvolt / c_light;
// static const double TeV_c = teraelectronvolt / c_light;
// static const double PeV_c = petaelectronvolt / c_light;
// static const double EeV_c = exaelectronvolt / c_light;

// }  // namespace mks

#endif
