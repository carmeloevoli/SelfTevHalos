// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef UTILITIES_H_
#define UTILITIES_H_

namespace CRWAVES {

double beta(const double& p);
double larmor_radius(const double& p, const double& B);
double inverse_larmor_radius(const double& k, const double& B);
double I_of_alpha(const double& alpha, const double& x_min, const double& x_cutoff);
double source_evolution(const double& t_now, const double& t_decay);
double source_profile_1D(const double& z, const double& size);
double source_profile_3D(const double& z, const double& size);
double source_spectrum(const double& p, const double& alpha, const double& p_min, const double& p_cutoff);

}  // namespace CRWAVES

#endif /* UTILITIES_H_ */
