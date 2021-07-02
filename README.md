[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-1807.09263-blue.svg)](https://arxiv.org/abs/1807.09263)

# Self-generated cosmic-ray confinement in TeV halos

## About

Code to compute self-consistently the self-generation of Alfven waves by pairs leaving pulsars. The code solves the time-dependent propagation and diffusion of cosmic rays from a time-dependent source in a regime with both time- and distance-dependent diffusion. We applied this code to show that pulsars can produce localized regions where diffusion is inhibited by 2–3 orders of magnitude, consistent with TeV halo observations. 

Read more: [1807.09263](https://arxiv.org/abs/1807.09263)

## Install

Quick procedure for installing the code:

```sh
mkdir build
cd build
cmake ..
make -j
```

It requires the `GSL` C library available at https://www.gnu.org/software/gsl/

