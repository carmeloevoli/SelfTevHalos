// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

Waves::Waves(const Params& params) : par(params) {
  std::cout << "... main class constructed\n";
}

Waves::~Waves() {
  std::cout << "... main class destructed\n";
}

}  // namespace CRWAVES