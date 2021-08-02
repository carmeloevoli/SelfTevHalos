// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef TGRID2D_H_
#define TGRID2D_H_

#include <iostream>
#include <string>
#include <vector>

namespace CRWAVES {

template <typename T>
class TGrid2D {
  std::vector<T> grid;
  size_t Nx, Nz;

 public:
  TGrid2D() : Nx(0), Nz(0) {}

  TGrid2D(size_t Nx, size_t Nz) { set_grid_size(Nx, Nz); }

  void set_grid_size(size_t Nx, size_t Nz) {
    this->Nx = Nx;
    this->Nz = Nz;
    grid.resize(Nx * Nz);
  }

  /** Accessor / Mutator */
  T &get(size_t ix, size_t iz) { return grid[ix * Nz + iz]; }

  T &get(const size_t &i) { return grid[i]; }

  /* Min / Max */
  T max_value() const { return *max_element(grid.begin(), grid.end()); }

  T min_value() const { return *min_element(grid.begin(), grid.end()); }

  /** Accessor */
  const T &get(size_t ix, size_t iz) const {
    // assert(ix >= 0 && ix < Nx);
    // assert(iz >= 0 && iz < Nz);
    return grid[ix * Nz + iz];
  }

  const T &get(const size_t &i) const { return grid[i]; }

  T getValue(size_t ix, size_t iz) { return grid[ix * Nz + iz]; }

  /** Return a reference to the grid values */
  std::vector<T> &getGrid() { return grid; }

  void clearGrid() {
    grid.clear();
    return;
  }

  size_t getNx() const { return Nx; }

  size_t getNz() const { return Nz; }

  size_t get_size() const { return (Nx * Nz); }

  void fill(const T &value) { std::fill(grid.begin(), grid.end(), value); }

  void copy(const std::vector<T> &v) { std::copy(v.begin(), v.end(), grid.begin()); }

  void show_grid(const std::string &name, const T &units) const {
    std::cout << name << " grid in : " << min_value() / units << " ... " << max_value() / units << " with "
              << get_size() << " elements."
              << "\n";
  }
};

}  // namespace CRWAVES

#endif  // TGRID2D_H_
