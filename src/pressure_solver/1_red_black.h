#pragma once

#include "pressure_solver/0_pressure_solver.h"

#include <cmath>

class RedBlack: public PressureSolver
{
public:
  //! use constructor of parent class pressure solver
  using PressureSolver::PressureSolver;
  
  //! solve the pressure poisson equation
  void solve();
  
};
