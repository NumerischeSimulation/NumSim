#pragma once

#include "pressure_solver/0_pressure_solver.h"
#include "data_transfer/data_transfer.h"

#include <cmath>

#include <mpi.h>

class RedBlack: public PressureSolver
{
public:
  //! use constructor of parent class pressure solver
  using PressureSolver::PressureSolver;
  
  //! solve the pressure poisson equation
  void solve();

  double calculateResidual();
  
};
