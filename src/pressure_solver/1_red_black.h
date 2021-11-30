#pragma once

#include "pressure_solver/0_pressure_solver.h"
#include "data_transfer/data_transfer.h"

#include <cmath>
#include <math.h>
#include <mpi.h>

class RedBlack: public PressureSolver
{
public:
  //! use constructor of parent class pressure solver
  RedBlack(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning);
  
  //! solve the pressure poisson equation
  void solve();

  double calculateResidual();

  //! sets the ghost layers in the discretization correctly
  void pExchangeVertical();
  void pExchangeHorizontal();

};
