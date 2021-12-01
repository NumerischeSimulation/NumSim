#pragma once

#include "pressure_solver/0_pressure_solver.h"
#include "partitioning/partitioning.h"

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

  //! helper that exchanges the ghost layers for p
  void exchange(int rankCorrespondent, //! rank of the other process with which to exchange 
                int indexToSend, //! index of the slice to send
                int indexFromReceive, //! index of the slice where to write the incoming data
                char direction, //! x or y direction
                bool ToFrom); //! send and then receive, else receive then send

protected:
  std::shared_ptr<Partitioning> partitioning_;
  
};
