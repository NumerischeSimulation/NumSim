#include "pressure_solver.h"

#include <memory>

PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations)

{
  discretization_ = discretization;
  epsilon = epsilon;
  maximumNumberOfIterations_ = maximumNumberOfIterations;

}

  void PressureSolver::setBoundaryValues()
  {
     
      //p_i,-1 = p_i,0. p_i,n+1 = p_i,n.
      for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
      {
        // bottom
        discretization_->p(i, discretization_->pJBegin())  = discretization_->p(i, discretization_->pJBegin() +1);
        // top
        discretization_->p(i, discretization_->pJEnd() -1) = discretization_->p(i, discretization_->pJEnd() -2);
      }

       // prioritise left and right boundaries
       // p_-1,j = p_0,j. p_n+1,j = p_n,j.
      for ( int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
      {
        // left
        discretization_->p(discretization_->pIBegin(),j) =  discretization_->p(discretization_->pIBegin()+1,j);
        // right
        discretization_->p(discretization_->pIEnd() -1,j) =  discretization_->p(discretization_->pIEnd() -2,j);
      }

       
  }     
