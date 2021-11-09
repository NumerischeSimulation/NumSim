pragma once

#include "computation.h"
#include "pressure_solver/pressure_solver.h"
#include "discretization/1_discretization.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"

#include <iostream>  // for cout
#include <memory>

void Computation::initialize(int argc, char *argv[])
{
    // parse the parameters
    settings_= Settings::Settings();
    settings_.loadFromFile(argv[1])
    settings_.printSettings();

    // calculate
    std::array<double, 2> meshWidth_;
    for (int i = 0; i++; i < 2)
    {
        meshWidth_[i] =  settings_.physicalSize[i] / settings_.nCells[i];
    }

    // initialize
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }
    
    if (settings_.pressureSolver == "SOR")
    {
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else if (settings_.pressureSolver == "GaussSeidel")
    {
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    } else 
    {
        std::cout << "The name of the pressure solver is not understood" << settings_.pressureSolver << std::endl;
        throw;
    }   

    // misc
    double dt_ = settings_.maximumDt;
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    outputWriterText_ = std::make_unique<OutputWriter>(discretization_);
}

void Computation::runSimulation()
{
    double currentTime = 0;

    // the steps correspond to the steps in our algorithm in the overleaf or docs/numsim-algos.tex
    // step 1: set the boundary values
    applyBoundaryValues();

    while (currentTime < settings_.endTime)
    {
        // step 2: compute time step width
        computeTimeStepWidth();

        // step 4: calculate F, G with first setting the boundary conditions of F, G (step 3)
        computePreliminaryVelocities()

        // step 5: compute the right hand side of the pressure equation
        computeRightHandSide()

        // step 6: solve the pressure equation
        computePressure()

        // step 7: calculate the final velocities
        computeVelocities()

        // step 8: reset boundary values
        applyBoundaryValues(); 
        // may be wrong in the reference solution 
        // should then be set at the beginning of the iteration only

        // step 9: write output
        currentTime += dt_;
        outputWriterParaview_->writeFile(currentTime);
        outputWriterText_->writeFile(currentTime);
    }

    // write output
    outputWriterParaview_->writePressureFile();
    outputWriterText_->writePressureFile();
}

void Computation::computeTimeStepWidth()
{
    // boundary from diffusion
    if (meshWidth_[0] == meshWidth_[1])
    {
        double boundary_diffusion = (settings_.re * meshWidth_[0] * meshWidth_[1]) / 4;
    } else
    {
        double h2x = meshWidth_[0] * meshWidth_[0];
        double h2y = meshWidth_[1] * meshWidth_[1];
        double boundary_diffusion = (settings_.re / 2) * (h2x * h2y) * (1/(h2x + h2y));
    }

    // boundary from convection
    double u_max = std::max(std::abs(discretization_->u.min()), std::abs(discretization_->u.max()));
    double v_max = std::max(std::abs(discretization_->v.min()), std::abs(discretization_->v.max()));
    double boundary_convection_u = meshWidth_[0] / u_max;
    double boundary_convection_v = meshWidth_[1] / v_max;

    // together
    double min_dt = std::min({boundary_diffusion, boundary_convection_u, boundary_convection_v});

    // security factor
    dt_ =  std::min({min_dt * settings_.tau, settings_.maximumDt});
}
    
void Computation::applyBoundaryValues()
{
    // bottom and top
    // set u
    for ( int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        // bottom
        discretization_->u(i, discretization_->uJBegin() +1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() +2);
        
        // top
        discretization_->u(i, discretization_->uJEnd() -1) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() -2);
    }
    // set v
    for ( int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
    {
        // bottom
        discretization_->u(i, discretization_->uJBegin() +1) = settings_.dirichletBcBottom[1];
        
        // top
        discretization_->u(i, discretization_->uJEnd() -1) = settings_.dirichletBcTop[1];
    }

    // sides - they take priority then
    // set u
    for ( int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
    {
        // left
        discretization_->u(discretization_->uIBegin() +1,j) = settings_.dirichletBcLeft[0];

        // top
        discretization_->u(discretization_->uIEnd() -1,j) = settings_.dirichletBcRight[0];
    }
    // set v
    for ( int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        // left
        discretization_->v(discretization_->vIBegin() +1,j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() +2, j);

        // top
        discretization_->v(discretization_->vIEnd() -1,j) = 2 * settings_.dirichletBcRight[1] - discretization_->u(discretization_->vIEnd() -2, j);
    }

    // set p
    pressureSolver_.setBoundaryValues()

    // set f,g boundary conditions the same as u,v
    // bottom and top
    // set f
    for ( int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
    {
        // bottom
        discretization_->f(i, discretization_->uJBegin() +1) = discretization_->u(i, discretization_->uJBegin() +1);
        
        // top
        discretization_->f(i, discretization_->uJEnd() -1) = discretization_->u(i, discretization_->uJEnd() -1);
    }
    // set g
    for ( int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
    {
        // bottom
        discretization_->g(i, discretization_->uJBegin() +1) = discretization_->u(i, discretization_->uJBegin() +1);
        
        // top
        discretization_->g(i, discretization_->uJEnd() -1) = discretization_->u(i, discretization_->uJEnd() -1);
    }

    // sides
    // calculating first the top/bottom then the sides gives the sides priority as required by the hint
    // set f
    for ( int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
    {
        // left
        discretization_->f(discretization_->uIBegin() +1,j) = discretization_->u(discretization_->uIBegin() +1,j);

        // top
        discretization_->f(discretization_->uIEnd() -1,j) = discretization_->u(discretization_->uIEnd() -1,j);
    }
    // set g
    for ( int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        // left
        discretization_->g(discretization_->vIBegin() +1,j) = discretization_->v(discretization_->vIBegin() +1,j);

        // top
        discretization_->g(discretization_->vIEnd() -1,j) = discretization_->v(discretization_->vIEnd() -1,j);
    }

}
    
void Computation::computePreliminaryVelocities()
{
    // calculate F
    for ( int j = discretization_->uJBegin() +1; j < discretization_->uJEnd() -1; j++)
    { 
        for ( int i = discretization_->uIBegin() +1; i < discretization_->uIEnd() -1; i++)
        {
            double diffusion = discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j);
            double convection = - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j);
            double sum = diffusion + convection + settings_.g[0];
            discretization_->f(i,j) = discretization_->u(i, j) + dt_ * (1 / settings_.re) * sum;
        }
    }

    // calculate G
    for ( int j = discretization_->vJBegin() +1; j < discretization_->vJEnd() -1; j++)
    { 
        for ( int i = discretization_->vIBegin() +1; i < discretization_->vIEnd() -1; i++)
        {
            double diffusion = discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j);
            double convection = - discretization_->computeDuvDx(i,j) - discretization_->computeDv2Dy(i,j);
            double sum = diffusion + convection + settings_.g[1];
            discretization_->g(i,j) = discretization_->v(i, j) + dt_ * (1 / settings_.re) * sum;
        }
    }
}
    
void Computation::computeRightHandSide()
{
    for ( int j = discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
    { 
        for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
        {
            double difference_f = (discretization_->f(i,j) - discretization_->f(i-1,j)) / meshWidth_[0];
            double difference_g = (discretization_->g(i,j) - discretization_->g(i,j-1)) / meshWidth_[1];
            discretization_->rhs(i,j) = (1 / dt_) * (difference_f + difference_g);
        }
    }
}
    
void Computation::computePressure()
{
    pressureSolver_.solve();
}

void Computation::computeVelocities()
{
    // calculate u
    for ( int j = discretization_->uJBegin() +1; j < discretization_->uJEnd() -1; j++)
    { 
        for ( int i = discretization_->uIBegin() +1; i < discretization_->uIEnd() -1; i++)
        {
            discretization_->u(i,j) = discretization_->f(i,j) - dt_ * discretization_->computeDpDx(i, j);
        }
    }

    // calculate v
    for ( int j = discretization_->vJBegin() +1; j < discretization_->vJEnd() -1; j++)
    { 
        for ( int i = discretization_->vIBegin() +1; i < discretization_->vIEnd() -1; i++)
        {
            discretization_->v(i,j) = discretization_->g(i,j) - dt_ * discretization_->computeDpDy(i, j);
        }
    }
}
