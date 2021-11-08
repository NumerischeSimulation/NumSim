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
        meshWidth_[i] = = settings_.physicalSize[i] / settings_.nCells[i];
    }

    // initialize
    if (settings_.useDonorCell)
    {
        std::cout << "DonorCell is not programmed yet!" << settings_.useDonorCell << std::endl;
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
    double dt_ = 0;
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    outputWriterText_ = std::make_unique<OutputWriter>(discretization_);
}

void Computation::runSimulation()
{
    double currentTime = 0;

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
    double boundary_diffusion = (settings_.re * settings_.physicalSize[0] * settings_.physicalSize[1]) / 4;

    // boundary from convection
    double u_max = std::max_element((discretization_->u_).begin(), (discretization_->u_).end());
    double v_max = std::max_element((discretization_->v_).begin(), (discretization_->v_).end());
    double boundary_convection_u = settings_.physicalSize[0] / u_max;
    double boundary_convection_v = settings_.physicalSize[1] / v_max;

    // together
    double min_dt = std::min({boundary_diffusion, boundary_convection_u, boundary_convection_v});

    // security factor
    dt_ =  min_dt * settings_.tau;
}
    
void Computation::applyBoundaryValues()
{
    return;
}
    
void Computation::computePreliminaryVelocities()
{
    return;
}
    
void Computation::computeRightHandSide()
{
    return;
}
    
void Computation::computePressure()
{
    return;
}

void Computation::computeVelocities()
{
    return;
}