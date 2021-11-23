#pragma once

#include "0_computation.h"

#include <memory>
#include <iostream>  // for cout
#include <cmath>
#include <algorithm>

class ComputationParallel: public Computation
{
public:    
    //! initialize the computation object, parse the settings from file that is given as the only command line argument
    void initialize(int argc, char *argv[]);

    //! run the whole simulation until t_end 
    void runSimulation();

protected:
    //! compute the time step width dt from maximum velocities 
    void computeTimeStepWidth();

    //!  set boundary values of u and v to correct values
    void applyBoundaryValues();

    //! solve the Poisson equation for the pressure 
    void computePressure();

    //! compute the new velocities, u,v, from the preliminary velocities, F,G and the pressure, p
    void computeVelocities();

    std::unique_ptr<Partitioning> partitioning_;

    std::unique_ptr<OutputWriterParaviewParallel> OutputWriterParaviewParallel_;
    std::unique_ptr<OutputWriterTextParallel> OutputWriterTextParallel_;

};
