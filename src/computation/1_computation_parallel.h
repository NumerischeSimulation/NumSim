#pragma once

#include "0_computation.h"
#include "discretization/1_discretization.h"
#include "data_transfer/0_data_transfer.h"

#include <memory>
#include <iostream>  // for cout
#include <cmath>
#include <algorithm>
#include <mpi.h>

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
    void applyBoundaryValuesLeft();
    void applyBoundaryValuesRight();
    void applyBoundaryValuesTop();
    void applyBoundaryValuesBottom();

    std::unique_ptr<Partitioning> partitioning_;

    std::unique_ptr<OutputWriterParaviewParallel> OutputWriterParaviewParallel_;
    std::unique_ptr<OutputWriterTextParallel> OutputWriterTextParallel_;

};
