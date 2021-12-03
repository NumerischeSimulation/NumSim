#pragma once

#include "0_computation.h"
#include "discretization/1_discretization.h"
#include "pressure_solver/1_red_black.h"
#include "partitioning/partitioning.h"

#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"

#include <memory>
#include <iostream>  // for cout
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include <vector>

class ComputationParallel: public Computation
{
public:    
    //! initialize the computation object, parse the settings from file that is given as the only command line argument
    void initialize(int argc, char *argv[]); // overwrite

    //! run the whole simulation until t_end 
    void runSimulation(); // overwrite

    //! test communications of u/v as well as setting of boundary values
    void communicationTest();


protected:
    //! compute the time step width dt from maximum velocities 
    void computeTimeStepWidthParallel(double currentTime);

    //!  set boundary values of u and v to correct values
    void applyBoundaryValues();  // overwrite of method by Computation
    void applyBoundaryValuesLeft();
    void applyBoundaryValuesRight();
    void applyBoundaryValuesTop();
    void applyBoundaryValuesBottom();

    //! sets the ghost layers in the discretization correctly
    void uvExchangeVertical();
    void uvExchangeHorizontal();

    //! exchanges the ghost layers for u, v
    void exchange(int rankCorrespondent, //! rank of the other process with which to exchange 
                  int indexToSend, //! index of the slice to send or NULL
                  int indexFromReceive, //! index of the slice where to write the incoming data or NULL
                  char direction, //! x or y direction
                  char variable, //! u or v
                  bool ToFrom); //! send and then receive, else receive then send

    std::shared_ptr<Partitioning> partitioning_;

    std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaviewParallel_;
    std::unique_ptr<OutputWriterTextParallel> outputWriterTextParallel_;

};
