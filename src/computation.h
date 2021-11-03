#pragma once

#include "settings.h"
#include "discretization/1_discretization.h"
#include "pressure_solver/pressure_solver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"

#include <memory>

class Computation{
public:
    void initialize(int argc, char *argv[]);
    void runSimulation();
private:
    void computeTimeStepWidth();
    void applyBoundaryValues();
    void computePreliminaryVelocities();
    void computeRightHandSide();
    void computePressure();
    void computeVelocities();

    Settings settings_;
    std::shared_ptr<Discretization> distcretization_;
    std::unique_ptr<PressureSolver> pressureSolver_;
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;
    std::array<double, 2> meshWidth_;
    double dt_;

}