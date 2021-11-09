#pragma once

#include "sor.h"
#include "pressure_solver.h"

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega)
{
discretization_ = discretization;
epsilon = epsilon;
maximumNumberOfIterations_ = maximumNumberOfIterations;
omega_ = omega;
}