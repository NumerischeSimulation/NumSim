#pragma once

#include <memory>

#include "discretization/1_discretization.h"

class PressureSolver {
public:
    PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations );
    virtual void solve() = 0;

    void setBoundaryValues();

protected:
    // set boundary values

protected:
    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    int maximumNumberOfIterations_;
};