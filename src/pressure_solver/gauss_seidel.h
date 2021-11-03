#pragma once

#include "pressure_solver/pressure_solver.h"

class GaussSeidel: public PressureSolver{
public:
    void solve();
}