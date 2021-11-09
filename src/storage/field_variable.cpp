#pragma once

#include <memory>
#include "storage/field_variable.h"
#include "storage/array2D.h"

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
Array2D(size), 
origin_(origin), 
meshWidth_(meshWidth)
{
}


double FieldVariable::interpolateAt(double x, double y) const 
{    
    const double dx = meshWidth_[0]; // mesh width in x dir. 
    const double dy = meshWidth_[1]; // mesh width in y dir. 

    // indicies of (cell (i,j) in which the point (x,y) lies
    const int iLeftEdge  = (int) floor((x - origin_[0]) / dx); 
    const int jLowerEdge = (int) floor((y - origin_[1]) / dy); 

    // get values at corner points
    const double f_lowerLeft  = Array2D::operator()(iLeftEdge,     jLowerEdge);
    const double f_upperLeft  = Array2D::operator()(iLeftEdge,     jLowerEdge + 1); 
    const double f_lowerRight = Array2D::operator()(iLeftEdge + 1, jLowerEdge); 
    const double f_upperRight = Array2D::operator()(iLeftEdge + 1, jLowerEdge + 1);

    // relative position of x and y in the cell
    // one cell: |<-xr1-> x <-xr2->|
    //           |<--    dx      ->|

    const double xr1 = x % dx;   // relative position of x from left edge
    const double yr1 = y % dy;   // relative poistion of y from lower edge
    const double xr2 = dx - xr1; // distance right_edge - x
    const double yr2 = dy - yr1; // distance upper edge - y

    // bilinear interpolation
    f_intp = (f_lowerLeft * xr2 * yr2 
        + f_lowerRight * xr1 * yr2
        + f_upperLeft  * xr2 * yr1
        + f_upperRight * xr1 * xr1) / (dx * dy);

return f_intp;

}
