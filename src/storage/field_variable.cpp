#include "storage/field_variable.h"

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
Array2D(size), 
origin_(origin), 
meshWidth_(meshWidth)
{};


double FieldVariable::interpolateAt(double x, double y) const 
{    
    const double dx = meshWidth_[0]; // mesh width in x dir. 
    const double dy = meshWidth_[1]; // mesh width in y dir. 

    // indicies of (cell (i,j) in which the point (x,y) lies
    int iLeftEdge  = (int) std::floor((x - origin_[0]) / dx); 
    int jLowerEdge = (int) std::floor((y - origin_[1]) / dy); 

    // shift right and upper boundaries so that they don't use cells outside of the grid
    if (iLeftEdge <= 0)
    {
        iLeftEdge = 0; // shift it one column to the right
    }
    if (jLowerEdge <= 0) 
    {
        jLowerEdge = 0; // shift it up
    }
    // relative position of x and y in the cell
    // one cell: |<-xr1-> x <-xr2->|
    //           |<--    dx      ->|
    const double xr1 = x  - (meshWidth_[0]*iLeftEdge + origin_[0]);   // relative position of x from left edge
    const double yr1 = y  - (meshWidth_[1]*jLowerEdge + origin_[1]);   // relative poistion of y from lower edge
    const double xr2 = dx - xr1; // distance right_edge - x
    const double yr2 = dy - yr1; // distance upper edge - y

    // get values at corner points
    assert(0 <= iLeftEdge && 0 <= jLowerEdge);
    const double f_lowerLeft  = Array2D::operator()(iLeftEdge,     jLowerEdge);
    const double f_upperLeft  = Array2D::operator()(iLeftEdge,     jLowerEdge + 1); 
    const double f_lowerRight = Array2D::operator()(iLeftEdge + 1, jLowerEdge); 
    const double f_upperRight = Array2D::operator()(iLeftEdge + 1, jLowerEdge + 1);

    // bilinear interpolation
    const double f_intp = (f_lowerLeft * xr2 * yr2 
                        + f_lowerRight * xr1 * yr2
                        + f_upperLeft  * xr2 * yr1
                        + f_upperRight * xr1 * xr1) / (dx * dy);

return f_intp;

}
