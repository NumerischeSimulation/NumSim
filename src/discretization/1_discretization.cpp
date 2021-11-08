
#include "0_staggered_grid.h"

// Contructor
Discretization::Discretization(std::array<int,2> nCells, std::array<double,2> meshWidth) :
    StaggeredGrid(nCells, meshWidth)
{
}


//! compute the 2nd derivative ∂^2 u / ∂x^2 at right edge
double computeD2uDx2(int i, int j) const
{
    const double hx = this->meshWidth_[0] ; // mesh width in x direction

    // return 2nd derivative ddu/dxx at right edge of cell i,j
    return 1./(hx*hx) * (u(i-1,j) - 2.*u(i,j) + u(i-1,j));

}

//! compute the 2nd derivative ∂^2 u / ∂y^2
double computeD2uDy2(int i, int j) const 
{
    const double hy = this->meshWidth_[1] ; // mesh width in y direction

    // return 2nd derivative ddu/dyy at right edge of cell i,j
    return 1./(hy*hy) * (u(i,j+1)-2.*u(i,j)+u(i,j-1));
}

//! compute the 2nd derivative ∂^2 v / ∂x^2
double computeD2vDx2(int i, int j) const 
{
    const double hx = this->meshWidth_[0] ; // mesh width in x direction

    // return 2nd derivative at top edge of cell
    return 1./(hx*hx) * (v(i+1,j)-2.*v(i,j)+v(i-1,j));
}

//! compute the 2nd derivative ∂^2 v / ∂y^2
double computeD2vDy2(int i, int j) const 
{
    const double hy = this->meshWidth_[1] ; // mesh width in y direction

    // return 2nd derivative at top edge of cell
    return 1./(hy*hy) * (v(i,j+1)-2.*v(i,j)+v(i,j-1));
}

//! compute the 1st derivative ∂p / ∂x
double computeDpDx(int i, int j) const 
{
    const double hx = this->meshWidth_[0] ; // mesh width in x direction

    // return 1st derivative at right edge of cell 
    return 1./hx = (p(i+1,j) - p(i,j));
}

//! compute the 1st derivative ∂p / ∂y
double computeDpDy(int i, int j) const 
{
    const double hy = this->meshWidth_[1] ; // mesh width in y direction

    // return 1st derivative at top edge of cell i,j
    return 1./hy = (p(i,j+1) - p(i,j));
}

