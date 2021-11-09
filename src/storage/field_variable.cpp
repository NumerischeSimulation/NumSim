#pragma once

#include <memory>
#include "storage/array2D.h"

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
size_=size,
origin_=origin,
meshWidth_=meshWidth {

}


double FieldVariable::interpolateAt(double x, double y) const 
{    
const double hx = meshWidth_[0]; // mesh width in x dir. 
const double hy = meshWidth_[1]; // mesh width in y dir. 

// indicies of (cell (i,j) in which the point (x,y) lies
const double iLeftEdge  = floor(x / hx); 
const double iRightEdge = ceil(x / hx); 

const double jLowerEdge = floor(y / hy); 
const double jUpperEdge = ceil(y / hy); 

// edge points of cell i,j
const double[2] edgeLowerLeft  = {iLeftEdge, jLowerEdge};
const double[2] edgeUpperLeft  = {iLeftEdge, jUpperEdge};

const double[2] edgeLowerRight = {iRightEdge, jLowerEdge};
const double[2] edgeUpperRight = {iRightEdge, jUpperEdge};

// values at the points
const double edgeLowerLeft_value = data_[(int)((edgeLowerLeft[0]-1) * (size_[0])) + edgeLowerLeft[1])];
const double edgeUpperLeft_value = data_[(int)((edgeUpperLeft[0]-1) * (size_[0])) + edgeUpperLeft[1])];
const double edgeLowerRight_value = data_[(int)((edgeLowerRight[0]-1) * (size_[0])) + edgeLowerRight[1])];
const double edgeUpperRight_value = data_[(int)((edgeUpperRight[0]-1) * (size_[0])) + edgeUpperRight[1])];

// relative position of x and y in the cell
const double dx = x % hx; // relative position of x from left edge
const double dy = x % hy; // relative poistion of y from lower edge

// bilinear interpolation
const double i_intp_LowerEdge = (iLeftEdge - x )/(iRightEdge - iLeftEdge) * edgeLowerLeft_value + (x - iLeftEdge)/(iRightEdge - iLeftEdge) * edgeLowerRight_value;
const double i_intp_UpperEdge = (iLeftEdge - x )/(iRightEdge - iLeftEdge) * edgeUpperLeft_value + (x - iLeftEdge)/(iRightEdge - iLeftEdge) * edgeUpperRight_value;

const double j_intp_LeftEdge  = 

return ;

}
