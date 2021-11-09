#include "gauss_seidel.h"
#include <cmath>

void GaussSeidel::solve()
{ 
    int nCellsx;
    int nCellsy;
    nCellsx = discretization_->nCells()[0];
    nCellsy = discretization_->nCells()[1];



    //sell size
    const double dy = discretization_->dy();
    const double dx = discretization_->dx();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    int iteration = 0;
    double pxx;
    double pyy;
    double res;
    double resij;
    double ressum;
    int i;
    int j;

    // iterate through grid 
    
    while( iteration < maximumNumberOfIterations_ && res > dx2 * dy2 * epsilon_)
    {
        for ( int j = discretization_->pJBegin() +1 ; j < discretization_->pJEnd() -1; j++)
        { 
            for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
            {
            discretization_->p(i, j) = (dx2 * dy2) / ( 2 * dx2 * dy2) *(pxx + pyy - discretization_->rhs(i, j));
            pxx = (discretization_->p(i+1, j) - discretization_->p(i-1, j)) / (2 *dx2);
            pyy = (discretization_->p(i, j+1) - discretization_->p(i, j-1)) / (2 * dy2);

            // calculate residual 
            resij = discretization_->rhs(i, j) - pxx - pyy;   
            ressum = ressum + resij * resij;  

            }
           
        }
        
        iteration +=1;

        //calculate residual

        res = sqrt(dx2* dy2 * ressum )/(nCellsx * nCellsy);

        //set new boundary values
        setBoundaryValues();
    }
       

    return solve();
}

