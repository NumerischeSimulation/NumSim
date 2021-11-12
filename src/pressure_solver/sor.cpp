#include "sor.h"

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations),
    omega_(omega)
{
}

void SOR::solve()
{ 
    int nCellsx = discretization_->nCells()[0] -2; // inner cells
    int nCellsy = discretization_->nCells()[1] -2; // inner cells

    //sell size
    const double dy = discretization_->dy();
    const double dx = discretization_->dx();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    int iteration = 0;
    double res = 10000000.; // random large value

    //set new boundary values
    setBoundaryValues();

    // iterate through grid 
     while( iteration < maximumNumberOfIterations_ && res > std::pow(epsilon_,2))
    {
        
        // one solver iteration
        for ( int j = discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
        { 
            for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
            {
            double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
            double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);    
            discretization_->p(i, j) = (1 - omega_) * discretization_->p(i, j) + (omega_ * dx2 * dy2) / ( 2 * (dx2 + dy2)) *(sum_x + sum_y - discretization_->rhs(i, j));
            }
        }
        
        // compute residuum after complete iteration of solver
        double ressum = 0;
        
        for ( int j = discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
        { 
            for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
            {
                // calculate residual 
                double pxx = (discretization_->p(i+1, j) - 2*discretization_->p(i,j) + discretization_->p(i-1, j)) / (dx2);
                double pyy = (discretization_->p(i, j+1) - 2*discretization_->p(i,j) + discretization_->p(i, j-1)) / (dy2);

                double resij = discretization_->rhs(i, j) - pxx - pyy;   
                ressum = ressum + resij * resij;  
            }
        }
        
        //calculate residual
        res = ressum/(nCellsx * nCellsy);

        iteration +=1;
        
        //set new boundary values
        setBoundaryValues();
    }
    
    std::cout << "SOR: " << iteration << " with a residuum of " << res << " from target " << std::pow(epsilon_,2) << std::endl;
       
}
