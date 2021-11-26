#include "1_red_black.h"

double RedBlack::calculateResidual()
  {
    int nCellsx = discretization_->nCells()[0] -4; // inner cells
    int nCellsy = discretization_->nCells()[1] -4; // inner cells

    //sell size
    double dy = discretization_->dy();
    double dx = discretization_->dx();
    double dx2 = pow(dx,2);
    double dy2 = pow(dy,2);

   
    double res_local = 0.0;

        for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i++)
        { 
            for ( int j= discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
            {
                // calculate residual 
                double pxx = (discretization_->p(i-1, j) - 2.0 *discretization_->p(i,j) + discretization_->p(i+1, j)) / (dx2);
                double pyy = (discretization_->p(i, j-1) - 2.0 *discretization_->p(i,j) + discretization_->p(i, j+1)) / (dy2);

                double resij = discretization_->rhs(i, j) - pxx - pyy;   
                res_local = res_local + (pow(resij,2));
            }
        }
        
        //calculate residual
         res_local = res_local/(nCellsx * nCellsy);
         double res_global;
          MPI_Allreduce(&res_local, &res_global, 1, MPI_DOUBLE, MPI_SUM,
              MPI_COMM_WORLD);
         return res_global;

  }


void RedBlack::solve
{
    int nCellsx = discretization_->nCells()[0] -4;
    int nCellsy = discretization_->nCells()[1] -4;



    //sell size
    double dy = discretization_->dy();
    double dx = discretization_->dx();
    double dx2 = pow(dx,2);
    double dy2 = pow(dy,2);
    double factor = (dx2 * dy2) / ( 2.0 * (dx2 + dy2));

    int iteration = 0;
    
    //initial residual
    double res = calculateResidual();

    // iterate through grid 
    
    
    while( iteration < maximumNumberOfIterations_ && res > pow(epsilon_,2))
    {
        // Black Solver
        // one half solver iteration
        for ( int j = discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
        {
            if( j % 2 == 0)
            {
                for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }

            }
            else
            {
                for ( int i = discretization_->pIBegin() +2; i < discretization_->pIEnd() -1; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }
            }
        }
        //TODO signature 
        DataTransfer::pExchangeHorizontal();
        DataTransfer::pExchangeVertical();
        
        // Red Solver
        // one half solver iteration
        for ( int j = discretization_->pJBegin() +1; j < discretization_->pJEnd() -1; j++)
        {
            if( j % 2 == 0)
            {
                for ( int i = discretization_->pIBegin() +2; i < discretization_->pIEnd() -1; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }

            }
            else
            {
                for ( int i = discretization_->pIBegin() +1; i < discretization_->pIEnd() -1; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }
            }
        }

        //TODO signature 
        DataTransfer::pExchangeHorizontal();
        DataTransfer::pExchangeVertical();
    
    
        iteration +=1;
        
        //set new boundary values
        setBoundaryValues();

        res = calculateResidual();
    }




}