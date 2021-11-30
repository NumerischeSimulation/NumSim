#include "1_red_black.h"

RedBlack::RedBlack(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning):
  PressureSolver(discretization, epsilon, maximumNumberOfIterations),
  partitioning_(partitioning)
{
}

double RedBlack::calculateResidual()
{
    // cell size
    double dy = discretization_->dy();
    double dx = discretization_->dx();
    double dx2 = pow(dx,2);
    double dy2 = pow(dy,2);

   
    double res_local = 0.0;

    for ( int i = 0; i < discretization_->nCells()[0]; i++)
    { 
        for ( int j = 0; j < discretization_->nCells()[1]; j++)
        {
            // calculate residual 
            double pxx = (discretization_->p(i-1, j) - 2.0 *discretization_->p(i,j) + discretization_->p(i+1, j)) / (dx2);
            double pyy = (discretization_->p(i, j-1) - 2.0 *discretization_->p(i,j) + discretization_->p(i, j+1)) / (dy2);

            double resij = discretization_->rhs(i, j) - pxx - pyy;   
            res_local = res_local + (pow(resij,2));
        }
    }

    //calculate residual
    res_local = res_local/(discretization_->nCells()[0] * discretization_->nCells()[1]);
    double res_global;
    MPI_Allreduce(&res_local, &res_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
    return res_global;
}


void RedBlack::solve()
{
    // cell size
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
        for ( int j = 0; j < discretization_->nCells()[1]; j++)
        {
            // if is an even row, we start with the first column, else second
            if( j % 2 == 0)
            {
               
                for ( int i = 0; i < discretization_->nCells()[0]; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }

            }
            else
            {
                for ( int i = 1; i < discretization_->nCells()[0]; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }
            }
        }
        // apply the horizontal exchange first such that the boundary values are set correctly
        pExchangeHorizontal();
        pExchangeVertical();
        
        // Red Solver
        // one half solver iteration
        for ( int j = 0; j < discretization_->nCells()[1]; j++)
        {
            if( j % 2 == 0)
            {
                for ( int i = 1; i < discretization_->nCells()[0]; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }

            }
            else
            {
                for ( int i = 0; i < discretization_->nCells()[0]; i = i+2)
                {
                    double sum_x = (discretization_->p(i+1, j) + discretization_->p(i-1, j)) / (dx2);
                    double sum_y = (discretization_->p(i, j+1) + discretization_->p(i, j-1)) / (dy2);
                    discretization_->p(i, j) = factor * (sum_x + sum_y - discretization_->rhs(i, j));
                }
            }
        }

        // apply the horizontal exchange first such that the boundary values are set correctly
        pExchangeHorizontal();
        pExchangeVertical();
    
    
        iteration +=1;
        
        //set new boundary values
        setBoundaryValues();

        res = calculateResidual();
    }
}

void RedBlack::pExchangeHorizontal()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even processes: send left, receive left, send right, receive right
    if ((partitioning_->ownRankNo() % 2) == 0)
    {
        // left
        if (partitioning_->ownPartitionContainsLeftBoundary())
        {
            setBoundaryValuesLeft();
        } else 
        {
            // p
            std::array<double, discretization_->nCells()[1]+4> p_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                p_to_other_ghost_left[j+2] = discretization_->p(0, j); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_left;

            MPI_Send(&p_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownLeftNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(-1, j) = other_ghost_to_p_left[j+2]; // p_{-1}
            }
        }

        // right
        if (partitioning_->ownPartitionContainsRightBoundary())
        {
            setBoundaryValuesRight();
        } else
        {
            // p
            std::array<double, discretization_->nCells()[1]+4> p_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                p_to_other_ghost_right[j+2] = discretization_->p(discretization_->nCells()[1] -1, j); // p_{n-1}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_right;

            MPI_Send(&p_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownRightNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(discretization_->nCells()[1], j) = other_ghost_to_p_right[j+2]; // p_{n}
            }
        }
    }
    // the uneven processes: receive right, send right, receive left, send left
    if ((partitioning_->ownRankNo() % 2) == 1)
    {
        if (partitioning_->ownPartitionContainsRightBoundary())
        {
            setBoundaryValuesRight();
        } else
        {
            // p
            std::array<double, discretization_->nCells()[1]+4> p_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                p_to_other_ghost_right[j+2] = p(discretization_->nCells()[1] -1, j); // p_{n-1}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_right;

            MPI_Recv(&other_ghost_to_p_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownRightNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(discretization_->nCells()[1], j) = other_ghost_to_p_right[j+2]; // p_{n}
            }
        }
        if (partitioning_->ownPartitionContainsLeftBoundary())
        {
            setBoundaryValuesLeft();
        } else 
        {
            // p
            std::array<double, discretization_->nCells()[1]+4> p_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                p_to_other_ghost_left[j+2] = discretization_->p(0, j); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_left;

            MPI_Recv(&other_ghost_to_p_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownLeftNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_->ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(-1, j) = other_ghost_to_p_left[j+2]; // p_{-1}
            }
        }
    }
}

void RedBlack::pExchangeVertical()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even row processes: send top, receive top, send bottom, receive bottom
    if ((partitioning_->ownRankCoordinate()[1] % 2) == 0)
    {
        // top
        if (partitioning_->ownPartitionContainsTopBoundary())
        {
            setBoundaryValuesTop();
        } else 
        {
            // p
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_top[i+2] = discretization_->p(i, discretization_->nCells()[0] -1); // p_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_p_top;

            MPI_Send(&p_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownTopNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, discretization_->nCells()[0]) = other_ghost_to_p_top[i+2]; // p_{n}
            }
        }
        // bottom
        if (partitioning_->ownPartitionContainsBottomBoundary())
        {
            setBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_bottom[i+2] = discretization_->p(i, 0); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_bottom;

            MPI_Send(&p_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownBottomNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, -1) = other_ghost_to_p_bottom[j+2]; // p_{-1}
            }
        }
    }
    // the uneven row processes: receive bottom, send bottom, receive top, send top, 
    if ((partitioning_->ownRankCoordinate()[1] % 2) == 1)
    {
        // bottom
        if (partitioning_->ownPartitionContainsBottomBoundary())
        {
            setBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_bottom[i+2] = discretization_->p(i, 0); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_bottom;

            MPI_Recv(&other_ghost_to_p_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownBottomNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, -1) = other_ghost_to_p_bottom[j+2]; // p_{-1}
            }
        }
        // top
        if (partitioning_->ownPartitionContainsTopBoundary())
        {
            setBoundaryValuesTop();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_top[i+2] = discretization_->p(i, discretization_->nCells()[0] -1); // p_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_p_top;

            MPI_Recv(&other_ghost_to_p_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownTopNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_->ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, discretization_->nCells()[0]) = other_ghost_to_p_top[i+2]; // p_{n}
            }
    }
}