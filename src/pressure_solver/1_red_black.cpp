#include "1_red_black.h"

RedBlack::RedBlack(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning):
  PressureSolver(discretization, epsilon, maximumNumberOfIterations),
  partitioning_(partitioning)
{
}

double RedBlack::calculateResidual()
{

    //std::cout << "Calculating Residual " << " (" << partitioning_->ownRankNo() << ")" << std::endl;

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
    double res_global;
    MPI_Allreduce(&res_local, &res_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    res_global = res_global/(partitioning_->nCellsGlobal()[0] *partitioning_->nCellsGlobal()[1]);

    // std::cout << "The global residual is: " << res_global << " and the local " << res_local << " (" << partitioning_->ownRankNo() << ")" << std::endl;
  
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

    //std::cout << "Starting with RedBlack iterations... " << std::endl;
    
    // iterate through grid 
    while( iteration < maximumNumberOfIterations_ && res > pow(epsilon_,2))
    {
        // std::cout << "Iter: " << iteration << std::endl;
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
        //std::cout << "Starting exchange... " << std::endl;
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

        res = calculateResidual();

        
    }

    //std::cout << "RedBlack: " << iteration << " with a residuum of " << res << " from target " << std::pow(epsilon_,2) << std::endl;
}




// for testing of boundary, communications and red black pattern

/*
void RedBlack::solve()  // communicationTest
{
        // Black Solver
        // one half solver iteration
        for (int j = 0; j < discretization_->nCells()[1]; j++)
        {
            // if is an even row, we start with the first column, else second
            if (j % 2 == 0)
            {
                for (int i = 0; i < discretization_->nCells()[0]; i = i + 2)
                {
                    discretization_->p(i, j) = partitioning_->ownRankNo()+1;
                }
            }
            else
            {
                for (int i = 1; i < discretization_->nCells()[0]; i = i + 2)
                {
                    discretization_->p(i, j) = partitioning_->ownRankNo()+1;
                }
            }
        }
        
        // apply the horizontal exchange first such that the boundary values are set correctly
        // std::cout << "Starting exchange... " << std::endl;
        pExchangeHorizontal();
        pExchangeVertical();
    
} 
*/

void RedBlack::pExchangeHorizontal()
{
    //std::cout << "pExchangeHorizontal" << std::endl;

    // the even processes: send left, receive left, send right, receive right
    if ((partitioning_->ownRankCoordinate()[0] % 2) == 0)
    {
        // left
        if (partitioning_->ownPartitionContainsLeftBoundary())
        {
            setBoundaryValuesLeft();
        } else 
        {
            // p
            // send p_0, receive p_-1
            exchange(partitioning_->ownLeftNeighbour(), 
                     0, -1, 
                     'x', true);
        }

        // right
        if (partitioning_->ownPartitionContainsRightBoundary())
        {
            setBoundaryValuesRight();
        } else
        {
            // p
            // send p_n-1, receive p_n
            exchange(partitioning_->ownRightNeighbour(), 
                     discretization_->nCells()[0] -1, discretization_->nCells()[0],
                     'x', true);
        }
    }
    // the uneven processes: receive right, send right, receive left, send left
    if ((partitioning_->ownRankCoordinate()[0] % 2) == 1)
    {
        if (partitioning_->ownPartitionContainsRightBoundary())
        {
            setBoundaryValuesRight();
        } else
        {
            // p
            // receive p_n, send p_n-1
            exchange(partitioning_->ownRightNeighbour(), 
                     discretization_->nCells()[0] -1, discretization_->nCells()[0],
                     'x', false);
        }
        if (partitioning_->ownPartitionContainsLeftBoundary())
        {
            setBoundaryValuesLeft();
        } else 
        {
            // p
            // receive p_-1, send p_0
            exchange(partitioning_->ownLeftNeighbour(), 
                     0, -1,
                     'x', false);
        }
    }
}

void RedBlack::pExchangeVertical()
{
    //std::cout << "pExchangeVertical" << std::endl;  

    // the even row processes: send top, receive top, send bottom, receive bottom
    if ((partitioning_->ownRankCoordinate()[1] % 2) == 1)
    {
        // top
        if (partitioning_->ownPartitionContainsTopBoundary())
        {
            setBoundaryValuesTop();
        } else 
        {
            // p
            // send p_n-1, receive p_n
            exchange(partitioning_->ownTopNeighbour(), 
                     discretization_->nCells()[1] -1, discretization_->nCells()[1], 
                     'y', true);
        }
        // bottom
        if (partitioning_->ownPartitionContainsBottomBoundary())
        {
            setBoundaryValuesBottom();
        } else
        {
            // p
            // send p_0, receive p_-1
            exchange(partitioning_->ownBottomNeighbour(), 
                     0, -1, 
                     'y', true);
        }
    }
    // the uneven row processes: receive bottom, send bottom, receive top, send top, 
    if ((partitioning_->ownRankCoordinate()[1] % 2) == 0)
    {
        // bottom
        if (partitioning_->ownPartitionContainsBottomBoundary())
        {
            setBoundaryValuesBottom();
        } else
        {
            // p
            // receive p_-1, send p_0
            exchange(partitioning_->ownBottomNeighbour(), 
                     0, -1, 
                     'y', false);
        }
        // top
        if (partitioning_->ownPartitionContainsTopBoundary())
        {
            setBoundaryValuesTop();
        } else 
        {
            // p
            // receive p_n, send p_n-1
            exchange(partitioning_->ownTopNeighbour(), 
                     discretization_->nCells()[1]-1, discretization_->nCells()[1], 
                     'y', false);
        }
    }
}

void RedBlack::exchange(int rankCorrespondent, int indexToSend, int indexFromReceive, char direction, bool ToFrom)
{
    int ownRankNo = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
    
    // std::cout << "Exchanges " << ownRankNo << " to " <<  rankCorrespondent << " with data slices " << indexToSend << " to " << indexFromReceive << " " << std::endl
    //          << " in " << direction << " with " << ToFrom << " (" << partitioning_->ownRankNo() << ")" << std::endl;
    // index to or from can be mpi null pointer

    // initialize variables
    int nValuesCommunication = 0;
    int offset = 0;

    // get constant variables for each case
    if (direction == 'x')
    {
        nValuesCommunication = discretization_->pJEnd() - discretization_->pJBegin();
        offset = std::abs(discretization_->pJBegin());
    } else if (direction == 'y')
    {
        nValuesCommunication = discretization_->pIEnd() - discretization_->pIBegin();
        offset = std::abs(discretization_->pIBegin());
    }

    // get slice to communicate for each case
    std::vector<double> toOtherGhost(nValuesCommunication, 5);
    if (direction == 'x')
    {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
        {
            toOtherGhost[j+offset] = discretization_->p(indexToSend, j);
        }
    } else if (direction == 'y')
    {
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            toOtherGhost[i+offset] = discretization_->p(i, indexToSend);
        }
    }

    // send-receive or receive-send depending on the direction
    std::vector<double> otherGhostFrom(nValuesCommunication, 7);
    if (ToFrom)
    {
        // send
        MPI_Send(toOtherGhost.data(),
            nValuesCommunication,
            MPI_DOUBLE,
            rankCorrespondent,
            0,
            MPI_COMM_WORLD
            );

        // receive
        MPI_Recv(otherGhostFrom.data(), 
            nValuesCommunication,
            MPI_DOUBLE,
            rankCorrespondent,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
            );
    } else 
    {
        // receive
        MPI_Recv(otherGhostFrom.data(), 
            nValuesCommunication,
            MPI_DOUBLE,
            rankCorrespondent,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
            );

        // send
        MPI_Send(toOtherGhost.data(),
            nValuesCommunication,
            MPI_DOUBLE,
            rankCorrespondent,
            0,
            MPI_COMM_WORLD
            );
    }

    // write slice to correct index
    if (direction == 'x')
    {
        for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++)
        {
            discretization_->p(indexFromReceive, j) = otherGhostFrom[j+offset];
        }
    } else if (direction == 'y')
    {
        for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++)
        {
            discretization_->p(i, indexFromReceive) = otherGhostFrom[i+offset];
        }
    }
}
