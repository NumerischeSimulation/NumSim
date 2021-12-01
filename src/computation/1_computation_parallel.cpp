#include "computation/1_computation_parallel.h"

void ComputationParallel::initialize(int argc, char *argv[])
{
    // parse the parameters
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();

    // start MPI
    MPI_Init(&argc, &argv);

    // get own rank number and total number of ranks
    int ownRankNo = 0;
    int nRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    // get own partion depending on where in the domain the process lies
    partitioning_ = std::make_shared<Partitioning>(ownRankNo, nRanks, settings_.nCells);


    // calculate
    for (int i = 0; i < 2; i++)
    {
        meshWidth_[i] =  settings_.physicalSize[i] / settings_.nCells[i];
        std::cout << "computed mesh width " << i << ": " << meshWidth_[i] << " " << settings_.physicalSize[i] << " " << settings_.nCells[i] << std::endl;
    }

    // initialize
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(partitioning_->nCellsLocal(), meshWidth_, partitioning_->ownPartitionNeighbours(), settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(partitioning_->nCellsLocal(), meshWidth_, partitioning_->ownPartitionNeighbours());
    }
    
    pressureSolver_ = std::make_unique<RedBlack>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_); 

    // misc
    dt_ = 0.;

    // initialize output writer
    outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
    outputWriterTextParallel_     = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);
}

void ComputationParallel::runSimulation()
{
    double currentTime = 0;

    std::cout << "+++++++++++++++++++++++" << std::endl;
    std::cout << "Starting at time: " << currentTime << std::endl;
    std::cout << "+++++++++++++++++++++++" << std::endl;

    // the steps correspond to the steps in our algorithm in the overleaf or docs/numsim-algos.tex
    while (currentTime < settings_.endTime)
    {
        // step 1: set the boundary values / exchange the final velocities at the borders
        applyBoundaryValues();
        std::cout << "Applied boundary values for u/v and F/G." << std::endl;

        // step 2: compute time step width
        computeTimeStepWidthParallel(currentTime);
        
        currentTime += dt_;
    
        std::cout << "+++++++++++++++++++++++" << std::endl;
        std::cout << "current Time: " << currentTime << std::endl;
        std::cout << "+++++++++++++++++++++++" << std::endl;    

        // step 4: calculate F, G with first setting the boundary conditions of F, G (step 3)
        computePreliminaryVelocities();

        std::cout << "Computed preliminary velocities" << std::endl;

        // step 5: compute the right hand side of the pressure equation
        computeRightHandSide();

        std::cout << "Computed right hand side" << std::endl;

        // step 6: solve the pressure equation
        computePressure();

        std::cout << "Computed presure" << " from process number " << partitioning_->ownRankNo() << std::endl;

        // step 7: calculate the final velocities
        computeVelocities();

        std::cout << "Computed velocities" << std::endl;

        // step 9: write output
        if (std::floor(currentTime) == currentTime)
        {
            outputWriterParaview_->writeFile(currentTime);
            outputWriterText_->writeFile(currentTime);
        }
    }
    
    // end the MPI-session
    MPI_Finalize();
}

void ComputationParallel::computeTimeStepWidthParallel(double currentTime)
{
    // compute timestep for each subdomain
    computeTimeStepWidth();

    // use the minimum as global timestep
    double dt_global;
    double dt_local = dt_;
    MPI_Allreduce(&dt_local, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    // if necessary adapt so that every full second is reached
    if (std::floor(currentTime + dt_) == std::floor(currentTime) + 1)
    {
        dt_global = (double) (std::floor(currentTime) + 1) - currentTime; // currentTime hits exactly next second
    }

    // or if necessary set currentTime + dt_ to endTime 
    if (currentTime + dt_global > settings_.endTime)
    {
        dt_global = settings_.endTime - currentTime;

        std::cout << std::endl;
        std::cout << "Final time step!" << std::endl;
    }

    dt_ = dt_global;
}
    
void ComputationParallel::applyBoundaryValues()
{
    // first exchange horizontal, then vertical for the corner cells
    uvExchangeHorizontal();

    uvExchangeVertical();
}

void ComputationParallel::applyBoundaryValuesLeft()
{
    // u
    for ( int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
    {
        discretization_->u(discretization_->uIBegin(),j)   = settings_.dirichletBcLeft[0];
        }
    // v
    for ( int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        discretization_->v(discretization_->vIBegin() ,j)  = 2. * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin() + 1, j);
    }
    // f
    for ( int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
    {
        discretization_->f(discretization_->uIBegin(),j) = discretization_->u(discretization_->uIBegin(),j);
    }
    // g
    for ( int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        discretization_->g(discretization_->vIBegin(),j) = discretization_->v(discretization_->vIBegin(),j);
    }
}

void ComputationParallel::applyBoundaryValuesRight()
{
    // u
    for ( int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
    {
        discretization_->u(discretization_->uIEnd() -1 ,j) = settings_.dirichletBcRight[0];
    }
    // v
    for ( int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        discretization_->v(discretization_->vIEnd() - 1,j) = 2. * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd()  - 2, j);
    }
    // f
    for ( int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
    {
        discretization_->f(discretization_->uIEnd() -1,j) = discretization_->u(discretization_->uIEnd() -1,j);
    }
    // g
    for ( int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    {
        discretization_->g(discretization_->vIEnd() -1,j) = discretization_->v(discretization_->vIEnd() -1,j);
    }
}

void ComputationParallel::applyBoundaryValuesBottom()
{
    // bottom, set boundaries only in domain as corners belong to sides, computational domain begins at idx 0
    
    for ( int i = 0; i < discretization_->nCells()[0]; i++)
    {
        // u 
        discretization_->u(i, -1)  = 2. * settings_.dirichletBcBottom[0] - discretization_->u(i, 0);
        // v
        discretization_->v(i, -1)  = settings_.dirichletBcBottom[1];
        // f
        discretization_->f(i, -1)  = discretization_->u(i, 0);
        // g
        discretization_->g(i, -1)  = discretization_->v(i, 0);
    }
}

void ComputationParallel::applyBoundaryValuesTop()
{
    // set boundaries only in domain as corners belong to size, domain begins at idx 0
    
    for ( int i = 0; i < discretization_->nCells()[0]; i++)
    {
        // u
        discretization_->u(i, discretization_->nCells()[1]) = 2. * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->nCells()[1] - 1);
        // v
        discretization_->v(i, discretization_->nCells()[1]) = settings_.dirichletBcTop[1];
        // f
        discretization_->f(i, discretization_->nCells()[1]) = discretization_->u(i, discretization_->nCells()[1] - 1);
        // g
        discretization_->g(i, discretization_->nCells()[1]) = discretization_->v(i, discretization_->nCells()[1] - 1);
    }
}

void ComputationParallel::uvExchangeHorizontal()
{
    // the even processes: send left, receive left (u inner, u outer), send right (u inner, u outer), receive right
    if ((partitioning_->ownRankNo() % 2) == 0)
    {
        // left
        if (partitioning_->ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // u
            // send u_0, receive u_-2
            exchange(partitioning_->ownLeftNeighbour(),
                     0, -2, 
                     'x', 'u', true);

            // v
            // send v_0, receive v_-1
            exchange(partitioning_->ownLeftNeighbour(),
                     0, -1, 
                     'x', 'v', true);
        }

        // right
        if (partitioning_->ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // u
            // send u_n-2, receive u_n
            exchange(partitioning_->ownRightNeighbour(),
                     discretization_->nCells()[0] -2, discretization_->nCells()[0], 
                     'x', 'u', true);

            // v
            // send v_n-1, receive v_n
            exchange(partitioning_->ownRightNeighbour(),
                     discretization_->nCells()[0] -1, discretization_->nCells()[0], 
                     'x', 'v', true);
        }
    }
    // the uneven processes: receive right, send right, receive left (u inner, u outer), send left
    if ((partitioning_->ownRankNo() % 2) == 1)
    {
        if (partitioning_->ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // u
            // receive u_n-2, send u_n
            exchange(partitioning_->ownRightNeighbour(),
                     discretization_->nCells()[0], discretization_->nCells()[0] -2, 
                     'x', 'u', false);

            // v
            // receive v_n-1, send v_n
            exchange(partitioning_->ownRightNeighbour(),
                     discretization_->nCells()[0], discretization_->nCells()[0] -1,
                     'x', 'v', false);
        }
        if (partitioning_->ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // u
            // receive u_0, send u_-2
            exchange(partitioning_->ownLeftNeighbour(),
                     -2, 0, 
                     'x', 'u', false);

            // v
            // receive v_0, send_-1
            exchange(partitioning_->ownLeftNeighbour(),
                     -1, 0,
                     'x', 'v', false);

        }
    }
}

void ComputationParallel::uvExchangeVertical()
{

    // the even processes: send top, receive top, send bottom (v inner, v outer), receive bottom
    if ((partitioning_->ownRankNo() % 2) == 0)
    {
        // top
        if (partitioning_->ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // u
            // send u_n-1, receive u_n
            exchange(partitioning_->ownTopNeighbour(),
                     discretization_->nCells()[1] -1, discretization_->nCells()[1], 
                     'y', 'u', true);

            // v
            // send v_n-2, receive v_n
            exchange(partitioning_->ownTopNeighbour(),
                     discretization_->nCells()[1] -2, discretization_->nCells()[1], 
                     'y', 'v', true);

        }
        // bottom
        if (partitioning_->ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            // send u_0, receive u_-1
            exchange(partitioning_->ownBottomNeighbour(),
                     0, -1, 
                     'y', 'u', true);

            // v
            // send v_0, receive v_-2
            exchange(partitioning_->ownBottomNeighbour(),
                     0, -2, 
                     'y', 'v', true);
        }
    }
    // the uneven processes: receive bottom, send bottom, receive top (v inner, v outer), send top, 
    if ((partitioning_->ownRankNo() % 2) == 1)
    {
        // bottom
        if (partitioning_->ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            // receive u_0, send u_-1
            exchange(partitioning_->ownBottomNeighbour(),
                     -1, 0, 
                     'y', 'u', false);

            // v
            // receive v_0, send v_-2
            exchange(partitioning_->ownBottomNeighbour(),
                     -2, 0, 
                     'y', 'v', false);
        }
        // top
        if (partitioning_->ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // u
            // receive u_n-1, send u_n
            exchange(partitioning_->ownTopNeighbour(),
                     discretization_->nCells()[1], discretization_->nCells()[1] -1, 
                     'y', 'u', false);

            // v
            // receive v_n-2, send v_n
            exchange(partitioning_->ownTopNeighbour(),
                     discretization_->nCells()[1], discretization_->nCells()[1] -2, 
                     'y', 'v', false);
        }
    }
}

void ComputationParallel::exchange(int rankCorrespondent, int indexToSend, int indexFromReceive, char direction, char variable, bool ToFrom)
{
    // index to or from can be NULL

    // initialize variables
    int nValuesCommunication = 0;
    int offset = 0;

    // get constant variables for each case
    if (direction == 'x')
    {
        if (variable == 'u')
        {
            nValuesCommunication = discretization_->uJEnd() - discretization_->uJBegin();
            offset = std::abs(discretization_->uJBegin());
        } else if (variable == 'v')
        {
            nValuesCommunication = discretization_->vJEnd() - discretization_->vJBegin();
            offset = std::abs(discretization_->vJBegin());
        } else 
        {
            throw;
        }
    } else if (direction == 'y')
    {
        if (variable == 'u')
        {
            nValuesCommunication = discretization_->uIEnd() - discretization_->uIBegin();
            offset = std::abs(discretization_->uIBegin());
        } else if (variable == 'v')
        {
            nValuesCommunication = discretization_->vIEnd() - discretization_->vIBegin();
            offset = std::abs(discretization_->vIBegin());
        } else 
        {
            throw;
        }
    }

    // get slice to communicate for each case
    std::vector<double> toOtherGhost(nValuesCommunication, 0);
    if (direction == 'x')
    {
        if (variable == 'u')
        {
            if (indexToSend != NULL){
                for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
                {
                    toOtherGhost[j+offset] = discretization_->u(indexToSend, j);
                }
            }
        } else if (variable == 'v')
        {
            if (indexToSend != NULL){
                for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
                {
                    toOtherGhost[j+offset] = discretization_->v(indexToSend, j);
                }
            }
        } else 
        {
            throw;
        }
    } else if (direction == 'y')
    {
        if (variable == 'u')
        {
            if (indexToSend != NULL){
                for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
                {
                    toOtherGhost[i+offset] = discretization_->u(i, indexToSend);
                }
            }
        } else if (variable == 'v')
        {
            if (indexToSend != NULL){
                for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
                {
                    toOtherGhost[i+offset] = discretization_->v(i, indexToSend);
                }
            }
        } else 
        {
            throw;
        }
    }

    // send-receive or receive-send depending on the direction
    std::vector<double> otherGhostFrom(nValuesCommunication, 0);
    if (ToFrom)
    {
        // send
        if (indexToSend != NULL)
        {
            MPI_Send(toOtherGhost.data(),
                nValuesCommunication,
                MPI_DOUBLE,
                rankCorrespondent,
                0,
                MPI_COMM_WORLD
                );
        }

        // receive
        if (indexFromReceive != NULL)
        {
            MPI_Recv(otherGhostFrom.data(), 
                nValuesCommunication,
                MPI_DOUBLE,
                rankCorrespondent,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE
                );
        }
    } else 
    {
        // receive
        if (indexFromReceive != NULL)
        {
            MPI_Recv(otherGhostFrom.data(), 
                nValuesCommunication,
                MPI_DOUBLE,
                rankCorrespondent,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE
                );
        }

        // send
        if (indexToSend != NULL)
        {
            MPI_Send(toOtherGhost.data(),
                nValuesCommunication,
                MPI_DOUBLE,
                rankCorrespondent,
                0,
                MPI_COMM_WORLD
                );
        }
    }

    // write slice to correct index
    if (direction == 'x')
    {
        if (variable == 'u')
        {
            if (indexToSend != NULL){
                for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
                {
                    discretization_->u(indexFromReceive, j) = otherGhostFrom[j+offset];
                }
            }
        } else if (variable == 'v')
        {
            if (indexToSend != NULL){
                for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
                {
                    discretization_->v(indexFromReceive, j) = otherGhostFrom[j+offset];
                }
            }
        } else 
        {
            throw;
        }
    } else if (direction == 'y')
    {
        if (variable == 'u')
        {
            if (indexToSend != NULL){
                for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
                {
                    discretization_->u(i, indexFromReceive) = otherGhostFrom[i+offset];
                }
            }
        } else if (variable == 'v')
        {
            if (indexToSend != NULL){
                for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
                {
                    discretization_->v(i, indexFromReceive) = otherGhostFrom[i+offset];
                }
            }
        } else 
        {
            throw;
        }
    }
}
