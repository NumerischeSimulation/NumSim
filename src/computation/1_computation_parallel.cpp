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
    partitioning_ = Partitioning(ownRankNo, nRanks, settings_.nCells)


    // calculate
    for (int i = 0; i < 2; i++)
    {
        meshWidth_[i] =  settings_.physicalSize[i] / settings_.nCells[i];
        std::cout << "computed mesh width " << i << ": " << meshWidth_[i] << " " << settings_.physicalSize[i] << " " << settings_.nCells[i] << std::endl;
    }

    // initialize
    if (settings_.useDonorCell)
    {
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);
    }
    
    pressureSolver_ = std::make_unique<RedBlack>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_); 

    // misc
    dt_ = 0.;

    // initialize output writer
    outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_);
    outputWriterTextParallel_     = std::make_unique<OutputWriterTextParallel>(discretization_);
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

        std::cout << "Computed presure" << " from process number " << partitioning_.ownRankNo() << std::endl;

        // step 7: calculate the final velocities
        computeVelocities();

        std::cout << "Computed velocities" << std::endl;

        // step 9: write output
        if (std:floor(currentTime) == currentTime)
        {
            outputWriterParaview_->writeFile(currentTime);
            outputWriterText_->writeFile(currentTime);
        }
    }
    
    // end the MPI-session
    MPI_Finalize()
}

void ComputationParallel::computeTimeStepWidthParallel(double currentTime)
{
    // compute timestep for each subdomain
    computeTimeStepWidth();

    // use the minimum as global timestep
    double dt_global;
    MPI_Allreduce(&dt_local, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    // if necessary adapt so that every full second is reached
    if (std:floor(currentTime + dt) == std::floor(currentTime) + 1)
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
    
    for ( int i = 0; i < discretization_->nCells[0]; i++)
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
    
    for ( int i = 0; i < discretization_->nCells[0]; i++)
    {
        // u
        discretization_->u(i, discretization_->nCells[1]) = 2. * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->nCells[1] - 1);
        // v
        discretization_->v(i, discretization_->nCells[1]) = settings_.dirichletBcTop[1];
        // f
        discretization_->f(i, discretization_->nCells[1]) = discretization_->u(i, discretization_->nCells[1] - 1);
        // g
        discretization_->g(i, discretization_->nCells[1]) = discretization_->v(i, discretization_->nCells[1] - 1);
    }
}

void ComputationParallel::uvExchangeHorizontal()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even processes: send left, receive left, send right, receive right
    if ((partitioning_.ownRankNo() % 2) == 0)
    {
        // left
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_left[j+2] = discretization_->u(0, j); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_left;

            MPI_Send(&u_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(-2, j) = other_ghost_to_u_left[j+2]; // u_{-2}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_left[j+2] = discretization_->v(0, j); // v_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_left;

            MPI_Send(&v_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(-1, j) = other_ghost_to_v_left[j+2]; // v_{-1}
            }
        }

        // right
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_right[j+2] = discretization_->u(discretization_->nCells()[1] -2, j); // u_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_right;

            MPI_Send(&u_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(discretization_->nCells()[1], j) = other_ghost_to_u_right[j+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_right[j+2] = discretization_->v(discretization_->nCells()[1] -1, j); // v_{n-1}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_right;

            MPI_Send(&v_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(discretization_->nCells()[1], j) = other_ghost_to_v_right[j+2]; // v_{n}
            }
        }
    }
    // the uneven processes: receive right, send right, receive left, send left
    if ((partitioning_.ownRankNo() % 2) == 1)
    {
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_right[j+2] = discretization_->u(discretization_->nCells()[1] -2, j); // u_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_right;

            MPI_Recv(&other_ghost_to_u_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(discretization_->nCells()[1], j) = other_ghost_to_u_right[j+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_right[j+2] = discretization_->v(discretization_->nCells()[1] -1, j); // v_{n-1}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_right;

            MPI_Recv(&other_ghost_to_v_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(discretization_->nCells()[1], j) = other_ghost_to_v_right[j+2]; // v_{n}
            }
        }
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_left[j+2] = discretization_->u(0, j); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_left;

            MPI_Recv(&other_ghost_to_u_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(-2, j) = other_ghost_to_u_left[j+2]; // u_{-2}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_left[j+2] = discretization_->v(0, j); // v_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_left;

            MPI_Recv(&other_ghost_to_v_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(-1, j) = other_ghost_to_v_left[j+2]; // v_{-1}
            }
        }
    }
}

ComputationParallel::uvExchangeVertical()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even row processes: send top, receive top, send bottom, receive bottom
    if ((partitioning_.ownRankCoordinate()[1] % 2) == 0)
    {
        // top
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_top[i+2] = discretization_->u(i, discretization_->nCells()[0] -1); // u_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_u_top;

            MPI_Send(&u_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, discretization_->nCells()[0]) = other_ghost_to_u_top[i+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_top[i+2] = discretization_->v(i,discretization_->nCells()[0] -2); // v_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_top;

            MPI_Send(&v_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, discretization_->nCells()[0]) = other_ghost_to_v_top[i+2]; // v_{n}
            }
        }
        // bottom
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_bottom[i+2] = discretization_->u(i, 0); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_bottom;

            MPI_Send(&u_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, -1) = other_ghost_to_u_bottom[j+2]; // u_{-1}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_bottom[j+2] = discretization_->v(i, 0); // v_{0}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_v_bottom;

            MPI_Send(&v_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, -2) = other_ghost_to_v_bottom[i+2]; // v_{-2}
            }
        }
    }
    // the uneven row processes: receive bottom, send bottom, receive top, send top, 
    if ((partitioning_.ownRankCoordinate()[1] % 2) == 1)
    {
        // bottom
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_bottom[i+2] = discretization_->u(i, 0); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_bottom;

            MPI_Recv(&other_ghost_to_u_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, -1) = other_ghost_to_u_bottom[j+2]; // u_{-1}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_bottom[j+2] = discretization_->v(i, 0); // v_{0}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_v_bottom;

            MPI_Recv(&other_ghost_to_v_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, -2) = other_ghost_to_v_bottom[i+2]; // v_{-2}
            }
        }
        // top
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_top[i+2] = discretization_->u(i, discretization_->nCells()[0] -1); // u_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_u_top;

            MPI_Recv(&other_ghost_to_u_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, discretization_->nCells()[0]) = other_ghost_to_u_top[i+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_top[i+2] = discretization_->v(i,discretization_->nCells()[0] -2); // v_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_top;

            MPI_Recv(&other_ghost_to_v_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, discretization_->nCells()[0]) = other_ghost_to_v_top[i+2]; // v_{n}
            }
        }
    }
}

void ComputationParallel::uvExchangeHorizontal()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even processes: send left, receive left, send right, receive right
    if ((partitioning.ownRankNo() % 2) == 0)
    {
        // left
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_left[j+2] = u(0, j); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_left;

            MPI_Send(&u_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(-2, j) = other_ghost_to_u_left[j+2]; // u_{-2}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_left[j+2] = v(0, j); // v_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_left;

            MPI_Send(&v_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(-1, j) = other_ghost_to_v_left[j+2]; // v_{-1}
            }
        }

        // right
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_right[j+2] = u(discretization_->nCells()[1] -2, j); // u_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_right;

            MPI_Send(&u_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(discretization_->nCells()[1], j) = other_ghost_to_u_right[j+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_right[j+2] = v(discretization_->nCells()[1] -1, j); // v_{n-1}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_right;

            MPI_Send(&v_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(discretization_->nCells()[1], j) = other_ghost_to_v_right[j+2]; // v_{n}
            }
        }
    }
    // the uneven processes: receive right, send right, receive left, send left
    if ((partitioning.ownRankNo() % 2) == 1)
    {
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_right[j+2] = u(discretization_->nCells()[1] -2, j); // u_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_right;

            MPI_Recv(&other_ghost_to_u_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(discretization_->nCells()[1], j) = other_ghost_to_u_right[j+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_right[j+2] = v(discretization_->nCells()[1] -1, j); // v_{n-1}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_right;

            MPI_Recv(&other_ghost_to_v_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(discretization_->nCells()[1], j) = other_ghost_to_v_right[j+2]; // v_{n}
            }
        }
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[1]+4> u_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                u_to_other_ghost_left[j+2] = u(0, j); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_left;

            MPI_Recv(&other_ghost_to_u_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->u(-2, j) = other_ghost_to_u_left[j+2]; // u_{-2}
            }

            // v
            std::array<double, discretization_->nCells()[1]+4> v_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                v_to_other_ghost_left[j+2] = v(0, j); // v_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_left;

            MPI_Recv(&other_ghost_to_v_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->v(-1, j) = other_ghost_to_v_left[j+2]; // v_{-1}
            }
        }
    }
}

ComputationParallel::uvExchangeVertical()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even processes: send top, receive top, send bottom, receive bottom
    if ((partitioning.ownRankNo() % 2) == 0)
    {
        // top
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_top[i+2] = u(i, discretization_->nCells()[0] -1); // u_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_u_top;

            MPI_Send(&u_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, discretization_->nCells()[0]) = other_ghost_to_u_top[i+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_top[i+2] = v(i,discretization_->nCells()[0] -2); // v_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_top;

            MPI_Send(&v_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, discretization_->nCells()[0]) = other_ghost_to_v_top[i+2]; // v_{n}
            }
        }
        // bottom
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_bottom[i+2] = u(i, 0); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_bottom;

            MPI_Send(&u_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_u_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, -1) = other_ghost_to_u_bottom[j+2]; // u_{-1}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_bottom[j+2] = v(i, 0); // v_{0}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_v_bottom;

            MPI_Send(&v_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_v_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, -2) = other_ghost_to_v_bottom[i+2]; // v_{-2}
            }
        }
    }
    // the uneven processes: receive bottom, send bottom, receive top, send top, 
    if ((partitioning.ownRankNo() % 2) == 1)
    {
        // bottom
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_bottom[i+2] = u(i, 0); // u_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_u_bottom;

            MPI_Recv(&other_ghost_to_u_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, -1) = other_ghost_to_u_bottom[j+2]; // u_{-1}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_bottom[j+2] = v(i, 0); // v_{0}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_v_bottom;

            MPI_Recv(&other_ghost_to_v_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, -2) = other_ghost_to_v_bottom[i+2]; // v_{-2}
            }
        }
        // top
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> u_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                u_to_other_ghost_top[i+2] = u(i, discretization_->nCells()[0] -1); // u_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_u_top;

            MPI_Recv(&other_ghost_to_u_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Send(&u_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->u(i, discretization_->nCells()[0]) = other_ghost_to_u_top[i+2]; // u_{n}
            }

            // v
            std::array<double, discretization_->nCells()[0]+4> v_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                v_to_other_ghost_top[i+2] = v(i,discretization_->nCells()[0] -2); // v_{n-2}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_v_top;

            MPI_Recv(&other_ghost_to_v_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Send(&v_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->v(i, discretization_->nCells()[0]) = other_ghost_to_v_top[i+2]; // v_{n}
            }
        }
    }
}

void ComputationParallel::pExchangeHorizontal()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even processes: send left, receive left, send right, receive right
    if ((partitioning.ownRankNo() % 2) == 0)
    {
        // left
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // p
            std::array<double, discretization_->nCells()[1]+4> p_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                p_to_other_ghost_left[j+2] = p(0, j); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_left;

            MPI_Send(&p_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(-1, j) = other_ghost_to_p_left[j+2]; // p_{-1}
            }
        }

        // right
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // p
            std::array<double, discretization_->nCells()[1]+4> p_to_other_ghost_right;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                p_to_other_ghost_right[j+2] = p(discretization_->nCells()[1] -1, j); // p_{n-1}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_right;

            MPI_Send(&p_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(discretization_->nCells()[1], j) = other_ghost_to_p_right[j+2]; // p_{n}
            }
        }
    }
    // the uneven processes: receive right, send right, receive left, send left
    if ((partitioning.ownRankNo() % 2) == 1)
    {
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
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
                     partitioning_.ownRightNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_right, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownRightNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(discretization_->nCells()[1], j) = other_ghost_to_p_right[j+2]; // p_{n}
            }
        }
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // p
            std::array<double, discretization_->nCells()[1]+4> p_to_other_ghost_left;
            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                p_to_other_ghost_left[j+2] = p(0, j); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_left;

            MPI_Recv(&other_ghost_to_p_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_left, 
                     discretization_->nCells()[1]+4,
                     MPI_DOUBLE,
                     partitioning_.ownLeftNeighbour(),
                     )

            for (int j = -2; j < discretization_->nCells()[1] +2; j++)
            {
                discretization_->p(-1, j) = other_ghost_to_p_left[j+2]; // p_{-1}
            }
        }
    }
}

ComputationParallel::pExchangeVertical()
{
    // the +4 in nCells+4 is from the consistent number of halo-cells in the staggered grid

    // the even processes: send top, receive top, send bottom, receive bottom
    if ((partitioning.ownRankNo() % 2) == 0)
    {
        // top
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // p
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_top[i+2] = p(i, discretization_->nCells()[0] -1); // p_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_p_top;

            MPI_Send(&p_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, discretization_->nCells()[0]) = other_ghost_to_p_top[i+2]; // p_{n}
            }
        }
        // bottom
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_bottom[i+2] = p(i, 0); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_bottom;

            MPI_Send(&p_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Recv(&other_ghost_to_p_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, -1) = other_ghost_to_p_bottom[j+2]; // p_{-1}
            }
        }
    }
    // the uneven processes: receive bottom, send bottom, receive top, send top, 
    if ((partitioning.ownRankNo() % 2) == 1)
    {
        // bottom
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_bottom;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_bottom[i+2] = p(i, 0); // p_{0}
            }
            std::array<double, discretization_->nCells()[1]+4> other_ghost_to_p_bottom;

            MPI_Recv(&other_ghost_to_p_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_bottom, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownBottomNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, -1) = other_ghost_to_p_bottom[j+2]; // p_{-1}
            }
        }
        // top
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // u
            std::array<double, discretization_->nCells()[0]+4> p_to_other_ghost_top;
            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                p_to_other_ghost_top[i+2] = p(i, discretization_->nCells()[0] -1); // p_{n-1}
            }
            std::array<double, discretization_->nCells()[0]+4> other_ghost_to_p_top;

            MPI_Recv(&other_ghost_to_p_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )
            MPI_Send(&p_to_other_ghost_top, 
                     discretization_->nCells()[0]+4,
                     MPI_DOUBLE,
                     partitioning_.ownTopNeighbour(),
                     )

            for (int i = -2; i < discretization_->nCells()[0] +2; i++)
            {
                discretization_->p(i, discretization_->nCells()[0]) = other_ghost_to_p_top[i+2]; // p_{n}
            }
    }
}