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
    
    pressureSolver_ = std::make_unique<RedBlack>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations); 

    // misc
    dt_ = 0.;

    // initialize output writer
    outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_);
    outputWriterTextParallel_     = std::make_unique<OutputWriterTextParallel>(discretization_);
}

void ComputationParallel::runSimulation()
{
    double currentTime = 0;
    int sec = 0;

    std::cout << "+++++++++++++++++++++++" << std::endl;
    std::cout << "Starting at time: " << currentTime << std::endl;
    std::cout << "+++++++++++++++++++++++" << std::endl;

    // the steps correspond to the steps in our algorithm in the overleaf or docs/numsim-algos.tex
    while (currentTime < settings_.endTime)
    {
        // step 1: set the boundary values
        applyBoundaryValues();
        std::cout << "Applied boundary values for u/v and F/G." << std::endl;

        // step 2: compute time step width
        computeTimeStepWidth();

        std::cout << "Computed time step width: " << dt_ << std::endl;

        // set the dt such that the simulation stops exactly on endTime
        if (currentTime + dt_ > settings_.endTime)
        {
            dt_ = settings_.endTime - currentTime;

            std::cout << std::endl;
            std::cout << "Final time step!" << std::endl;
        } 
        // set the dt such that the simulation steps exactly over full multiples of 1s (for the output)
        if (std:floor(currentTime + dt) == sec +1.)
        {
            dt_ = (double) sec - currentTime;
        }
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
        if (std:floor(currentTime + dt) == sec +1.)
        {
            outputWriterParaview_->writeFile(currentTime);
            outputWriterText_->writeFile(currentTime);
            sec = sec +1;
        }

    }
    
    // end the MPI-session
    MPI_Finalize()
}

void ComputationParallel::computeTimeStepWidth()
{
    double boundary_diffusion = 0.;

    // boundary from diffusion
    if (meshWidth_[0] == meshWidth_[1])
    {
        boundary_diffusion = (settings_.re * meshWidth_[0] * meshWidth_[1]) / 4.;
    } else
    {
        double h2x = meshWidth_[0] * meshWidth_[0];
        double h2y = meshWidth_[1] * meshWidth_[1];
        boundary_diffusion = (settings_.re / 2.) * (h2x * h2y) * (1./(h2x + h2y));
      
        std::cout << settings_.re << " " << h2x << " " << h2y << std::endl;
    }

    // calculate max absolute velocities
    double u_max = 0.;
    for ( int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++)
    { 
        for ( int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++)
        {
            double value = std::abs(discretization_->u(i,j));
            if (value > u_max)
            {
                u_max = value;
            }
        }
    }

    double v_max = 0.;
    for ( int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++)
    { 
        for ( int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++)
        {
            double value = std::abs(discretization_->v(i,j));
            if (value > v_max)
            {
                v_max = value;
            }
        }
    }

    // boundary from convection
    double boundary_convection_u = meshWidth_[0] / u_max;
    double boundary_convection_v = meshWidth_[1] / v_max;

    // together
    double min_dt = std::min({boundary_diffusion, boundary_convection_u, boundary_convection_v});

    std::cout << "dt boundaries - diffusion: " << boundary_diffusion << " convection_u: " << boundary_convection_u << " convection_v: " << boundary_convection_v << std::endl;

    // security factor
    double dt_local =  std::min(min_dt * settings_.tau, settings_.maximumDt);

    // minimize by communication with the other processes
    double dt_global;
    MPI_Allreduce(&dt_local, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    dt_ = dt_global;
}
    
void ComputationParallel::applyBoundaryValues()
{
    // first horizontal pairwise exchange
    // the even processes: left, then right
    if ((partitioning.ownRankNo() % 2) == 0)
    {
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // uExchangeVertical
            // vExchangeVertical
            // pExchangeVertical
        }
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // uExchangeVertical
            // vExchangeVertical
            // pExchangeVertical
        }
    }
    // the uneven processes: right, then left
    if ((partitioning.ownRankNo() % 2) == 1)
    {
        if (partitioning_.ownPartitionContainsRightBoundary())
        {
            applyBoundaryValuesRight();
        } else
        {
            // uExchangeVertical
            // vExchangeVertical
            // pExchangeVertical
        }
        if (partitioning_.ownPartitionContainsLeftBoundary())
        {
            applyBoundaryValuesLeft();
        } else 
        {
            // uExchangeVertical
            // vExchangeVertical
            // pExchangeVertical
        }
    }

    // then vertical pairwise exchange
    // the even processes: top, then bottom
    if ((partitioning.ownRankNo() % 2) == 0)
    {
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // uExchangeHorizontal
            // vExchangeHorizontal
            // pExchangeHorizontal
        }
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // uExchangeHorizontal
            // vExchangeHorizontal
            // pExchangeHorizontal
        }
    }
    // the uneven processes: bottom, then top
    if ((partitioning.ownRankNo() % 2) == 1) 
    {
        if (partitioning_.ownPartitionContainsBottomBoundary())
        {
            applyBoundaryValuesBottom();
        } else
        {
            // uExchangeHorizontal
            // vExchangeHorizontal
            // pExchangeHorizontal
        }
        if (partitioning_.ownPartitionContainsTopBoundary())
        {
            applyBoundaryValuesTop();
        } else 
        {
            // uExchangeHorizontal
            // vExchangeHorizontal
            // pExchangeHorizontal
        }
    }
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
    // u
    for ( int i = 0; i < discretization_->uIEnd(); i++)
    {
        discretization_->u(i, discretization_->uJBegin())   = 2. * settings_.dirichletBcBottom[0] - discretization_->u(i, discretization_->uJBegin() +1);
    }
    
    // v
    for ( int i = 0; i < discretization_->vIEnd(); i++)
    {
        discretization_->v(i, discretization_->vJBegin())  = settings_.dirichletBcBottom[1];
    }
    
    // F
    for ( int i = 0; i < discretization_->uIEnd(); i++)
    {
        discretization_->f(i, discretization_->uJBegin()) = discretization_->u(i, discretization_->uJBegin());
    }

    // set g
    for ( int i = 0; i < discretization_->vIEnd(); i++)
    {
        discretization_->g(i, discretization_->vJBegin()) = discretization_->v(i, discretization_->vJBegin());
    }

}

void ComputationParallel::applyBoundaryValuesTop()
{
    // set boundaries only in domain as corners belong to size, domain begins at idx 0
    
    // u
    for ( int i = 0; i < discretization_->uIEnd()-1; i++)
    {
        discretization_->u(i, discretization_->uJEnd() - 1) = 2. * settings_.dirichletBcTop[0] - discretization_->u(i, discretization_->uJEnd() -2);
    }
    
    // set v
    for ( int i = 0; i < discretization_->vIEnd()-1; i++)
    {
        discretization_->v(i, discretization_->vJEnd() - 1) = settings_.dirichletBcTop[1];
    }

    // set f
    for ( int i = 0; i < discretization_->uIEnd()-1; i++)
    {
        discretization_->f(i, discretization_->uJEnd() -1) = discretization_->u(i, discretization_->uJEnd() -1);
    }
    // set g
    for ( int i = 0; i < discretization_->vIEnd()-1; i++)
    {
        discretization_->g(i, discretization_->vJEnd() -1) = discretization_->v(i, discretization_->vJEnd() -1);
    }
}
