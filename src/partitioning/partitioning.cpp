#include "partitioning.h"

Partitioning::Partitioning(int ownRankNo, int nRanks, std::array<int,2> nCells) :
ownRankNo_(ownRankNo),
nRanks_(nRanks),
nCellsGlobal_(nCells)
{
    // divide computational domain into subdomains
    factorizeSubdomains();  // sets nSubdomains
    
    int n_subd = nSubdomains_[0]; // number of subdomains in i direction or nCells[0]
    int m_subd = nSubdomains_[1]; // number of subdomains in j direction or nCells[1]
    
    // map ownRank -> subdomain
    int process_column = ownRankNo_ % n_subd ;
    int process_row = std::floor(ownRankNo_ / n_subd);

    // save as rank coordinate
    ownRankCoordinate_ = {process_column, process_row};

    // get neighbors
    // bottom border
    if (process_row == 0)
    {
        partitionNeighbours_[0] = MPI_PROC_NULL;
    } else
    {
        partitionNeighbours_[0] = ownRankNo_ - n_subd;
    }
    // left border
    if (process_column == 0)
    {
        partitionNeighbours_[1] = MPI_PROC_NULL;
    } else
    {
        partitionNeighbours_[1] = ownRankNo_ -1;
    }
    // top border
    if (process_row == (m_subd - 1))
    {
        partitionNeighbours_[3] = MPI_PROC_NULL;
    } else
    {
        partitionNeighbours_[3] = ownRankNo_ + n_subd;
    }
    // right border
    if (process_column == (n_subd - 1))
    {
        partitionNeighbours_[3] = MPI_PROC_NULL;
    } else
    {
        partitionNeighbours_[3] = ownRankNo_ +1;
    }

    nCellsLocal_ = {nCellsGlobal_[0] / n_subd, nCellsGlobal_[1] / m_subd};

    std::cout << "computed nCellsLocal: " << nCellsLocal_[0] << nCellsLocal_[1]  << std::endl;
    
    // check if nCellsLocal is int, else throw an error
    if (nCellsLocal_[0] != std::floor(nCellsLocal_[0]) || nCellsLocal_[1] != std::floor(nCellsLocal_[1])) 
    {
        std::cout << "nCellsLocal is not an array of integers, computed partition not suited for global domain size." << std::endl;
        throw;
    }
    
    nodeOffset_ = {nCellsLocal_[0] * process_column, nCellsLocal_[1] * process_row};
}

void Partitioning::factorizeSubdomains() {

    // init factorization with maximum cost
    int cost_opt = nCellsGlobal_[0]*(nCellsGlobal_[1]-1) + (nCellsGlobal_[0]-1)*nCellsGlobal_[1]; // communication cost - number inner edges 
    
    // save best combinations here (minimizes commuication costs)
    int n_opt = 1;         // number of subdomains in i direction
    int m_opt = nRanks_;   // number of subdomains in j direction

    int m = nRanks_; // temporary factor
    // iterate over all possible factorizations 
    for ( int n = 1; n < nRanks_+1; n++)    
    {
        if ((nRanks_ % n) == 0) // if nRanks can be devided in n * (nRanks_/n)
        {
            m = nRanks_ / n;
            int cost = nCellsGlobal_[0]*(m-1) + nCellsGlobal_[1]*(n-1); // number of innner edges with partition n x m
            
            if ( cost < cost_opt) {
                n_opt = n; 
                m_opt = m; 
                cost_opt = cost;
            }
        }
    }

    // save partition
    nSubdomains_[0] = n_opt;
    nSubdomains_[1] = m_opt;
    
    std::cout << "Computed optimal partition in subdomains: " << nSubdomains_[0] << nSubdomains_[1]  << " with costs " << cost_opt << std::endl;
}

std::array<int, 2> Partitioning::nodeOffset() {
    return nodeOffset_;
}

std::array<int, 2> Partitioning::nCellsGlobal() {
    return nCellsGlobal_;
}

int Partitioning::ownRankNo() {
    return ownRankNo_;
}

std::array<int, 2> Partitioning::ownRankCoordinate() {
    return ownRankCoordinate_;
}

// returns whether partition contains boundary
bool Partitioning::ownPartitionContainsBottomBoundary() {
    return partitionNeighbours_[0] == MPI_PROC_NULL;
}

bool Partitioning::ownPartitionContainsLeftBoundary() {
    return partitionNeighbours_[1] == MPI_PROC_NULL;
}

bool Partitioning::ownPartitionContainsTopBoundary() {
    return partitionNeighbours_[2] == MPI_PROC_NULL;
}

bool Partitioning::ownPartitionContainsRightBoundary() {
    return partitionNeighbours_[3] == MPI_PROC_NULL;
}

// returns the rank numbers
int Partitioning::ownBottomNeighbour() {
    return partitionNeighbours_[0];
}

int Partitioning::ownLeftNeighbour() {
    return partitionNeighbours_[1];
}

int Partitioning::ownTopNeighbour() {
    return partitionNeighbours_[2];
}

int Partitioning::ownRightNeighbour() {
    return partitionNeighbours_[3];
}
