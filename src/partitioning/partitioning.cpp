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
    std::cout << "computed nCellsLocal: " << nCellsLocal_ << std::endl;
    
    // check if nCellsLocal is int, else throw an error
    if (nCellsLocal_ != floor(nCellsLocal)) 
    {
        std::cout << "nCellsLocal is not an array of integers, computed partition not suited for global domain size." << std::endl;
        throw;
    }
    
    nodeOffset_ = {nCellsLocal_[0] * process_column, nCellsLocal_[1] * process_row}
}

Partitioning::factorizeSubdomains() {

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
            cost = nCellsGlobal_[0]*(m-1) + nCellsGlobal_[1]*(n-1); // number of innner edges with partition n x m
            
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
    
    std::cout << "Computed optimal partition in subdomains: " << nSubdomains_ << " with costs " >> cost_opt << std::endl;
}
