#include "partitioning.h"

Partitioning::Partitioning(int ownRankNo, int nRanks, std::array<int,2> nCells) :
ownRankNo_(ownRankNo),
nRanks_(nRanks),
nCellsGlobal_(nCells)
{
    // map ownRank -> subdomain

    // TODO: more than two rows
    // factorize N
    // minimize cost for each factorization

    // only for two rows
    int n_subd = 2; // number of subdomains in i direction or nCells[0]
    int m_subd = 2; // number of subdomains in j direction or nCells[1]
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
    nodeOffset_ = {nCellsLocal_[0] * process_column, nCellsLocal_[1] * process_row}
}
