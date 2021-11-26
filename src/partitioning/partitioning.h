#pragma once

#include <memory>
#include <cmath>
#include <mpi.h>

class Partitioning
{
public: 

    //! construct partitioning
    Partitioning(int ownRankNo, int nRanks, std::array<int,2> nCells);

    //! return the global node offset of lower left corner cell
    const std::array<int,2> nodeOffset() const;

    //! number of global cells
    const std::array<int, 2> nCellsGlobal() const;

    //! return own rank number
    const int ownRankNo() const;

    //! return true if partition constains right boundary
    const bool ownPartitionContainsRightBoundary() const;
    //! return true if partition contains left boundary
    const bool ownPartitionContainsLeftBoundary() const;
    //! return true if partition  constains top boundary
    const bool ownPartitionContainsTopBoundary() const;
    //! return true if partition contains bottom boundary
    const bool ownPartitionContainsBottomBoundary() const;

    //! return rank number of right neighbour, -1 if boundary
    const int ownRightNeighbour() const;
    //! return rank number of left neighbour, -1 if boundary
    const int ownLeftNeighbour() const;
    //! return rank number of top neighbour, -1 if boundary
    const int ownTopNeighbour() const;
    //! return rank number of bottom neighbour, -1 if boundary
    const int ownBottomNeighbour() const;

protected:

    //! global node offset of lower left cell
    const std::array<int, 2> nodeOffset_;
    //! number of global cells
    const std::array<int, 2> nCellsGlobal_;
    //! number of local cells
    const std::array<int, 2> nCellsLocal_;

    //! own rank
    const int ownRankNo_;
    //! global number of ranks
    const int nRanks_;

    //! return rank number of partition neighbours (>0) 
    //! or MPI_PROC_NULL if edge lies on domain boundary
    //! array order <bottom, left, top, right> 
    const std::array<int, 4> partitionNeighbours_;
};
