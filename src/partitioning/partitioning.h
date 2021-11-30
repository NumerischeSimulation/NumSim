#pragma once

#include <memory>
#include <cmath>
#include <mpi.h>
#include <iostream>

class Partitioning
{
    public: 

        //! construct partitioning
        Partitioning(int ownRankNo, int nRanks, std::array<int,2> nCells);

        //! return the global node offset of lower left corner cell
        std::array<int,2> nodeOffset() ;

        //! number of global cells
        std::array<int, 2> nCellsGlobal() ;

        //! return own rank number
        int ownRankNo() ;

        //! return own 2d rank coordinate on grid of subdomains
        std::array<int,2> ownRankCoordinate() ;

        //! return true if partition constains right boundary
        bool ownPartitionContainsRightBoundary();
        //! return true if partition contains left boundary
        bool ownPartitionContainsLeftBoundary();
        //! return true if partition  constains top boundary
        bool ownPartitionContainsTopBoundary();
        //! return true if partition contains bottom boundary
        bool ownPartitionContainsBottomBoundary();

        //! return rank number of right neighbour, -1 if boundary
        int ownRightNeighbour();
        //! return rank number of left neighbour, -1 if boundary
        int ownLeftNeighbour();
        //! return rank number of top neighbour, -1 if boundary
        int ownTopNeighbour();
        //! return rank number of bottom neighbour, -1 if boundary
        int ownBottomNeighbour();

    protected:

        //! global node offset of lower left cell
        std::array<int, 2> nodeOffset_;
        //! number of global cells
        const std::array<int, 2> nCellsGlobal_;
        //! number of local cells
        std::array<int, 2> nCellsLocal_;

        //! own rank
        const int ownRankNo_;
        //! global number of ranks
        const int nRanks_;

        //! own rank coordinates on 2D grid of subdomains
        std::array<int, 2> ownRankCoordinate_;

        //! return rank number of partition neighbours (>0) 
        //! or MPI_PROC_NULL if edge lies on domain boundary
        //! array order <bottom, left, top, right> 
        std::array<int, 4> partitionNeighbours_;

    private: 
        //! finds the optimal factorization of subdomains given nRanks_
        void factorizeSubdomains() ; 

        //! number of subdomains in each direction
        std::array<int, 2> nSubdomains_;
};
