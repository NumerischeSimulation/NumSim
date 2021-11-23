#pragma once

class DataTransfer
{
    DataTransfer(std::array<int,2> nLocalCells);

    //! exchange ghost layer vertically to the left and then to the right
    void uExchangeVertical(FieldVariable localULeft, FieldVariable localURight); // TODO

    protected:
    int nLocalCellsX_;
    int nLocalCellsY_;
};