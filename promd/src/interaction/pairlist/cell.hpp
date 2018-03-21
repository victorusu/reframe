#ifndef CELL_HPP
#define CELL_HPP

#include "../../main.hpp"

struct Cell {

    int id;

    int nAtoms;                             // number of atoms in this cell
    // int nCgs;                               // number of charge groups in this cell
    int nNeighCells;                        // number of neighbors of this cell
    // int nGhostCells;                        // number of neighbors of this cell
    std::vector<int> atomList;              // list of atom indices
    // std::vector<int> cgList;                // list of charge group indices
    std::vector<int> neighCellListId;       // neigh cell list
    // std::vector<int> mpiRanksToCommunicate; // neigh cell list
    // std::vector<int> ghostCells;         // ghost cell list

    // bool updateNeighList;

    bool operator==(Cell cell)
    {
        if(this->id == cell.id)
            return true;
        else
            return false;
    }

    bool operator==(int cellId)
    {
        if(this->id == cellId)
            return true;
        else
            return false;
    }

    template<typename FileType>
    void print(FileType &outfile)
    {
        outfile.printf("    cell:\n");
        outfile.printf("  natoms: %d\n", nAtoms);
        outfile.printf("atomList:\n");

        // return;
        for(int i = 0; i < this->nAtoms; i++)
            outfile.printf("atomList[%d]: %d\n", i, atomList[i]);

        outfile.printf("neighCellListId: %d\n", this->neighCellListId.size());
        for(size_t i = 0; i < this->neighCellListId.size(); i++)
            outfile.printf("neighCellListId[%06zu]: %d\n", i, neighCellListId[i]);
    }


};

bool operator==(int cellId, Cell cell);

#endif //CELL_HPP