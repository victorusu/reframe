#ifndef PAIRLIST_PERATOM_DOMAIN_CPU_HPP
#define PAIRLIST_PERATOM_DOMAIN_CPU_HPP

#include "../../main.hpp"
#include "../../domain/domain.hpp"
#include "../../utilities/stdout.hpp"

#include "cell.hpp"
#include "../../configuration/conf.hpp"
#include "../../simulation/box.hpp"
#include "../../simulation/simparam.hpp"
#include "../../topology/topology.hpp"

#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

struct NeighList {
    int nAtoms;
    int atomId;
    std::vector<int> listIds;    // list of atom indices
};

struct PairList {

    std::vector<Cell> cells;

    std::vector<NeighList> vdWNeigh;
    std::vector<NeighList> coulNeigh;
    std::vector<NeighList> vdWCoulNeigh;

    // std::vector<Cell> clist;
    // std::vector<int> headList;
    // std::vector<int> cList;

    // std::vector<number> posX;
    // std::vector<number> posY;
    // std::vector<number> posZ;
    number maxDisplacement;
    RVector refPosition;

    bool doMaxDisplacement;

    VHR_ALWAYS_INLINE
    int estimateNNeighAtoms(const number nAtoms, Box &simBox, const number rlist)
    {
        // 4.2 is approx 4./3. * M_PI
        // This aprroximation is not a problem per se, for we are
        // just estimating the number of neighbours
        number cutVol = 4.2 * rlist * rlist * rlist;

        number boxVol = simBox.volume();
        number density = nAtoms / boxVol;

        return ceil(density * cutVol);
    }

    VHR_ALWAYS_INLINE
    int estimateNAtomsPerCells(const int nAtoms, const int ncell)
    {
        // int nidx = 2*nAtoms / ncell + 2;
        // nidx = ((nidx/2) + 1) * 2;
        // return nidx;
        const int nidx = 2*nAtoms / ncell + 2;
        return ((nidx/2) + 1) * 2;
    }

    VHR_ALWAYS_INLINE
    void addAtomJToNeighListOfI(const int ii, const int jj, const Topology &topology, NeighList & neighList)
    {
        // startTimer(addAtomJToNeighListOfIPairlistTimer);

        // if(notexcluded(topology.excl, ii, jj-ii-1) && notexcluded(topology.excl, jj, ii-jj-1)) {
        int i = ii;
        int j = jj;
        if(i > j) {
            std::swap(i,j);
        }
        if((j-i > 31) || notexcluded(topology.excl, i, j-i-1)) {

            const int nAtoms = neighList.nAtoms;

            if(nAtoms < neighList.listIds.size()) {
                // cout.printf("firstIF, address %d\n", &neighList);
                neighList.listIds[nAtoms] = jj;
                neighList.nAtoms++;
            }
            else {
                // cout.printf("secondIF, address %d\n", &neighList);
                neighList.listIds.push_back(jj);
                neighList.nAtoms++;
            }
        }
        // else {
        //     cerr.printf("atom %d is excluded from %d\n", i, j);
        // }

        // addToTime(addAtomJToNeighListOfIPairlistTimer, addAtomJToNeighListOfIPairlistTime);
    }

    /*
     * This function should always be implemented as thread safe
     *
     */
    VHR_INLINE
    void updateAtomNeighListOfLocalCells(const Topology &topology)
    {
        startTimer(updateAtomNeighListOfLocalCellsPairlistTimer);

        const int nCells = this->cells.size();
        int localCellsIndex;

        #pragma omp parallel for schedule(dynamic) private(localCellsIndex) shared(topology)
        for(localCellsIndex = 0; localCellsIndex < nCells; localCellsIndex++) {

            Cell *c1 = this->cells.data() + localCellsIndex;

            // skip empty cells
            // this should be a rare case
            // should happen for surface tension calculations
            if(c1->nAtoms < 1)
                continue;

            int i, j, k, ii, jj;
            bool iiHasVdW, iiHasCoul, jjHasVdW, jjHasCoul;

            for (i = 0; i < c1->nAtoms; i++) {

                ii = c1->atomList[i];

                iiHasVdW = hasVdW(ii, topology);
                iiHasCoul= hasCoul(ii, topology);

                NeighList &vdWCoulNeigh = this->vdWCoulNeigh[ii];
                NeighList &vdWNeigh = this->vdWNeigh[ii];
                NeighList &coulNeigh = this->coulNeigh[ii];

                if(iiHasVdW && iiHasCoul) {

                    cerr.printf("we should not have vdw and coul\n");

                    // self interaction of atoms in cell
                    for(j = i + 1; j < c1->nAtoms; j++) {

                        jj = c1->atomList[j];

                        // we check for vdw and coul interactions
                        jjHasVdW = hasVdW(jj, topology);
                        jjHasCoul = hasCoul(jj, topology);

                        // but j can have only one of the two
                        // so we prompt j also
                        if(jjHasVdW && jjHasCoul)
                            addAtomJToNeighListOfI(ii, jj, topology, vdWCoulNeigh);
                        else if(jjHasVdW)
                            addAtomJToNeighListOfI(ii, jj, topology, vdWNeigh);
                        else if(jjHasCoul)
                            addAtomJToNeighListOfI(ii, jj, topology, coulNeigh);
                    }

                    // interaction with atoms in neighCells
                    for (std::vector<int>::iterator it = c1->neighCellListId.begin(); it != c1->neighCellListId.end(); it++) {

                        const Cell *c2 = this->cells.data() + *it;

                        // skip empty cells
                        // this should be a rare case
                        // should happen for surface tension calculations
                        if(c2->nAtoms < 1)
                            continue;

                        for(j = 0; j < c2->nAtoms; j++) {

                            jj = c2->atomList[j];

                            // we check for vdw and coul interactions
                            jjHasVdW = hasVdW(jj, topology);
                            jjHasCoul = hasCoul(jj, topology);

                            // but j can have only one of the two
                            // so we prompt j also
                            if(jjHasVdW && jjHasCoul)
                                addAtomJToNeighListOfI(ii, jj, topology, vdWCoulNeigh);
                            else if(jjHasVdW)
                                addAtomJToNeighListOfI(ii, jj, topology, vdWNeigh);
                            else if(jjHasCoul)
                                addAtomJToNeighListOfI(ii, jj, topology, coulNeigh);
                        }
                    }

                } else if(iiHasVdW) {

                    // self interaction of atoms in cell
                    for(j = i + 1; j < c1->nAtoms; j++) {

                        jj = c1->atomList[j];

                        // we check only for vdW interactions
                        // because ii has no coul anyway
                        jjHasVdW = hasVdW(jj, topology);

                        // but jj can have only have vdw interactions
                        // to interact with ii
                        if(jjHasVdW)
                            addAtomJToNeighListOfI(ii, jj, topology, vdWNeigh);

                    }

                    // interaction with atoms in neighCells
                    for (std::vector<int>::iterator it = c1->neighCellListId.begin(); it != c1->neighCellListId.end(); it++) {

                        const Cell *c2 = this->cells.data() + *it;

                        // skip empty cells
                        // this should be a rare case
                        // should happen for surface tension calculations
                        if(c2->nAtoms < 1)
                            continue;

                        for(j = 0; j < c2->nAtoms; j++) {

                            jj = c2->atomList[j];

                            // we check only for vdW interactions
                            // because ii has no coul anyway
                            jjHasVdW = hasVdW(jj, topology);

                            // but jj can have only have vdw interactions
                            // to interact with ii
                            if(jjHasVdW)
                                addAtomJToNeighListOfI(ii, jj, topology, vdWNeigh);
                        }
                    }

                } else if(iiHasCoul) {

                    cerr.printf("we should not have coul\n");

                    // self interaction of atoms in cell
                    for(j = i + 1; j < c1->nAtoms; j++) {

                        jj = c1->atomList[j];

                        // we check only for coul interactions
                        // because ii has no vdw anyway
                        jjHasCoul = hasCoul(jj, topology);

                        // but jj can have only have coul interactions
                        // to interact with ii
                        if(jjHasCoul)
                            addAtomJToNeighListOfI(ii, jj, topology, coulNeigh);
                    }

                    // interaction with atoms in neighCells
                    for (std::vector<int>::iterator it = c1->neighCellListId.begin(); it != c1->neighCellListId.end(); it++) {

                        const Cell *c2 = this->cells.data() + *it;

                        // skip empty cells
                        // this should be a rare case
                        // should happen for surface tension calculations
                        if(c2->nAtoms < 1)
                            continue;

                        for(j = 0; j < c2->nAtoms; j++) {

                            jj = c2->atomList[j];

                            // we check only for could interactions
                            // because ii has no vdw anyway
                            jjHasCoul = hasCoul(jj, topology);

                            // but jj can have only have coul interactions
                            // to interact with ii
                            if(jjHasCoul)
                                addAtomJToNeighListOfI(ii, jj, topology, coulNeigh);
                        }
                    }
                }
            }
        }

        addToTime(updateAtomNeighListOfLocalCellsPairlistTimer, updateAtomNeighListOfLocalCellsPairlistTime);
    }

    VHR_ALWAYS_INLINE
    void resetLocalCellsAndTheAtomNeighList()
    {
        startTimer(resetLocalCellsAndTheAtomNeighListPairlistTimer);

        size_t i, size = this->vdWNeigh.size();

        for (i = 0; i < size; i++) {
            this->vdWNeigh[i].nAtoms = 0;
            // this->vdWNeigh[i].listIds.clear();
        }

        size = this->coulNeigh.size();
        for (i = 0; i < size; i++) {
            this->coulNeigh[i].nAtoms = 0;
            // this->coulNeigh[i].listIds.clear();
        }

        size = this->vdWCoulNeigh.size();
        for (i = 0; i < size; i++) {
            this->vdWCoulNeigh[i].nAtoms = 0;
            // this->vdWCoulNeigh[i].listIds.clear();
        }

        size = this->cells.size();
        for(i = 0; i < size; i++) {
            this->cells[i].nAtoms = 0;
            // this->cells[i].atomList.clear();
        }

        addToTime(resetLocalCellsAndTheAtomNeighListPairlistTimer, resetLocalCellsAndTheAtomNeighListPairlistTime);
    }

    VHR_ALWAYS_INLINE
    void addAtomIToCellList(const int i, const int cellId)
    {
        // startTimer(addAtomIToCellListPairlistTimer);

        // Adding the atom to this cell
        size_t nAtoms = this->cells[cellId].nAtoms;

        // TODO
        // check the order of this if
        // the most frequent should come first
        if(nAtoms < this->cells[cellId].atomList.size()) {
            // cout.printf("firstIF\n");
            this->cells[cellId].atomList[nAtoms] = i;
            nAtoms++;
        }
        else {
            // cout.printf("secondIF\n");
            this->cells[cellId].atomList.push_back(i);
            nAtoms++;
        }

        this->cells[cellId].nAtoms = nAtoms;

        // addToTime(addAtomIToCellListPairlistTimer, addAtomIToCellListPairlistTime);
    }

    // TODO
    // Parallelize this routine. Because we removed domain decomposition, we lost the parallelization of
    // this routine. We should compute the cells in parallel the reduce the parallel cells
    VHR_INLINE
    bool placeAtomsIntoAllCells(const Configuration &conf, Box &simBox, Topology &topology)
    {
        startTimer(placeAtomsIntoAllCellsPairlistTimer);

        int i;
        const int nAtoms = conf.nAtoms;
        int cellId;

        IVec cellIJK;
        RVec boxBy2 = 0.5 * simBox.len;
        RVec delta = simBox.len / dd.nCells;


        // Placing the atoms into cells
        for (i = 0; i < nAtoms; i++) {

            cellIJK = floor((simBox.pbc(conf.current->x[i]) + boxBy2) / delta);
            // cerr.printf("atom %d cellIJK %d %d %d\n", i, cellIJK[XX], cellIJK[YY], cellIJK[ZZ]);
            // cerr.flush();

            cellId = dd.cellIJKToCellIndex(cellIJK[0], cellIJK[1], cellIJK[2]);

            addAtomIToCellList(i, cellId);
            // if(cellId > 63) {
            //     cerr.printf("atom %d at cellId %d\n", i, cellId);
            // }


        } // for (i = 0; i < conf.nAtoms; i++) {

        addToTime(placeAtomsIntoAllCellsPairlistTimer, placeAtomsIntoAllCellsPairlistTime);

        return true;
    }


    bool create(const Configuration &conf, Box &simBox, Topology &topology, const bool firstTime)
    {
        startTimer(pairListCreateTimer);

        int i;
        const int nCells = dd.nTotalCells();

        if(firstTime) {

            const int nAtoms = conf.nAtoms;
            // This is ready for MPI
            dd.getLocalCellsWithNeighInfo(this->cells);
            for(i = 0; i < nCells; i++)
                this->cells[i].nAtoms = 0;

            vdWNeigh.resize(nAtoms);
            coulNeigh.resize(nAtoms);
            vdWCoulNeigh.resize(nAtoms);

            if(this->doMaxDisplacement)
                refPosition = conf.current->x;

            // int nAtomsPerAtomList = estimateNAtomsPerCells(nAtoms, nCells);
            // nAtomsPerAtomList *= nAtomsPerAtomList;

            // for(i = 0; i < nAtoms; i++) {
            //     this->vdWNeigh[i].nAtoms = 0;
            //     this->vdWNeigh[i].listIds.resize(nAtomsPerAtomList);

            //     this->coulNeigh[i].nAtoms = 0;
            //     this->coulNeigh[i].listIds.resize(nAtomsPerAtomList);

            //     this->vdWCoulNeigh[i].nAtoms = 0;
            //     this->vdWCoulNeigh[i].listIds.resize(nAtomsPerAtomList);
            // }

        }

        placeAtomsIntoAllCells(conf, simBox, topology);

        // updateAtomNeighListOfLocalCells(topology);

        addToTime(pairListCreateTimer, pairListCreateTime);

        return true;
    }

    bool update(Configuration &conf, Box &simBox, Topology &topology)
    {
        startTimer(pairListUpdateTimer);

        int i;
        const int nCells = dd.nTotalCells();

        resetLocalCellsAndTheAtomNeighList();

        // cerr.printf("placeAtomsIntoAllCells\n");
        // cerr.flush();
        placeAtomsIntoAllCells(conf, simBox, topology);
        // cerr.printf("placeAtomsIntoAllCells\n");
        // cerr.flush();

        // cerr.printf("updateAtomNeighListOfLocalCells\n");
        // cerr.flush();
        updateAtomNeighListOfLocalCells(topology);
        // cerr.printf("updateAtomNeighListOfLocalCells\n");
        // cerr.flush();

        if(this->doMaxDisplacement)
            refPosition = conf.current->x;

        addToTime(pairListUpdateTimer, pairListUpdateTime);

        return true;
    }


    template<typename FileType>
    void print(FileType &outfile, int step)
    {
        long int total = 0;

        outfile.printf("        STEP: %d\n", step);
        outfile.printf("  localCells:\n");
        for(size_t i = 0; i < cells.size(); i++) {
            outfile.printf("cells[%d]:\n", i);
            cells[i].print(outfile);
            outfile.printf("----------\n\n");
        }
        outfile.printf("end of cells\n");

        outfile.printf("vdWNeigh:\n");
        for(size_t i = 0; i < vdWNeigh.size(); i++) {

            outfile.printf("check\n");
            outfile.printf("vdWNeigh[%d]: %d    %d\n", i, vdWNeigh[i].nAtoms,  vdWNeigh[i].atomId);

            std::sort(vdWNeigh[i].listIds.begin(), vdWNeigh[i].listIds.end());

            for(size_t j = 0; j < vdWNeigh[i].listIds.size(); j++) {
                outfile.printf("vdWNeigh[%d].listIds[%zu]: %d\n", i, j, vdWNeigh[i].listIds[j]);
                total++;
            }
        }
        outfile.printf("end of vdWNeigh: %d\n", total);

        // outfile.printf("coulNeigh:\n");
        // total = 0;
        // for(size_t i = 0; i < coulNeigh.size(); i++) {
        //     outfile.printf("coulNeigh[%d]: %d\n", i, coulNeigh[i].nAtoms);

        //     std::sort(coulNeigh[i].listIds.begin(), coulNeigh[i].listIds.end());

        //     for(size_t j = 0; j < coulNeigh[i].listIds.size(); j++) {
        //         outfile.printf("coulNeigh[%d].listIds[%zu]: %d\n", i, j, coulNeigh[i].listIds[j]);
        //         total++;
        //     }
        // }
        // outfile.printf("end of coulNeigh: %d\n", total);

        // outfile.printf("vdWCoulNeigh:\n");
        // total = 0;
        // for(size_t i = 0; i < vdWCoulNeigh.size(); i++) {
        //     outfile.printf("vdWCoulNeigh[%d]: %d\n", i, vdWCoulNeigh[i].nAtoms);

        //     std::sort(vdWCoulNeigh[i].listIds.begin(), vdWCoulNeigh[i].listIds.end());

        //     for(size_t j = 0; j < vdWCoulNeigh[i].listIds.size(); j++) {
        //         outfile.printf("vdWCoulNeigh[%d].listIds[%zu]: %d\n", i, j, vdWCoulNeigh[i].listIds[j]);
        //         total++;
        //     }
        // }
        // outfile.printf("end of vdWCoulNeigh: %d\n", total);

    }


    private:

    VHR_ALWAYS_INLINE
    bool hasVdW(const int ii, const Topology &topology)
    {
        const int type = topology.atomTypes[ii];

        // const int iac = ((type+1)*type)/2 + type;
        // TODO
        // Test this
        const int iac = ((type+3)*type)/2;
        const LJParameters lj = topology.ljParameters[iac];
        if(lj.c12 > 0 || lj.c6 > 0)
            return true;

        return false;
    }

    VHR_ALWAYS_INLINE
    bool hasCoul(const int ii, const Topology &topology)
    {
        const int charge = topology.charges[ii];

        // TODO
        // Test this
        if(charge != 0)
            return true;

        return false;
    }

};


#endif //PAIRLIST_PERATOM_DOMAIN_CPU_HPP
