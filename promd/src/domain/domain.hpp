#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "../utilities/config.hpp"
#include "../utilities/timer.hpp"
#include "../utilities/math.hpp"

#include "../simulation/box.hpp"
#include "../interaction/pairlist/cell.hpp"

typedef bool MPI_Comm;
#define MPI_COMM_WORLD false

#include <algorithm>
#include <utility>
#include <unistd.h>

struct DomainDecomposition {

    int ierr;
    int rank, nprocs;
    size_t len;

    char hostname[200];

    IVec nCells;
    IVec shifts;

    // std::bitset<27> localCellsBitSet;

    // std::vector<int> neighIds;

    // std::vector<Cell> localCells;

    void init(int &argc, char* argv[])
    {
        nprocs=1;

        rank = 0;
        len = 200;
        gethostname(hostname, len);

        nCells[0] = 0;
        nCells[1] = 0;
        nCells[2] = 0;

        shifts[0] = 0;
        shifts[1] = 0;
        shifts[2] = 0;
    }

    inline
    bool master() const
    {
        // return rank == 0;
        return true;
    }


    inline
    int nEstimatedGhostCells() const
    {
        return (shifts[0] * nCells[0]) * (shifts[1] * nCells[1]) * (shifts[2] * nCells[2]);
    }

    inline
    int nTotalCells() const
    {
        return (nCells[0] * nCells[1] * nCells[2]);
    }

    inline
    IVec totalCells() const
    {
        return nCells;
    }

    inline
    int nTotalCellsX() const
    {
        return nCells[0];
    }

    inline
    int nTotalCellsY() const
    {
        return nCells[1];
    }

    inline
    int nTotalCellsZ() const
    {
        return nCells[2];
    }

    void abort()
    {
        exit(1);
    }

    void barrier() const
    {
    }

    void finalize()
    {
        exit(0);
    }


    void setNumberProcs(const int nthreads = 0)
    {
        if(nthreads != 0) {
            this->nprocs = nthreads;
        }
        else {
            #if defined(_OPENMP)
            this->nprocs = omp_get_max_threads();
            #else
            this->nprocs = 1;
            #endif
        }
    }

    // inline
    bool splitBoxIntoCells(const number natoms, const number &cutoff, const Box &simBox)
    {
        // std::vector<int> primes, histogramOfPrimes;
        // const int ndiv = factorizeIntoPrimes(nprocs, primes, histogramOfPrimes);

        // nCells[0] = 1;
        // nCells[1] = 1;
        // nCells[2] = 1;

        // int possibleDecomposition[3];
        // possibleDecomposition[0] = 1;
        // possibleDecomposition[1] = 1;
        // possibleDecomposition[2] = 1;

        // splitBoxBasedOnNumberOfProcessors(simBox, cutoff, natoms, ndiv, primes.data(), histogramOfPrimes.data(), possibleDecomposition);

        // if(nprocs > 1 && (nCells[0] == 1  && nCells[1] == 1 && nCells[2] == 1)) {
        //     return false;
        // }

        // fprintf(stderr, "possibleDecomposition: %dx%dx%d\n", possibleDecomposition[0], possibleDecomposition[1], possibleDecomposition[2]);

        // const int nCutOffCellsX = floor(cellrat * simBox.x / cutoff);
        // const int nCutOffCellsY = floor(cellrat * simBox.y / cutoff);
        // const int nCutOffCellsZ = floor(cellrat * simBox.z / cutoff);

        // if(master()) {
        //     fprintf(stderr, "nCellsIdeal: %dx%dx%d\n", nCutOffCellsX, nCutOffCellsY, nCutOffCellsZ);
        //     fprintf(stderr, "     nCells: %dx%dx%d\n", nCells[0], nCells[1], nCells[2]);
        // }

        splitBoxesIntoCells(natoms, cutoff, simBox);

        int i;
        for(i = 0; i < DIM; i++) {
            if(nCells[i] == 1)
                shifts[i] = 0;
            else
                shifts[i] = ceil(nCells[i] * cutoff / simBox.len[i]);
        }

        // if(nTotalCellsX() == 1)
        //     shifts[0] = 0;
        // else
        //     shifts[0] = ceil(nTotalCellsX() * cutoff / simBox.len[XX]);

        // if(nTotalCellsY() == 1)
        //     shifts[1] = 0;
        // else
        //     shifts[1] = ceil(nTotalCellsY() * cutoff / simBox.x);

        // if(nTotalCellsZ() == 1)
        //     shifts[2] = 0;
        // else
        //     shifts[2] = ceil(nTotalCellsZ() * cutoff / simBox.x);

        return true;
    }

    inline
    void getGhostCells(const int cellId, std::vector<int> &ghostCells)
    {
        int i, j, k, adj;
        ghostCells.resize(0);
        ghostCells.reserve(nEstimatedGhostCells());

        for(i = -shifts[0]; i <= shifts[0]; i++) {
            for(j = -shifts[1]; j <= shifts[1]; j++) {
                for(k = -shifts[2]; k <= shifts[2]; k++) {

                    if(i == 0 && j == 0 && k == 0) {
                        continue;
                    }

                    adj = getAdjCellId(cellId, i, j, k);

                    if(adj > cellId) {
                        ghostCells.push_back(adj);
                    }

                }
            }
        }

        // Sorting list. This is helpfull in cache coerency
        std::sort(ghostCells.begin(), ghostCells.end());

        // // This should not necessary but I guess that it is safe to add this.
        // // Because it is sorted it is "easier" but it is still "costily" to remove any repeated element! :-)
        // neighIds.erase(std::unique(neighIds.begin(), neighIds.end()),neighIds.end());
    }

    inline
    int cellIJKToCellIndex(const int i, const int j, const int k) const
    {
        // Sollution taken from:
        // http://stackoverflow.com/questions/7367770/how-to-flatten-or-index-3d-array-in-1d-array
        return (k * nTotalCellsX() * nTotalCellsY()) + (j * nTotalCellsX()) + i;
    }

    inline
    void getLocalCellsWithNeighInfo(std::vector<Cell> &localCells)
    {
        int i, j, k, l;
        int lowerI, lowerJ, lowerK;


        // const int lowerCellIndex = lowerstCellIndexFromMPIRank(rank);
        // This is already prepared for MPI
        const int lowerCellIndex = 0;
        cellIndexToIJK(lowerCellIndex, lowerI, lowerJ, lowerK);

        const int higherI = lowerI + nCells[0];
        const int higherJ = lowerJ + nCells[1];
        const int higherK = lowerK + nCells[2];

        // localCells.resize(nLocalCells());
        localCells.resize(nTotalCells());

        // Doing in the reverse order for we get them sorted
        int counter = 0;
        for(k = lowerJ; k < higherK; k++) {
            for(j = lowerJ; j < higherJ; j++) {
                for(i = lowerI; i < higherI; i++) {

                    const int cellId = cellIJKToCellIndex(i, j, k);
                    localCells[counter].id = cellId;

                    // Getting neighbor info
                    cellNeigh(cellId, localCells[counter].neighCellListId);
                    int size = localCells[counter].neighCellListId.size();
                    localCells[counter].nNeighCells = size;

                    // Getting the mpi ranks info
                    // localCells[counter].mpiRanksToCommunicate.resize(0);
                    // localCells[counter].mpiRanksToCommunicate.reserve(nprocs);

                    // for(l = 0; l < size; l++) {
                    //     int mpiRank = getMPIRankThatHoldsCell(cellId);

                    //     if(mpiRank != rank)
                    //         localCells[counter].mpiRanksToCommunicate.push_back(mpiRank);
                    // }

                    // // Sorting list. This is helpfull in cache coerency
                    // std::sort(localCells[counter].mpiRanksToCommunicate.begin(), localCells[counter].mpiRanksToCommunicate.end());

                    // // Removing repetitions of mpiRanks
                    // localCells[counter].mpiRanksToCommunicate.erase(std::unique(localCells[counter].mpiRanksToCommunicate.begin(), localCells[counter].mpiRanksToCommunicate.end()), localCells[counter].mpiRanksToCommunicate.end());


                    // Sorting list. This is helpfull in cache coerency
                    std::sort(localCells[counter].neighCellListId.begin(), localCells[counter].neighCellListId.end());

                    // Removing repetitions of neighCells
                    localCells[counter].neighCellListId.erase(std::unique(localCells[counter].neighCellListId.begin(), localCells[counter].neighCellListId.end()), localCells[counter].neighCellListId.end());

                    counter++;
                }
            }
        }
    }

private:

    inline
    void cellNeigh(const int cellId, std::vector<int> &neighIds)
    {
        int i, j, k;
        neighIds.resize(0);
        neighIds.reserve(nEstimatedGhostCells());

        for(i = -shifts[0]; i <= shifts[0]; i++) {
            for(j = -shifts[1]; j <= shifts[1]; j++) {
                for(k = -shifts[2]; k <= shifts[2]; k++) {

                    if(i == 0 && j == 0 && k == 0)
                        continue;

                    const int adj = getAdjCellId(cellId, i, j, k);

                    if(adj > cellId)
                        // neighIds.push_back(getAdjCellId(cellId, i, j, k));
                        neighIds.push_back(adj);



                }
            }
        }

        // Sorting list. This is helpfull in cache coerency
        std::sort(neighIds.begin(), neighIds.end());

        // // This should not necessary but I guess that it is safe to add this.
        // // Because it is sorted it is "easier" but it is still "costily" to remove any repeated element! :-)
        // neighIds.erase(std::unique(neighIds.begin(), neighIds.end()),neighIds.end());
    }


    inline
    void getLocalCells(std::vector<int> &localCells)
    {
        // const int tid = rank + 1;
        int i, j, k;
        int lowerI, lowerJ, lowerK, higherI, higherJ, higherK;

        // const int lowerCellIndex = lowerstCellIndexFromMPIRank(rank);
        // This is already prepared for MPI
        const int lowerCellIndex = 0;

        cellIndexToIJK(lowerCellIndex, lowerI, lowerJ, lowerK);
        higherI = lowerI + nCells[0];
        higherJ = lowerJ + nCells[1];
        higherK = lowerK + nCells[2];

        localCells.resize(0);
        localCells.reserve((higherI - lowerI) * (higherJ - lowerJ) * (higherK - lowerK));

        // Doing in the reverse order for we get them sorted
        for(k = lowerJ; k < higherK; k++) {
            for(j = lowerJ; j < higherJ; j++) {
                for(i = lowerI; i < higherI; i++) {

                    localCells.push_back(cellIJKToCellIndex(i, j, k));

                }
            }
        }
    }

    inline
    void getGhostCells(std::vector<int> &localCells, std::vector< std::pair <int, int> > &ghostCellList)
    // void getGhostCellList(std::vector<int> &localCells, std::vector<int> &ghostCellList, std::vector<int> &mpiRankList)
    {
        size_t i, j;
        std::vector<int> neighIds(nTotalCells(), 0);

        for(i = 0; i < localCells.size(); i++) {

            cellNeigh(localCells[i], neighIds);

            for(j = 0; j < neighIds.size(); j++) {
                std::pair<int, int> foo = std::make_pair(getMPIRankThatHoldsCell(neighIds[j]), neighIds[j]);
                ghostCellList.push_back(foo);

                fprintf(stderr, "rank: %5d, localCells[%5zu]: %5d, neighIds[%5zu]: %5d, getMPIRankThatHoldsCell: %d\n", rank, i, localCells[i], j, neighIds[j], foo.first);
            }
        }

        // Sorting list. This is helpfull in cache coerency
        std::sort(ghostCellList.begin(), ghostCellList.end());

        ghostCellList.erase(std::unique(ghostCellList.begin(), ghostCellList.end()),ghostCellList.end());
    }

    inline
    void getGhostCells(std::vector<int> &localCells, std::vector<int> &ghostCellList)
    {
        size_t i, j;
        std::vector<int> neighIds(nTotalCells(), 0);

        for(i = 0; i < localCells.size(); i++) {

            cellNeigh(localCells[i], neighIds);

            for(j = 0; j < neighIds.size(); j++) {

                if(getMPIRankThatHoldsCell(neighIds[j]) != rank)
                    ghostCellList.push_back(neighIds[j]);
            }
        }

        // Sorting list. This is helpfull in cache coerency
        std::sort(ghostCellList.begin(), ghostCellList.end());

        ghostCellList.erase(std::unique(ghostCellList.begin(), ghostCellList.end()),ghostCellList.end());
    }

    inline
    void getInteractingCells(std::vector<int> &localCells, std::vector<int> &ghostCellList)
    {
        size_t i, j;
        std::vector<int> neighIds(nTotalCells(), 0);

        for(i = 0; i < localCells.size(); i++) {

            cellNeigh(localCells[i], neighIds);

            for(j = 0; j < neighIds.size(); j++) {
                ghostCellList.push_back(neighIds[j]);
            }
        }

        // Sorting list. This is helpfull in cache coerency
        std::sort(ghostCellList.begin(), ghostCellList.end());

        ghostCellList.erase(std::unique(ghostCellList.begin(), ghostCellList.end()),ghostCellList.end());
    }

    //!
    // Trial Division Algorithm
    // This function implements the simple and "expensive" trial division algorithm
    // in order to factorize numbers into primes.
    //
    // TODO
    // Probably. I am not sure. Because computers and people tend to use number of
    // cores and mpi process that are multiples of 2, 3, 5 and sometimes 7.
    // Instead of linearly divide the numbers by a linear sequence of divisors
    // the divisor could be taken from a prime number generator. More info at:
    // https://en.wikipedia.org/wiki/Sieve_of_Atkin
    //
    // Algorithm taken from:
    // https://en.wikipedia.org/wiki/Trial_division
    //
    int factorizeIntoPrimes(int n, std::vector<int> &primes, std::vector<int> &histogramOfPrimes)
    {
        int divisor, ndivisors;

        // The resize part should not be necessary for it is expected that:
        // 1. primes and histogramOfPrimes are local variables;
        // 2. this function is called only once in the program; and
        // 3. the next programmer will read the function body before using it.
        primes.resize(0);
        histogramOfPrimes.resize(0);

        // Because we have the histogramOfPrimes we do not need to store
        // a lot of memory. If we take size 5 the number has to be multiple of
        // 2 * 3 * 5 * 7 * 11 = 2310. Which is kinda insane, right?
        primes.reserve(5);
        histogramOfPrimes.reserve(5);

        divisor = 2;
        ndivisors = 0;
        while (n > 1) {
            while (n % divisor == 0) {
                if (ndivisors == 0 || primes[ndivisors-1] != divisor) {
                    ndivisors++;
                    primes.push_back(divisor);
                    histogramOfPrimes.push_back(1);
                }
                else
                    histogramOfPrimes[ndivisors-1]++;
                n /= divisor;
            }
            divisor++;
        }
        return ndivisors;
    }

    //!
    // MPI communication cost
    //
    // This functions computes the raw communication cost for a given possible box
    // decomposition. The cost does not take into account the particle density nor
    // the vector size.
    // The units are arbitrary.
    //
    // This function assumes that we will compute the nonbonded calculations using
    // Newton's thrid law.
    //
    // TODO
    // improve this function to be more realistic.
    //
    number costMPICommunication(const int natoms, const number cutoff, const Box &simBox, IVec &possibleDecomposition)
    {
        int  i, j, k;
        number commCost;

        const number M_PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

        const number numberOfMPICalls[3] = { possibleDecomposition[0] * cellrat * cutoff / simBox.len[XX],
                                             possibleDecomposition[1] * cellrat * cutoff / simBox.len[YY],
                                             possibleDecomposition[2] * cellrat * cutoff / simBox.len[ZZ] };

        const int nCellsx = floor(cellrat * simBox.len[XX] / cutoff);
        const int nCellsy = floor(cellrat * simBox.len[YY] / cutoff);
        const int nCellsz = floor(cellrat * simBox.len[ZZ] / cutoff);

        commCost = 0;

        if(nCellsx < possibleDecomposition[0])
            commCost = (float)possibleDecomposition[0]/(float)nCellsx;
        else if(nCellsy < possibleDecomposition[1])
            commCost = (float)possibleDecomposition[1]/(float)nCellsy;
        else if(nCellsz < possibleDecomposition[2])
            commCost = (float)possibleDecomposition[2]/(float)nCellsz;

        for(i = 0; i < 3; i++) {
            if (possibleDecomposition[i] > 1) {
                commCost += numberOfMPICalls[i];
                for(j = i + 1; j < 3; j++) {
                    if (possibleDecomposition[j] > 1) {
                        commCost += numberOfMPICalls[i]*numberOfMPICalls[j]*M_PI/4;
                        for(k = j + 1; k < 3; k++) {
                            if (possibleDecomposition[k] > 1) {
                                commCost += numberOfMPICalls[i]*numberOfMPICalls[j]*numberOfMPICalls[k]*M_PI/6;
                            }
                        }
                    }
                }
            }
        }

        // fprintf(stderr, "commCost: %f, pd[0]: %d, pd[1]: %d, pd[2]: %d\n", commCost, possibleDecomposition[0], possibleDecomposition[1], possibleDecomposition[2]);

        return commCost;
    }

    //!
    //
    //
    // TODO
    // Implement a non-recursive version of this function
    //
    //
    void splitBoxBasedOnNumberOfProcessors(const Box &simBox, const number cutoff, const int natoms, const int ndiv, int *primes, int *histogramOfPrimes, IVec &possibleDecomposition)
    {
        int x, y, i;
        float ce, co;

        if (ndiv == 0) {
            co = costMPICommunication(natoms, cutoff, simBox, nCells);
            ce = costMPICommunication(natoms, cutoff, simBox, possibleDecomposition);
            if(ce >= 0 && ((nCells[0] == 1  && nCells[1] == 1 && nCells[2] == 1) || ce < co)) {
                nCells[0] = possibleDecomposition[0];
                nCells[1] = possibleDecomposition[1];
                nCells[2] = possibleDecomposition[2];
            }
            return;
        }

        for(x = histogramOfPrimes[0]; x >= 0; x--) {

            //
            for(i = 0; i < x; i++)
                possibleDecomposition[0] *= primes[0];

            for(y = histogramOfPrimes[0] - x; y >= 0; y--) {

                //
                for(i = 0; i < y; i++)
                    possibleDecomposition[1] *= primes[0];

                for(i = 0; i < histogramOfPrimes[0] - x - y; i++)
                    possibleDecomposition[2] *= primes[0];

                splitBoxBasedOnNumberOfProcessors(simBox, cutoff, natoms, ndiv-1, primes+1, histogramOfPrimes+1, possibleDecomposition);

                for(i = 0; i < histogramOfPrimes[0] - x - y; i++)
                    possibleDecomposition[2] /= primes[0];

                for(i = 0; i < y; i++)
                    possibleDecomposition[1] /= primes[0];
            }

            //
            for(i = 0; i < x; i++)
                possibleDecomposition[0] /= primes[0];

        }
    }

    inline
    int lowerstCellIndexFromMPIRank(int mpiRank)
    {
        // Sollution taken from:
        // http://stackoverflow.com/questions/7367770/how-to-flatten-or-index-3d-array-in-1d-array
        const int nCellsx = nCells[0];
        const int nCellsy = nCells[1];

        int lowerI, lowerJ, lowerK;

        lowerK = mpiRank / (nCellsx * nCellsy);

        mpiRank -= (lowerK * nCellsx * nCellsy);
        lowerJ = mpiRank / nCellsx;
        lowerI = mpiRank % nCellsx;

        // Mapping the mpi cells into global cells
        lowerI *= nCells[0];
        lowerJ *= nCells[1];
        lowerK *= nCells[2];

        return cellIJKToCellIndex(lowerI, lowerJ, lowerK);
    }

    inline
    int getMPIRankThatHoldsCell(const int cellId)
    {
        // Sollution taken from:
        // http://stackoverflow.com/questions/7367770/how-to-flatten-or-index-3d-array-in-1d-array
        int i, j, k;
        cellIndexToIJK(cellId, i, j, k);

        // Mapping the mpi cells into global cells
        i /= nCells[0];
        j /= nCells[1];
        k /= nCells[2];

        return (k * nCells[0] * nCells[1]) + (j * nCells[0]) + i;
    }


    inline
    bool cellIndexToIJK(int cellId, int &i, int &j, int &k)
    {
        // Sollution taken from:
        // http://stackoverflow.com/questions/7367770/how-to-flatten-or-index-3d-array-in-1d-array

        const int nCellsx = nTotalCellsX();
        const int nCellsy = nTotalCellsY();

        k = cellId / (nCellsx * nCellsy);

        cellId -= (k * nCellsx * nCellsy);
        j = cellId / nCellsx;
        i = cellId % nCellsx;

        return true;
    }



    inline
    int pbcIndex(int i, const int dimSize)
    {
        while(i >= dimSize) i -= dimSize;
        while(i < 0) i += dimSize;

        return i;
    }

    inline
    int getAdjCellId(const int cellId, const int shiftI, const int shiftJ, const int shiftK)
    {
        // const int nTotalCells = (nCells[0] * nCells[0]) * (nCells[1] * nCells[1]) * (nCells[2] * nCells[2]);

        int ci, cj, ck, c0i, c0j, c0k;

        cellIndexToIJK(cellId, c0i, c0j, c0k);

        ci = pbcIndex(c0i + shiftI, nCells[0]);
        cj = pbcIndex(c0j + shiftJ, nCells[1]);
        ck = pbcIndex(c0k + shiftK, nCells[2]);

        // const int cid = cellIJKToCellIndex(ci, cj, ck);
        // return cid < nTotalCells ? cid : -1;

        return cellIJKToCellIndex(ci, cj, ck);;
    }

    inline
    void splitBoxesIntoCells(const number natoms, const number &cutoff, const Box &simBox)
    {
        // nCells[0] = floor(cellrat * simBox.x / cutoff / nCells[0]);
        // nCells[1] = floor(cellrat * simBox.y / cutoff / nCells[1]);
        // nCells[2] = floor(cellrat * simBox.z / cutoff / nCells[2]);

        nCells[0] = floor(cellrat * simBox.len[XX] / cutoff);
        nCells[1] = floor(cellrat * simBox.len[YY] / cutoff);
        nCells[2] = floor(cellrat * simBox.len[ZZ] / cutoff);

        if(nCells[0] == 0)
            nCells[0] = 1;

        if(nCells[1] == 0)
            nCells[1] = 1;

        if(nCells[2] == 0)
            nCells[2] = 1;
    }



};

extern DomainDecomposition dd;

#endif //DOMAIN_HPP