#include "timer.hpp"

// struct timeval;

// Timing variables
timeval totalProgramTimer,
        mdLoopTimer,
        pairListCreateTimer,
        pairListUpdateTimer,
        integrationTimer,
        kineticEnergyTimer,
        nonbondedComputeTimer,
        nonbondedFirstTimer,
        nonbondedSecondTimer,
        nonbondedReductionTimer,
        printTRJTimer;

double  totalProgramTime,
        mdLoopTime,
        pairListCreateTime,
        pairListUpdateTime,
        integrationTime,
        kineticEnergyTime,
        nonbondedComputeTime,
        nonbondedFirstTime,
        nonbondedSecondTime,
        nonbondedReductionTime,
        printTRJTime;


// pairlist
timeval addAtomJToNeighListOfIPairlistTimer,
        updateAtomNeighListOfLocalCellsPairlistTimer,
        resetLocalCellsAndTheAtomNeighListPairlistTimer,
        addAtomIToCellListPairlistTimer,
        placeAtomsIntoAllCellsPairlistTimer;

double  updateAtomNeighListOfLocalCellsPairlistTime,
        resetLocalCellsAndTheAtomNeighListPairlistTime,
        placeAtomsIntoAllCellsPairlistTime;

// bonded
timeval bondTimer;

double  bondTime;


