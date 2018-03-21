#ifndef __TIMER_H__
#define __TIMER_H__

#include <sys/time.h>
#include <cstdio>

#include "config.hpp"

// #include <ctime>

#if TIMEMS
#define timeUnit "ms"
#define nsdayconstant 86400.0
#else
#define timeUnit "s"
#define nsdayconstant 86.4
#endif


// Timing variables
extern
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


extern
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
extern
timeval addAtomJToNeighListOfIPairlistTimer,
        updateAtomNeighListOfLocalCellsPairlistTimer,
        resetLocalCellsAndTheAtomNeighListPairlistTimer,
        addAtomIToCellListPairlistTimer,
        placeAtomsIntoAllCellsPairlistTimer;

extern
double  updateAtomNeighListOfLocalCellsPairlistTime,
        resetLocalCellsAndTheAtomNeighListPairlistTime,
        placeAtomsIntoAllCellsPairlistTime;


// bonded
extern
timeval bondTimer;

extern
double  bondTime;

static
VHR_ALWAYS_INLINE
void startTimer(timeval &atimeval)
{
    // time(&atimeval);
    gettimeofday(&atimeval, NULL);
}

// time elapsed in ms
static
VHR_ALWAYS_INLINE
double getElapsed(timeval &atimeval)
{
    // timeval timerNow;
    // time(&timerNow);
    // return difftime(timerNow, atimeval);

    struct timeval timerNow, timerElapsed;

    gettimeofday(&timerNow, NULL);
    timersub(&timerNow, &atimeval, &timerElapsed);

#if TIMEMS
    return timerElapsed.tv_sec*1000.0+timerElapsed.tv_usec/1000.0;
#else
    return timerElapsed.tv_sec+timerElapsed.tv_usec/1000000.0;
#endif


}

static
VHR_ALWAYS_INLINE
void addToTime(timeval &atimeval, double &time)
{
    time += getElapsed(atimeval);
}

#endif // __TIMER_H__