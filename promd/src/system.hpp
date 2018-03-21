#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "main.hpp"

#include "imdfile.hpp"
#include "omdfile.hpp"
#include "conf.hpp"
#include "cnffile.hpp"
#include "pairlist.hpp"
#include "nonbonded.hpp"
#include "integrator.hpp"

#include "timer.hpp"

struct Replica
{

    SimParameters simParam;
    Configuration conf;
    Box simBox;
    PairList pairlist;
    Nonbonded nonbonded;
    Integrator integrator;



};


#endif //SYSTEM_HPP