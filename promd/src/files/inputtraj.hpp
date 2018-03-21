#ifndef INPUTTRAJ_HPP
#define INPUTTRAJ_HPP

#include "../main.hpp"
#include "../configuration/conf.hpp"
#include "../simulation/box.hpp"
#include "../simulation/simparam.hpp"
#include "../topology/topology.hpp"


#include "inputfile.hpp"

class InputTraj : public InputFile
{

public:

    InputTraj();

    virtual ~InputTraj();

    virtual bool readfile(Configuration &conf, Box &simBox, Topology &topology, SimParameters &simParam) = 0;

};

#endif //INPUTTRAJ_HPP