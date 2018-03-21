#ifndef OUTPUTTRAJ_HPP
#define OUTPUTTRAJ_HPP

#include "../main.hpp"
#include "../configuration/conf.hpp"
#include "../topology/topology.hpp"
#include "../simulation/box.hpp"

#include "outputfile.hpp"

class OutputTraj : public OutputFile
{

public:

    OutputTraj();

    virtual ~OutputTraj();

    virtual void writeTitle(const std::string title) = 0;

    virtual bool writeFile(const Configuration &conf, const Topology &topology, const Box &simBox, const number time, const int step) = 0;

    VHR_ALWAYS_INLINE
    void flush()
    {
    	fflush(file);
    }

};

#endif //OUTPUTTRAJ_HPP