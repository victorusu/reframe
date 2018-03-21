#ifndef GROFILEWRITER_HPP
#define GROFILEWRITER_HPP

#include "../main.hpp"

#include "outputtraj.hpp"

class GROFileWriter : public OutputTraj
{

public:

    GROFileWriter();

    virtual ~GROFileWriter();

    void writeTitle(const std::string title) {};

    // TODO implement the error check
    bool writeFile(const Configuration &conf, const Topology &topology, const Box &simBox, const number simTime, const int step)
    {
        startTimer(printTRJTimer);

        int i;
        const int nAtoms = conf.nAtoms;

        std::string title = conf.title;
        std::replace(title.begin(), title.end(), '\n', ' ');
        std::replace(title.begin(), title.end(), '\t', ' ');

        fprintf(file, "%s timestep: %d   time: %f\n", title.c_str(), step, simTime);
        fprintf(file, "%10d\n", nAtoms);

        for(i = 0; i < nAtoms; i++) {
            int resId = topology.resIds[i];
            fprintf(file, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ], conf.current->v[i][XX], conf.current->v[i][YY], conf.current->v[i][ZZ]);
            // fprintf(file, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f   %22.9f%22.9f%22.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ], conf.current->v[i][XX], conf.current->v[i][YY], conf.current->v[i][ZZ], conf.f[i][XX], conf.f[i][YY], conf.f[i][ZZ]);
            // fprintf(file, "%5d%-5s%5s%5d%15.9f%15.9f%15.9f%22.9f%22.9f%22.9f   %22.9f%22.9f%22.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ], conf.current->v[i][XX], conf.current->v[i][YY], conf.current->v[i][ZZ], conf.f[i][XX], conf.f[i][YY], conf.f[i][ZZ]);
        }
            // fprintf(file, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n", i+1, topology.resNames[i].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ],  conf.current->v[i][XX], conf.current->v[i][YY], conf.current->v[i][ZZ]);


        fprintf(file, " %9.5f%9.5f%9.5f\n", simBox.len[XX], simBox.len[YY], simBox.len[ZZ]);

        addToTime(printTRJTimer, printTRJTime);

        return true;
    }

};

#endif //GROFILEWRITER_HPP