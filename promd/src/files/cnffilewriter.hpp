#ifndef CNFFILEWRITER_HPP
#define CNFFILEWRITER_HPP

#include "../main.hpp"

#include "outputtraj.hpp"


class CNFFileWriter : public OutputTraj
{

public:

    CNFFileWriter();

    virtual ~CNFFileWriter();

    void writeTitle(const std::string title)
    {
        fprintf(file, "TITLE\n%s\nEND\n", title.c_str());
    };

    // TODO implement the error check
    bool writeFile(const Configuration &conf, const Topology &topology, const Box &simBox, const number simTime, const int step)
    {
        startTimer(printTRJTimer);

        int i;
        const int nAtoms = conf.nAtoms;

        // fprintf(file, "TITLE\n%s\nEND\n", conf.title.c_str());
        fprintf(file, "TIMESTEP\n%15d%15.9f\nEND\n", step, simTime);


        fprintf(file, "POSITION\n");
        fprintf(file, "# first 24 chars ignored\n");

        for(i = 0; i < nAtoms; i++) {
            int resId = topology.resIds[i];
            fprintf(file, "%5d %-5s %-6s%6d%15.9f%15.9f%15.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ]);
        }
        fprintf(file, "END\n");

        // fprintf(file, "LATTICESHIFTS\n");
        // for(i = 0; i < nAtoms; i++) {
        //     fprintf(file, "%10d%10d%10d\n", 0, 0, 0);
        // }
        // fprintf(file, "END\n");

        fprintf(file, "VELOCITY\n");
        fprintf(file, "# first 24 chars ignored\n");

        for(i = 0; i < nAtoms; i++) {
            int resId = topology.resIds[i];
            fprintf(file, "%5d %-5s %-6s%6d%15.9f%15.9f%15.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->v[i][XX], conf.current->v[i][YY], conf.current->v[i][ZZ]);
        }
        fprintf(file, "END\n");

        fprintf(file, "GENBOX\n");
        fprintf(file, "%5d\n%15.9f%15.9f%15.9f\n", 1, simBox.len[XX], simBox.len[YY], simBox.len[ZZ]);
        fprintf(file, "%15.9f%15.9f%15.9f\n", 90.0, 90.0, 90.0);
        fprintf(file, "%15.9f%15.9f%15.9f\n",  0.0,  0.0,  0.0);
        fprintf(file, "%15.9f%15.9f%15.9f\n",  0.0,  0.0,  0.0);
        fprintf(file, "END\n");

        // fprintf(file, "BOX\n");
        // fprintf(file, "%15.9f%15.9f%15.9f\n", simBox.len[XX], simBox.len[YY], simBox.len[ZZ]);
        // fprintf(file, "END\n");


        addToTime(printTRJTimer, printTRJTime);

        return true;
    }

};

#endif //CNFFILEWRITER_HPP