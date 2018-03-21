#ifndef TRFFILEWRITER_HPP
#define TRFFILEWRITER_HPP

#include "../main.hpp"

#include "outputtraj.hpp"


class TRFFileWriter : public OutputTraj
{

public:

    TRFFileWriter();

    virtual ~TRFFileWriter();

    void writeTitle(const std::string title)
    {
        fprintf(file, "TITLE\n%sEND\n", title.c_str());
    };

    // TODO implement the error check
    bool writeFile(const Configuration &conf, const Topology &topology, const Box &simBox, const number simTime, const int step)
    {
        startTimer(printTRJTimer);

        int i;
        const int nAtoms = conf.nAtoms;

        fprintf(file, "TIMESTEP\n%15d%15.9f\nEND\n", step, simTime);


        fprintf(file, "FREEFORCERED\n");

        for(i = 0; i < nAtoms; i++) {
            fprintf(file, "%18.9f%18.9f%18.9f\n", conf.f[i][XX], conf.f[i][YY], conf.f[i][ZZ]);

            if ((i + 1) % 10 == 0)
                fprintf(file, "#%10d\n", i+1);
        }
        fprintf(file, "END\n");

        fprintf(file, "CONSFORCERED\n");

        for(i = 0; i < nAtoms; i++) {
            fprintf(file, "%18.9f%18.9f%18.9f\n", 0.0, 0.0, 0.0);

            if ((i + 1) % 10 == 0)
                fprintf(file, "#%10d\n", i+1);
        }
        fprintf(file, "END\n");

        addToTime(printTRJTimer, printTRJTime);

        return true;
    }

};

#endif //TRFFILEWRITER_HPP