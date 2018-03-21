#ifndef TRCFILEWRITER_HPP
#define TRCFILEWRITER_HPP

#include "../main.hpp"

#include "outputtraj.hpp"


class TRCFileWriter : public OutputTraj
{

public:

    TRCFileWriter();

    virtual ~TRCFileWriter();

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


        fprintf(file, "POSITIONRED\n");

        for(i = 0; i < nAtoms; i++) {
            fprintf(file, "%15.9f%15.9f%15.9f\n", conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ]);

            if ((i + 1) % 10 == 0)
                fprintf(file, "#%10d\n", i+1);
        }
        fprintf(file, "END\n");

        fprintf(file, "GENBOX\n");
        fprintf(file, "%5d\n%15.9f%15.9f%15.9f\n", 1, simBox.len[XX], simBox.len[YY], simBox.len[ZZ]);
        fprintf(file, "%15.9f%15.9f%15.9f\n", 90.0, 90.0, 90.0);
        fprintf(file, "%15.9f%15.9f%15.9f\n",  0.0,  0.0,  0.0);
        fprintf(file, "%15.9f%15.9f%15.9f\n",  0.0,  0.0,  0.0);
        fprintf(file, "END\n");

        addToTime(printTRJTimer, printTRJTime);

        return true;
    }

};

#endif //TRCFILEWRITER_HPP