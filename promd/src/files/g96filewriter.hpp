#ifndef G96FILEWRITER_HPP
#define G96FILEWRITER_HPP

#include "../main.hpp"

#include "outputtraj.hpp"

#include "../interaction/pairlist/pairlist.hpp"

#include <sstream>
#include <iomanip>

class G96FileWriter : public OutputTraj
{

public:

    G96FileWriter();

    virtual ~G96FileWriter();

    void writeTitle(const std::string title) {};

    // TODO implement the error check
    bool writeFile(const Configuration &conf, const Topology &topology, const Box &simBox, const number simTime, const int step)
    {
        startTimer(printTRJTimer);

        int i;
        const int nAtoms = conf.nAtoms;

        fprintf(file, "TITLE\n%sTIMESTEP\n%15d%15.9f\nEND\n", conf.title.c_str(), step, simTime);

        fprintf(file, "POSITION\n");
        fprintf(file, "# first 24 chars ignored\n");

        for(i = 0; i < nAtoms; i++) {
            int resId = topology.resIds[i];
            fprintf(file, "%5d %-5s %-6s%6d%15.9f%15.9f%15.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ]);
        }
        fprintf(file, "END\n");

        fprintf(file, "VELOCITY\n");
        fprintf(file, "# first 24 chars ignored\n");

        for(i = 0; i < nAtoms; i++) {
            int resId = topology.resIds[i];
            fprintf(file, "%5d %-5s %-6s%6d%15.9f%15.9f%15.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->v[i][XX], conf.current->v[i][YY], conf.current->v[i][ZZ]);
        }
        fprintf(file, "END\n");

        fprintf(file, "BOX\n");
        fprintf(file, "%15.9f%15.9f%15.9f\n", simBox.len[XX], simBox.len[YY], simBox.len[ZZ]);
        fprintf(file, "END\n");


        addToTime(printTRJTimer, printTRJTime);

        return true;
    }

    bool writeCells(const Configuration &conf, const Topology &topology, const Box &simBox, const PairList &pairlist, const number simTime, const int step)
    {
        startTimer(printTRJTimer);

        int cellId, i, j;

        for(cellId = 0; cellId < pairlist.cells.size(); cellId++) {

            const Cell *c1 = pairlist.cells.data() + cellId;
            const int nAtoms = c1->nAtoms;

            std::stringstream ss;
            ss << "cellId_" << std::setfill('0') << std::setw(5) << cellId << "_step_" << std::setfill('0') << std::setw(9) << step << ".g96";
            this->setFileName(ss.str());

            this->open();

            fprintf(file, "TITLE\n%s\ncell id %d\nnAtoms %d\nTIMESTEP\n %15d%15.9f\nEND\n", conf.title.c_str(), cellId, nAtoms, step, simTime);

            fprintf(file, "POSITION\n");
            fprintf(file, "# first 24 chars ignored\n");
            for (j = 0; j < nAtoms; j++) {

                i = c1->atomList[j];

                int resId = topology.resIds[i];
                fprintf(file, "%5d %-5s %-6s%6d%15.9f%15.9f%15.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->x[i][XX], conf.current->x[i][YY], conf.current->x[i][ZZ]);

            }
            fprintf(file, "END\n");

            fprintf(file, "VELOCITY\n");
            fprintf(file, "# first 24 chars ignored\n");
            for (j = 0; j < nAtoms; j++) {

                i = c1->atomList[j];

                int resId = topology.resIds[i];
                fprintf(file, "%5d %-5s %-6s%6d%15.9f%15.9f%15.9f\n", resId+1, topology.resNames[resId].c_str(), topology.atomNames[i].c_str(), i+1, conf.current->v[i][XX], conf.current->v[i][YY], conf.current->v[i][ZZ]);
            }
            fprintf(file, "END\n");


            fprintf(file, "BOX\n");
            fprintf(file, "%15.9f%15.9f%15.9f\n", simBox.len[XX]/dd.nTotalCellsX(), simBox.len[YY]/dd.nTotalCellsY(), simBox.len[ZZ]/dd.nTotalCellsZ());
            fprintf(file, "END\n");

            this->close();

        }

        addToTime(printTRJTimer, printTRJTime);

        return true;
    }

};

#endif //G96FILEWRITER_HPP