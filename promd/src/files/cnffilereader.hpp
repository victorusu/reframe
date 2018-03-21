#ifndef CNFFILEREADER_HPP
#define CNFFILEREADER_HPP

#include "../main.hpp"
#include "../configuration/conf.hpp"
#include "../simulation/box.hpp"
#include "../topology/topology.hpp"

#include "blocksfile.hpp"

#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <cstdlib>
#include <exception>


class CNFFileReader : public BlocksFile
{

public:

    CNFFileReader();

    virtual ~CNFFileReader();

    bool readfile(Configuration &conf, Box &simBox, Topology &topology, SimParameters &simParam)
    {
        number x, y, z;

        std::string line;

        // reading the title
        if(!findBlock("TITLE")) {
            fprintf(stderr, "TITLE not found in file %s\n", getFileName().c_str());
            return false;
        }

        conf.title = "";
        while(readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
            conf.title.append(line + '\n');
        }

        // This looks for END not inside the title but as a block name END
        // if(!findBlock("END")) {
        //     fprintf(stderr, "TITLE block incomplete in file %s\n", getFileName().c_str());
        //     return false;
        // }

        if(!findBlock("POSITION")) {
            fprintf(stderr, "POSITION not found in file %s\n", getFileName().c_str());
            return false;
        }

        for(int i = 0; i < conf.nAtoms; i++) {
            line = readLineSkippingCommentedLines('#');

            try {
                x = atof(line.substr(23,15).c_str());
                y = atof(line.substr(38,15).c_str());
                z = atof(line.substr(53,15).c_str());
            } catch(std::exception e) {
                return false;
            }

            pushBack(conf.current->x, x, y, z);

            // conf.current->y.push_back(y);
            // conf.current->z.push_back(z);
        }
        line = readLineSkippingCommentedLines('#');
        if(line.find("END") != 0) {
            fprintf(stderr, "POSITION block incomplete in file %s\n", getFileName().c_str());
            return false;
        }

        if(findBlock("VELOCITY")) {
            for(int i = 0; i < conf.nAtoms; i++) {
                line = readLineSkippingCommentedLines('#');

                try {
                    x = atof(line.substr(23,15).c_str());
                    y = atof(line.substr(38,15).c_str());
                    z = atof(line.substr(53,15).c_str());
                } catch(std::exception e) {
                    return false;
                }

                pushBack(conf.current->v, x, y, z);

                // conf.current->vx.push_back(x);
                // conf.current->vy.push_back(y);
                // conf.current->vz.push_back(z);
                // conf.vxOld.push_back(x);
                // conf.vyOld.push_back(y);
                // conf.vzOld.push_back(z);
            }
            line = readLineSkippingCommentedLines('#');
            if(line.find("END") != 0) {
                fprintf(stderr, "VELOCITY block incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        } else {
            for(int i = 0; i < conf.nAtoms; i++) {
                pushBack(conf.current->v, 0.0, 0.0, 0.0);
                // conf.current->vx.push_back(0.0);
                // conf.current->vy.push_back(0.0);
                // conf.current->vz.push_back(0.0);
                // conf.vxOld.push_back(0.0);
                // conf.vyOld.push_back(0.0);
                // conf.vzOld.push_back(0.0);
            }
        }

        if(!findBlock("GENBOX")) {
            fprintf(stderr, "GENBOX not found in file %s\n", getFileName().c_str());
            return false;
        }

        // read box type
        line = readLineSkippingCommentedLines('#');
        // read box first line
        line = readLineSkippingCommentedLines('#');

        try {
            std::stringstream str(line);
            // x = atof(line.substr(0,15).c_str());
            // y = atof(line.substr(15,15).c_str());
            // z = atof(line.substr(30,15).c_str());
            str >> x;
            str >> y;
            str >> z;

        } catch(std::exception e) {
            return false;
        }

        simBox.len[XX] = x;
        simBox.len[YY] = y;
        simBox.len[ZZ] = z;

        // placeMolsIntoBox(conf, simBox);

        return true;
    }

};

#endif //CNFFILEREADER_HPP