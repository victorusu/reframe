#ifndef GROFILEREADER_HPP
#define GROFILEREADER_HPP

#include "../main.hpp"
#include "../configuration/conf.hpp"
#include "../simulation/box.hpp"
#include "../topology/topology.hpp"

#include "inputtraj.hpp"

#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>

class GROFileReader : public InputTraj
{

public:

    GROFileReader();

    virtual ~GROFileReader();

    bool readfile(Configuration &conf, Box &simBox, Topology &topology, SimParameters &simParam)
    {
        int n;
        number x, y, z, vx, vy, vz, fx, fy, fz;

        std::string line;

        // reading the title
        // getline(cnffile, line);
        line = readLine();
        conf.title = line;

        // reading the number of atoms in file
        // getline(cnffile, line);
        line = readLine();
        n = atoi(line.c_str());

        if(n != conf.nAtoms)
            return false;

        for(int i = 0; i < conf.nAtoms; i++) {
            line = readLine();
            // getline(cnffile, line);

            try {
                x = atof(line.substr(20,8).c_str());
                y = atof(line.substr(28,8).c_str());
                z = atof(line.substr(36,8).c_str());
            } catch(std::exception e) {
                return false;
            }

            try {
                vx = atof(line.substr(44,8).c_str());
                vy = atof(line.substr(52,8).c_str());
                vz = atof(line.substr(60,8).c_str());
            } catch(std::exception e) {
                vx = 0.0;
                vy = 0.0;
                vz = 0.0;
            }

            try {
                fx = atof(line.substr(64,8).c_str());
                fy = atof(line.substr(72,8).c_str());
                fz = atof(line.substr(80,8).c_str());
            } catch(std::exception e) {
                fx = 0.0;
                fy = 0.0;
                fz = 0.0;
            }

            pushBack(conf.current->x, x, y, z);

            pushBack(conf.current->v, vx, vy, vz);

            pushBack(conf.f, fx, fy, fz);

            // conf.current->x.push_back(x);
            // conf.current->y.push_back(y);
            // conf.current->z.push_back(z);

            // conf.current->vx.push_back(vx);
            // conf.current->vy.push_back(vy);
            // conf.current->vz.push_back(vz);
            // conf.vxcurrent.push_back(x);
            // conf.vycurrent.push_back(y);
            // conf.vzcurrent.push_back(z);

            // conf.fx.push_back(fx);
            // conf.fy.push_back(fy);
            // conf.fz.push_back(fz);
        }

        // reading the box
        // TODO correct this next lines
        line = readLine();
        std::stringstream sstr(line);
        sstr >> simBox.len[XX];
        sstr >> simBox.len[YY];
        sstr >> simBox.len[ZZ];

        // placeMolsIntoBox(conf, simBox);

        return true;
    }

};

#endif //GROFILEREADER_HPP