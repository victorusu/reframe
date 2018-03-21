#ifndef IMDFILEREADER_HPP
#define IMDFILEREADER_HPP

#include "../main.hpp"
#include "../domain/domain.hpp"
#include "../simulation/simparam.hpp"
#include "../simulation/box.hpp"

#include "inputfile.hpp"

#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <cstdlib>
#include <exception>


class IMDFileReader : public InputFile
{
#define NFIELDS 25

public:

    IMDFileReader();

    virtual ~IMDFileReader();

    // TODO re-implement properly
    bool readfile(SimParameters &simParam)
    {
        std::string line;
        std::vector<number> data;
        data.reserve(NFIELDS);

        while(!eof()) {
            line = readLineSkippingCommentedLines('#');
            if(line.length() > 0) {
                std::stringstream sstr(line);
                while(!sstr.eof()) {
                    number val;
                    sstr >> val;
                    data.push_back(val);
                }
            }
        }

        if(data.size() != NFIELDS) {
            fprintf(stderr, "Number of fields in the imd file (%zu) does not meet the expected number (%d)\n", data.size(), NFIELDS);
            return false;
        }

        simParam.nstlim   = data[0];
        simParam.dt       = data[1];

        simParam.com      = data[2];

        // We only support rectangular boxes!!!
        // if(simParam.com != 0)
        simParam.boxdof = 3;
        // else
        //     simParam.boxdof = 0;

        simParam.genvel   = data[3];
        simParam.ig       = data[4];
        simParam.tempi    = data[5];

        simParam.shake    = data[6];
        simParam.shaketol = data[7];

        simParam.ntt      = data[8];
        simParam.temp0    = data[9];
        simParam.taut     = data[10];

        simParam.ntp      = data[11];
        simParam.compressibility = data[12];
        simParam.taup     = data[13];

        simParam.refPressureTensor(XX, XX) = data[14];
        simParam.refPressureTensor(XX, YY) = 0.0;
        simParam.refPressureTensor(XX, ZZ) = 0.0;

        simParam.refPressureTensor(YY, XX) = 0.0;
        simParam.refPressureTensor(YY, YY) = data[14];
        simParam.refPressureTensor(YY, ZZ) = 0.0;

        simParam.refPressureTensor(ZZ, XX) = 0.0;
        simParam.refPressureTensor(ZZ, YY) = 0.0;
        simParam.refPressureTensor(ZZ, ZZ) = data[14];

        simParam.rcutvdw  = data[15];
        simParam.rcutcoul = data[16];
        simParam.rlist    = data[17];
        simParam.nstlist  = data[18];

        simParam.epscs    = data[19];
        simParam.epsrf    = data[20];

        simParam.ntpr     = data[21];
        simParam.ntwx     = data[22];
        simParam.ntwf     = data[23];

        simParam.property = data[24];

        return true;
    }

};

inline bool checkSimParam(SimParameters &simParam, Box &simBox)
{
    if(simBox.len[XX] == 0 || simBox.len[YY] == 0 || simBox.len[ZZ] == 0) {
        cerr.printf("BOX[0..2] == 0\n");
        return false;
    }

    if(simParam.ntt > 0 && simParam.taut <= 0) {
        cerr.printf("NTT > 0 and TAUT <= 0\n");
        return false;
    }

    // if(simParam.boltz <= 0) {
    //     std::cerr << "BOLTZ == 0\n" ;
    //     return false;
    // }

    if(2*simParam.rcutvdw > simBox.len[XX] || 2*simParam.rcutvdw > simBox.len[YY] || 2*simParam.rcutvdw > simBox.len[ZZ]) {
        cerr.printf("2*RCUT-VDW > BOX[0..2]\n");
        return false;
    }

    if(2*simParam.rcutcoul > simBox.len[XX] || 2*simParam.rcutcoul > simBox.len[YY] || 2*simParam.rcutcoul > simBox.len[ZZ]) {
        cerr.printf("2*RCUT-COUL > BOX[0..2]\n");
        return false;
    }

    if(simParam.rcutvdw > simParam.rlist) {
        cerr.printf("RCUT-VDW > RLIST. Twin cut-off not supported yet\n");
        return false;
    }
    else if(simParam.rcutcoul > simParam.rlist) {
        cerr.printf("RCUT-COUL > RLIST. Twin cut-off not supported yet\n");
        return false;
    }
    else if(simParam.rcutvdw == simParam.rlist && simParam.nstlist < -1) {
        cerr.printf("RCUT-VDW == RLIST AND NSTLIST < 0\n This will update the pairlist everystep. I am not sure if it is what you want.\n");
        return false;
    }
    else if(simParam.rcutcoul == simParam.rlist && simParam.nstlist < -1) {
        cerr.printf("RCUT-VDW == RLIST AND NSTLIST < 0\n This will update the pairlist everystep. I am not sure if it is what you want.\n");
        return false;
    }


    // if(simParam.ntpr == 0 || simParam.ntwx == 0) {
    //     std::cerr << "NPTR == 0 or NTWX == 0\n" ;
    //     return false;
    // }

    return true;
}


#endif //IMDFILEREADER_HPP