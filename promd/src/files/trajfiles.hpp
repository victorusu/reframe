#ifndef TRAJFILES_HPP
#define TRAJFILES_HPP

#include "../main.hpp"

#include "cnffilereader.hpp"
#include "grofilereader.hpp"

#include "grofilewriter.hpp"
#include "cnffilewriter.hpp"
#include "g96filewriter.hpp"
#include "trcfilewriter.hpp"

#include "trffilewriter.hpp"

#include <string>
#include <algorithm>


inline bool allocateInputTRAJFile(InputTraj **inputTRAJFile, std::string fileName)
{
    std::string ext = FileHandler::getFileNameExt(fileName);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if(ext == "cnf") {
        *inputTRAJFile = new CNFFileReader();
    }
    else if(ext == "gro") {
        *inputTRAJFile = new GROFileReader();
    }
    else {
        return false;
    }

    (*inputTRAJFile)[0].setFileName(fileName);
    return true;
}


inline bool allocateOutputTRAJFile(OutputTraj **outputTRAJFile, std::string fileName)
{
    std::string ext = FileHandler::getFileNameExt(fileName);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if(ext == "gro") {
        *outputTRAJFile = new GROFileWriter();
    }
    else if(ext == "g96") {
        *outputTRAJFile = new G96FileWriter();
    }
    else if(ext == "cnf") {
        *outputTRAJFile = new CNFFileWriter();
    }
    else if(ext == "trc") {
        *outputTRAJFile = new TRCFileWriter();
    }
    else {
        return false;
    }

    (*outputTRAJFile)[0].setFileName(fileName);
    return true;
}

#endif //TRAJFILES_HPP