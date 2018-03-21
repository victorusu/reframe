/*
 * InputFile.cpp
 *
 *  Created on: Oct 1, 2010
 *      Author: victor
 */

#include "inputfile.hpp"

InputFile::InputFile() :
    FileHandler(std::fstream::in)
{
}

InputFile::~InputFile()
{

}

bool InputFile::isStringInFile(std::string str)
{
    std::string line;

    line = readLine();
    while(line.find(str) == std::string::npos && !eof())
        line = readLine();

    if(line.find(str) != std::string::npos)
        return true;
    else
        return false;
}
