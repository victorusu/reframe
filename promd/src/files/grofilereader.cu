/*
 * grofilereader.cu
 *
 *  Created on: Dec 1, 2015
 *      Author: victor
 */

#include "grofilereader.hpp"

GROFileReader::GROFileReader()
{
    InputFile();
}

GROFileReader::~GROFileReader()
{
    // FileHandler::close();
}
