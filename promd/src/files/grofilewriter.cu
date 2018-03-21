/*
 * grofilereader.cu
 *
 *  Created on: Dec 1, 2015
 *      Author: victor
 */

#include "grofilewriter.hpp"

GROFileWriter::GROFileWriter()
{
    OutputFile();
}

GROFileWriter::~GROFileWriter()
{
    // FileHandler::close();
}
