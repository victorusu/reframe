/*
 * trcfilereader.cu
 *
 *  Created on: Dec 1, 2015
 *      Author: victor
 */

#include "trcfilewriter.hpp"

TRCFileWriter::TRCFileWriter()
{
    OutputFile();
}

TRCFileWriter::~TRCFileWriter()
{
    // FileHandler::close();
}
