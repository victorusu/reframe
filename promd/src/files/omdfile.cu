/*
 * omdfile.cu
 *
 *  Created on: Oct 6, 2010
 *      Author: victor
 */

#include "omdfile.hpp"

OMDFile::OMDFile() :
	OutputFile()
{

}

OMDFile::OMDFile(std::fstream::openmode openMode) :
    OutputFile(openMode)
{
}

OMDFile::~OMDFile()
{
    // OutputFile::close();
}
