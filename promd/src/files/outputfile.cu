/*
 * OutputFile.cpp
 *
 *  Created on: Oct 6, 2010
 *      Author: victor
 */

#include "outputfile.hpp"

OutputFile::OutputFile() :
    FileHandler(std::fstream::out)
{
	this->file = NULL;
}

OutputFile::OutputFile(std::fstream::openmode openMode) :
    FileHandler(openMode)
{
	this->file = NULL;
}

OutputFile::~OutputFile()
{
	// close();
    // FileHandler::close();
}
