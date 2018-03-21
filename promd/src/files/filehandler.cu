/*
 * FileHandler.cpp
 *
 *  Created on: Oct 1, 2010
 *      Author: victor
 */

#include "filehandler.hpp"

FileHandler::FileHandler(std::fstream::openmode openMode)
{
    this->openMode = openMode;
}

FileHandler::~FileHandler()
{
}
