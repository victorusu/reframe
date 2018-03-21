/*
 * imdfilereader.cu
 *
 *  Created on: Dec 1, 2015
 *      Author: victor
 */

#include "imdfilereader.hpp"

IMDFileReader::IMDFileReader()
{
    InputFile();
}

IMDFileReader::~IMDFileReader()
{
    // FileHandler::close();
}
