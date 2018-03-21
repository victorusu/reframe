/*
 * FileHandler.h
 *
 *  Created on: Oct 1, 2010
 *      Author: victor
 */

#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include "../main.hpp"
#include "../domain/domain.hpp"

#include <string>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>

#include <sys/stat.h>

#include <stdarg.h>

class FileHandler
{

protected:

    std::string fileName;

    std::fstream::openmode openMode;

// protected:

    // std::fstream file;
    // FILE *file;

public:

    FileHandler(std::fstream::openmode openMode);

    virtual ~FileHandler();

    std::string getFileName() const
    {
        return fileName;
    }

    std::string getFileNameWithoutExt()
    {
        size_t pos = fileName.rfind('.');

        return fileName.substr(0,pos);
    }

    static
    std::string getFileNameWithoutExt(std::string fileName)
    {
        size_t pos = fileName.rfind('.');

        return fileName.substr(0,pos);
    }

    std::string getFileNameExt()
    {
        size_t pos = fileName.rfind('.');

        return fileName.substr(pos+1,fileName.length());
    }

    static
    std::string getFileNameExt(std::string fileName)
    {
        size_t pos = fileName.rfind('.');

        return fileName.substr(pos+1,fileName.length());
    }

    void setFileName(std::string fileName)
    {
        this->fileName = fileName;
    }

};


#endif /* FILEHANDLER_H_ */
