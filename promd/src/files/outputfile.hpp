/*
 * OutputFile.h
 *
 *  Created on: Oct 6, 2010
 *      Author: victor
 */

#ifndef OUTPUTFILE_H_
#define OUTPUTFILE_H_

#include "filehandler.hpp"

#include <stdarg.h>

class OutputFile : public FileHandler
{

protected:

    FILE *file;

public:

    OutputFile();

    OutputFile(std::fstream::openmode openMode);

    virtual ~OutputFile();

    bool open() {

        char op = 'w';
        switch(openMode) {

            case std::fstream::app:
                op = 'a';
            break;
            case std::fstream::out:
            default:
                op = 'w';
            break;
        }
        this->file = fopen(fileName.c_str(), &op);

        return isOpen();
    }

    bool open(std::string fileName) {

        this->fileName = fileName;

        return open();
    }

    void close() {

        try {
            fclose(this->file);
        } catch(std::exception e) {
            return;
        }
        // return !this->isOpen();
    }

    bool isOpen() const
    {
        return this->file != NULL;
    }

    void flush()
    {
        fflush(file);
    }

    // FILE *buffer()
    // {
    // 	return this->file;
    // }

    // template <class T>
    // std::fstream& operator<<(T input)
    // {
    //     file << input;
    //     return file;
    // }



};

#endif /* OUTPUTFILE_H_ */
