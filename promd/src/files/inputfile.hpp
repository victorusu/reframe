/*
 * InputFile.h
 *
 *  Created on: Oct 1, 2010
 *      Author: victor
 */

#ifndef INPUTFILE_H_
#define INPUTFILE_H_

#include "filehandler.hpp"

class InputFile : public FileHandler
{

protected:

    std::fstream file;

public:

    InputFile();

    virtual ~InputFile();

    bool open() {
        file.open(fileName.c_str(), openMode);

        return isOpen();
    }

    bool open(std::string fileName) {

        this->fileName = fileName;
        return open();
    }

    void close() {

        try {
            file.close();
        } catch(std::exception e) {
            return;
        }
        // return !this->isOpen();
    }

    bool isOpen() const
    {
        return file.is_open();
    }

    std::string readLine()
    {
        std::string str;
        getline(file, str);

        return str;
    }

    std::string readLineSkippingCommentedLines(const char commentChar)
    {
        std::string line;
        bool found = false;

        while(!found && !eof()) {
            line = readLine();
            line.erase(0, line.find_first_not_of(' '));
            if(line.length() > 0) {
                if(line[0] != commentChar)
                    found = true;
            }
        }
        if(eof() || !found)
            return "";

        return line;
    }

    bool eof()
    {
        return file.eof();
    }

    bool isStringInFile(std::string str);

    bool rewind()
    {
        close();
        return open();
    }

    // template <class T>
    // T& operator>>(InputFile& input)
    // {
    //     // file << input;
    //     T result;
    //     input >> result;
    //     return result;
    // }

};

#endif /* INPUTFILE_H_ */
