#ifndef BLOCKSFILE_HPP
#define BLOCKSFILE_HPP

#include "../main.hpp"
#include "../simulation/simparam.hpp"
#include "../simulation/box.hpp"

#include "inputtraj.hpp"

#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <cstdlib>
#include <exception>


class BlocksFile : public InputTraj
{

public:

    BlocksFile();

    virtual ~BlocksFile();

protected:

    inline
    bool readLineSkippingCommentedLinesAndCheckEndOfBlock(std::string &line, const char commentChar, const std::string endBlockName)
    {
        bool found = false;

        while(!found && !eof()) {
            line = readLine();
            line.erase(0, line.find_first_not_of(' '));
            if(line.length() > 0) {
                if(line[0] != commentChar)
                    found = true;
            }
        }
        if(eof() || !found) {
            line = "";
            return false;
        }

        if(line == endBlockName)
            return false;

        return true;
    }

    inline
    bool findBlock(std::string blockName, const char commentChar = '#')
    {
        std::string line;
        rewind();

        line = readLineSkippingCommentedLines(commentChar);
        while(line.compare(blockName) != 0 && !eof()) {
            line = readLineSkippingCommentedLines(commentChar);
        }

        if(eof()) {
            fprintf(stderr, "Block %s not found in file %s\n", blockName.c_str(), getFileName().c_str());
            return false;
        }

        return true;
    }

};

#endif //BLOCKSFILE_HPP