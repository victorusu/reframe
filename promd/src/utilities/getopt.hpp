#ifndef GETOPT_HPP
#define GETOPT_HPP

#include "../main.hpp"
#include "stdout.hpp"

// #include "../domain/domain.hpp"
#include "../files/filehandler.hpp"

#include <string>
#include <sstream>
#include <typeinfo>
#include <algorithm>

#include <getopt.h>

class ArgumentOptions
{

public:

    enum Type {NUMBER, INT, BOOL, STRING};

    static char *typenames;

    std::string typeName;

    Type type;
    std::string defaultValue;
    std::vector<std::string> values;

    // std::string longName;
    std::string argumentName;

    std::string help;
    std::string description;

    bool required;
    bool holdsValue;
    bool multipleValued;

    ArgumentOptions(const std::string argumentName, const std::string description, const ArgumentOptions::Type type, std::string defaultValue,
                    const std::string help, const bool required, const bool holdsValue, const bool multipleValued);

    bool valid()
    {
        if(!defaultValue.empty() && !argumentName.empty() && !description.empty())
            return true;
        return false;
    }

    void print()
    {
        std::string str;

        switch(type) {
            case NUMBER:
                str = "number";
            break;
            case INT:
                str = "integer";
            break;
            case BOOL:
                str = "bool";
            break;
            case STRING:
                str = "string";
            break;
            default:
                str = "unknown";
            break;
        }

        cerr.printf("%15s %15s %-10s  %-30s\n", argumentName.c_str(), defaultValue.c_str(), str.c_str(), description.c_str());
    }

};

struct GetOpt {

    int argc;
    std::vector<std::string> argv;
    std::vector<ArgumentOptions> options;
    std::string help;

    // DomainDecomposition *dd;

    GetOpt(const int argc, char**argv)
    {
        int i;
        for(i = 0; i < argc; i++)
            this->argv.push_back(argv[i]);

        this->argc = argc;

        options.reserve(argc);
        help = "";

        *this << ArgumentOptions("-h", "this help", ArgumentOptions::BOOL, "false", "help", false, false, false);

        // this->dd = &dd;
    }

    GetOpt & operator << (ArgumentOptions opt)
    {
      options.push_back(opt);
      return *this;
    }

    void printUserArguments()
    {
        int i;
        for(i = 0; i < argc; i++)
            cerr.printf("%s ", argv[i].c_str());
        cerr.printf("\n");
        cerr.flush();
    }

    void printOptions()
    {
        size_t i;
        cerr.printf("%15s %15s %-10s  %-30s\n", "Argument Name", "Default Value", "Type", "Description");
        for(i = 0; i < options.size(); i++) {
            options[i].print();
        }
    }

    void printHelp()
    {
        cerr.printf("%s\n", help.c_str());
        printOptions();
    }

    std::string parse(std::string optionName)
    {
        std::string result = "";
        if(!userTypedOption(optionName)) {
            return findDefaultValueOfOption(optionName);
        }
        else
        {
            result = getUserArgument(optionName);
        }
        return result;
    }

private:

    ArgumentOptions findOption(const std::string& optionName)
    {
        size_t i;
        for(i = 0; i < options.size(); i++) {
            if(options[i].argumentName == optionName)
                return options[i];
        }
        return ArgumentOptions("", "", ArgumentOptions::STRING, "", "", false, false, false);
    }

    bool userTypedOption(const std::string& optionName)
    {
        return std::find(argv.begin(), argv.end(), optionName) != argv.end();
    }

    std::string getUserArgument(const std::string & optionName)
    {
        std::vector<std::string>::iterator itr = std::find(argv.begin(), argv.end(), optionName);

        if (itr != argv.end())
        {
            ArgumentOptions opt = findOption(optionName);
            if(!opt.holdsValue) {
                return *itr;
            }
            else if(opt.multipleValued) {
                // TODO
            }
            else if(++itr != argv.end()) {
                return *itr;
            }
        }
        return std::string();
    }

    std::string findDefaultValueOfOption(const std::string& optionName)
    {
        size_t i;
        for(i = 0; i < options.size(); i++) {
            if(options[i].argumentName == optionName)
                return options[i].defaultValue;
        }
        return "";
    }



    // template<typename T>
    // std::vector<T> parse(std::string optionName)
    // {
    //     for(i = 0; i < argc; i++) {
    //         if(!longArgumentOptions[i].holdsValue) {
    //             // T result = longArgumentOptions[i].values;
    //             T result =
    //         }
    //     }

    //     // T val = option
    //     int i;
    //     // int c = getopt_long(argc, argv, "abc:d:f:", longArgumentOptions, &optionsIndex);

    //     return std::vector<std::string>();
    // }


};


#endif //GETOPT_HPP
