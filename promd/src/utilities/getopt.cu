/*
 * getopt.cu
 *
 *  Created on: Dec 15, 2015
 *      Author: victor
 */
#include "getopt.hpp"

ArgumentOptions::ArgumentOptions(const std::string argumentName, const std::string description, const ArgumentOptions::Type type, const std::string defaultValue, const std::string help, const bool required, const bool holdsValue, const bool multipleValued)
    // : argumentName(argumentName), defaultValue(defaultValue), help(help), required(required), holdsValue(holdsValue), multipleValued(multipleValued)
{
    this->argumentName = argumentName;
    this->description = description;

    this->type = type;
    this->defaultValue = defaultValue;

    this->help = help;

    this->required = required;
    this->holdsValue = holdsValue;
    this->multipleValued = multipleValued;
}
