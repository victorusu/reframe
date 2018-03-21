#ifndef STDOUT_HPP
#define STDOUT_HPP

// #include "../domain/domain.hpp"

#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <string.h>

#include <stdarg.h>

struct Stderr {

    // DomainDecomposition *dd;

    // inline
    // void setDomainDecomposition(DomainDecomposition &dd)
    // {
    //     this->dd = &dd;
    // }

    inline
    int printf (const char *format, ...)
    {
        // if(this->dd->master()) {
            va_list arg;
            int done;

            va_start (arg, format);
            done = vfprintf(stderr, format, arg);
            va_end (arg);

            return done;
        // }
        // return 0;
    }

    inline
    void flush()
    {
        fflush(stderr);
    }

};

struct Stdout {

    // DomainDecomposition *dd;

    // inline
    // void setDomainDecomposition(DomainDecomposition &dd)
    // {
    //     this->dd = &dd;
    // }

    inline
    int printf (const char *format, ...)
    {
        // if(this->dd->master()) {
            va_list arg;
            int done;

            va_start (arg, format);
            done = vfprintf(stdout, format, arg);
            va_end (arg);

            return done;
        // }
        // return 0;
    }

    inline
    void flush()
    {
        fflush(stdout);
    }

};

extern Stderr cerr;
extern Stdout cout;

#endif //STDOUT_HPP