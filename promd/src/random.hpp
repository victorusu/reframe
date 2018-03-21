#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "main.hpp"

const int M      =  100000000;
const int M1     =  10000;
const int MULT   =  31415821;

#define abs(x)  (((x) > 0) ? (x) : -(x))

number a_random(int *seed)
{
    /*
    static int      irand, new = 1;
    */
    int irand;
    int irandh, irandl, multh, multl;
    number r;

    irand  = abs(*seed) % M;
    irandh = (int)(irand/M1);
    irandl = irand%M1;
    multh  = (int)(MULT/M1);
    multl  = MULT%M1;
    irand  = ((irandh*multl + irandl*multh)%M1)*M1 + irandl*multl;
    irand  = (irand+1)%M;
    r = (number)(irand/10.)*10./(number)M;
    if(r <= 0. || r > 1.)
        r = 0.;

    *seed = irand;
    return r;
}

#endif //RANDOM_HPP
