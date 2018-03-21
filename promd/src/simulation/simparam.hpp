#ifndef SIMPARAM_HPP
#define SIMPARAM_HPP

#include "../main.hpp"

struct SimParameters
{
    // time
    int nstlim;
    number dt;

    // degrees of freedom
    int com;
    int boxdof;

    // initial velocity
    int genvel;
    int ig;
    number tempi;

    // SHAKE
    int shake;
    number shaketol;

    // temperature
    int ntt;
    number temp0;
    number taut;

    // pressure
    int ntp;
    number compressibility;
    number taup;
    Tensor refPressureTensor;

    // nonbonded
    number rcutvdw;
    number rcutcoul;
    number rlist;
    int nstlist;

    // coulomb interaction
    number epscs;
    number epsrf;

    // printing
    int ntpr;
    int ntwx;
    int ntwf;

    // property
    int property;

    // physical constants
    number boltz;
    number vacuumpermittivity;
    number hbar;
    number speedoflight;


};



#endif //SIMPARAM_HPP