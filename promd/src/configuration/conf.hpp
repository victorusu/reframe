#ifndef CONF_HPP
#define CONF_HPP

#include "../main.hpp"  // For function prototypes and run-time constants.
#include "../simulation/box.hpp"

struct StateVariables
{
    RVector x;
    RVector v;
    // RVector vTG;
    // RVector comVTG;

    number kineticEnergy;
    number temperature;

    Tensor virialTensor;
    Tensor kineticEnergyTensor;


};

struct Configuration {

private:
    StateVariables stateA, stateB;

public:

    std::string title;
    int nAtoms;

    StateVariables * old;
    StateVariables * current;

    RVector f;

    void swap()
    {
        StateVariables * tmp = current;
        current = old;
        old = tmp;
    }

    Configuration()
    {
        old = &stateA;
        current = &stateB;
    }

    void reserve(int n, int nthreads)
    {
        old->x.reserve(n);

        old->v.reserve(n);

        current->x.reserve(n);

        current->v.reserve(n);

        f.reserve(nthreads*n);
    }

    void correctForceAllocation(int nthreads)
    {
#if defined(Thrust)
        return;
#else
        if(nthreads > 1) {
            int size = f.size();

            f.resize(nthreads*size, 0.0);
            // this->fy.resize(nthreads*size, 0.0);
            // this->fz.resize(nthreads*size, 0.0);
        }
#endif
    }

    void correctAllocation(int nthreads)
    {
        // correcting the force vector allocation
        // if(nthreads > 1) {
        //     int size = f.size();
        //     f.resize(nthreads*size);
        // }
        if(f.size() != nAtoms)
            f.resize(nthreads * nAtoms);

        if(current->x.size() != nAtoms)
            current->x.resize(nAtoms);

        if(old->x.size() != nAtoms)
            old->x.resize(nAtoms);

        if(current->v.size() != nAtoms)
            current->v.resize(nAtoms);

        if(old->v.size() != nAtoms)
            old->v.resize(nAtoms);

    }
};


#endif //CONF_HPP