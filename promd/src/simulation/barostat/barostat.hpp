#ifndef BAROSTAT_HPP
#define BAROSTAT_HPP

#include "../../main.hpp"
#include "../../domain/domain.hpp"

#include "../box.hpp"
#include "../simparam.hpp"

#include "../../configuration/conf.hpp"
#include "../../topology/topology.hpp"

#include <iostream>

class Barostat
{

public:

    bool isotropic;
    bool computeAndApply;
    bool compute;

    // number refPressureTensor[9];
    // number pressureTensor[9];
    Tensor pressureTensor;
    // Tensor refPressureTensor;

    // number virialTensor[9];

    number compressibility;
    number invTauP;
    number timestepSize;

public:

    Barostat()
    {
        pressureTensor = 0.0;
    }

    number pressure()
    {
        if(isotropic)
            // return (pressureTensor[XXXX] + pressureTensor[YYYY] + pressureTensor[ZZZZ]) / 3.0;
            return trace(pressureTensor) / 3.0;
        else {
            cerr.printf("Anisotropic or full anisotropic or semi anisotropic pressure scaling not implemented yet\n");
            dd.abort();
        }
        return 0.0;
    }

    // void convertAtomicToMolecularVirial(Configuration &conf, Box &simBox, Topology &topology)
    // {
    //     int  i, mol;

    //     // #pragma omp parallel for
    //     for(mol = 0; mol < topology.nMolecules; mol++) {

    //         int first = topology.molecules.molsBegin[mol];
    //         int last = topology.molecules.molsEnd[mol];

    //         RVec comX = topology.molecules.comX[mol];
    //         number totalMass = 0.0;

    //         // looping over atoms i in molecule mol
    //         for(i = first; i < last; i++) {

    //             RVec dist = conf.current->x[i] - comX;

    //             number mass = topology.masses[i];

    //         }
    //         // finished gathering molecule
    //         comX /= totalMass;

    //         const number halfTotalMass = 0.5 * totalMass;

    //         // saving the data

    //         // inner product
    //         const number transKin = halfTotalMass * trans(comV) * comV;

    //         topology.molecules.totalKin[mol] = totalKin;
    //         topology.molecules.transKin[mol] = transKin;
    //         topology.molecules.rotIntKin[mol] = totalKin - transKin;

    //         topology.molecules.comX[mol] = comX;

    //         // outer product
    //         topology.molecules.kineticTensor[mol] = halfTotalMass * comV * trans(comV);

    //         kineticEnergy += totalKin;
    //     }
    // }

    virtual void pressureCalculation(Configuration &conf, Box &simBox, Topology &topology) {}

    virtual void scaleBoxAndCoordinates(Configuration &conf, Box &simBox, const SimParameters &simParam, PairList &pairlist){};

};


class BerendsenBarostat : public Barostat
{
    void pressureCalculation(Configuration &conf, Box &simBox, Topology &topology)
    {
        // int i, mol, nmols = topology.nMolecules;
        if(isotropic) {

            // this is a atomic virial
            pressureTensor = ((2.0 * conf.current->kineticEnergyTensor) + conf.current->virialTensor) / simBox.volume();

            // pressureTensor = 0.0;

            // // looping over molecules
            // for(mol = 0; mol < nmols; mol++)
            //     pressureTensor += ((2.0 * topology.molecules.kineticTensor[mol]) + topology.molecules.virialTensor[mol]) / simBox.volume();

        } else {
            cerr.printf("Anisotropic or full anisotropic or semi anisotropic pressure scaling not implemented yet\n");
            dd.abort();
        }
    }

    void scaleBoxAndCoordinates(Configuration &conf, Box &simBox, const SimParameters &simParam, PairList &pairlist)
    {
        int i, nAtoms = conf.nAtoms;
        if(isotropic) {
            // TODO
            // if add positional restraints this should be modified

            // computing the total pressure
            number totalPressure = pressure();

            if(computeAndApply) {

                // computing the scale factor
                number mu = std::pow(1.0 - (this->compressibility * this->timestepSize * this->invTauP * (simParam.refPressureTensor(XX, XX) - totalPressure)), 1.0/3.0);

                // scaling the box
                simBox.scaleBy(mu);

                // scaling the pairlist also
                if(pairlist.doMaxDisplacement)
                    pairlist.maxDisplacement *= mu;

                // scaling the positions
                #pragma omp parallel for
                for (i = 0; i < nAtoms; i++) {
                    conf.current->x[i] *= mu;

                    // scaling the pairlist displacement also
                    if(pairlist.doMaxDisplacement)
                        pairlist.refPosition[i] *= mu;
                }
            }

        } else {
            cerr.printf("Anisotropic or full anisotropic or semi anisotropic pressure scaling not implemented yet\n");
            dd.abort();
        }
    }
};

VHR_INLINE
bool allocateBarostat(Barostat **barostat, const SimParameters &simParam)
{
    int i;
    switch(simParam.ntp) {

        case 2:
            *barostat = new BerendsenBarostat();

            // for(i = 0; i < 9; i++)
            //  (*barostat)[0].refPressureTensor[i] = simParam.refPressureTensor[i];
            // (*barostat)[0].refPressureTensor = simParam.refPressureTensor;

            (*barostat)[0].compressibility = simParam.compressibility;
            (*barostat)[0].invTauP = 1.0 / simParam.taup;
            (*barostat)[0].timestepSize = simParam.dt;

            (*barostat)[0].isotropic = true;
            (*barostat)[0].computeAndApply = true;
            (*barostat)[0].compute = true;
        break;

        case 0:
            *barostat = new Barostat();
            (*barostat)[0].isotropic = true;
            (*barostat)[0].computeAndApply = false;
            (*barostat)[0].compute = false;
        break;

        case 1:
            *barostat = new BerendsenBarostat();

            // for(i = 0; i < 9; i++)
            //     (*barostat)[0].refPressureTensor[i] = simParam.refPressureTensor[i];
            // (*barostat)[0].refPressureTensor = simParam.refPressureTensor;

            (*barostat)[0].compressibility = simParam.compressibility;
            (*barostat)[0].invTauP = 1.0 / simParam.taup;
            (*barostat)[0].timestepSize = simParam.dt;

            (*barostat)[0].isotropic = true;
            (*barostat)[0].computeAndApply = false;
            (*barostat)[0].compute = true;

        break;

        default:
            return false;
        break;
    }

    return true;
}

#endif //BAROSTAT_HPP
