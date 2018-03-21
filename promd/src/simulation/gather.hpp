#ifndef GATHER_HPP
#define GATHER_HPP

#include "../main.hpp"  // For function prototypes and run-time constants.
#include "../simulation/box.hpp"
#include "../configuration/conf.hpp"
#include "../topology/topology.hpp"

void gather(Configuration &conf, Box &simBox, Topology &topology)
{
    constexpr const number PARAMETIC_DISTANCE = 0.8;
    int  i, mol;

    // gathering
    for(mol = 0; mol < topology.nMolecules; mol++) {

        int first = topology.molecules.molsBegin[mol] + 1;
        int last = topology.molecules.molsEnd[mol];

        for(i = first; i < last; i++) {

            int j;
            number dist;

            for(j = 0; j < DIM; j++) {

                dist = conf.current->x[i][j] - conf.current->x[i-1][j];

                if(fabs(dist) > PARAMETIC_DISTANCE * simBox.len[j]) {
                    if(dist > 0) {
                        conf.current->x[i][j] -= simBox.len[j];
                    }
                    else {
                        conf.current->x[i][j] += simBox.len[j];
                    }
                }
            }
        }
    }
}

void computePerMoleculeCOMAndCOMV(Configuration &conf, Box &simBox, Topology &topology, const SimParameters &simParam, const number invBoltz)
{
    const number PARAMETIC_DISTANCE = 0.8;
    int  i, mol, nMols = topology.nMolecules;

    number kineticEnergy = 0.0;

    // TODO
    // paralelize this loop

    // looping over molecules
    for(mol = 0; mol < nMols; mol++) {

        int first = topology.molecules.molsBegin[mol] + 1;
        int last = topology.molecules.molsEnd[mol];

        RVec comX, comV;
        comX = 0.0;
        comV = 0.0;

        number totalKin = 0.0;
        number totalMass = 0.0;

        // looping over atoms i in molecule mol
        for(i = first; i < last; i++) {

            int j;
            number dist;

            number mass = topology.masses[i];
            totalMass += mass;

            // looping over dimension in order to gather atom i
            for(j = 0; j < DIM; j++) {

                dist = conf.current->x[i][j] - conf.current->x[i-1][j];

                if(fabs(dist) > PARAMETIC_DISTANCE * simBox.len[j]) {
                    if(dist > 0) {
                        conf.current->x[i][j] -= simBox.len[j];
                    }
                    else {
                        conf.current->x[i][j] += simBox.len[j];
                    }
                }
            }

            // finished gathering atom i in molecule mol
            comX += conf.current->x[i] * mass;
            comV += conf.current->v[i] * mass;
            totalKin += 0.5 * mass * trans(conf.current->v[i]) * conf.current->v[i];


        }
        // finished gathering molecule
        comX /= totalMass;
        comV /= totalMass;

        const number halfTotalMass = 0.5 * totalMass;

        // saving the data

        // inner product
        const number transKin = halfTotalMass * trans(comV) * comV;

        topology.molecules.totalKin[mol] = totalKin;
        topology.molecules.transKin[mol] = transKin;
        topology.molecules.rotIntKin[mol] = totalKin - transKin;

        topology.molecules.comX[mol] = comX;
        topology.molecules.comV[mol] = comV;

        // outer product
        topology.molecules.kineticTensor[mol] = halfTotalMass * comV * trans(comV);

        kineticEnergy += totalKin;
    }

    conf.current->kineticEnergy += kineticEnergy;
    conf.current->temperature  = 2.0 * kineticEnergy / ((3.0 * conf.nAtoms)-simParam.boxdof) * invBoltz;

}


#endif //GATHER_HPP