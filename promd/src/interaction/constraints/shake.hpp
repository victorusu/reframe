#ifndef SHAKE_HPP
#define SHAKE_HPP

#include "../../main.hpp"
#include "../../domain/domain.hpp"

#include "../../configuration/conf.hpp"
#include "../../simulation/box.hpp"
#include "../../topology/topology.hpp"

struct SHAKE
{
    std::vector<std::vector<int> > bondsPerMolecule;

    void prepare(Topology &topology)
    {
        int i;

        const int nBonds = topology.bonds.nBonds;
        const int nMols = topology.nMolecules;

        this->bondsPerMolecule.resize(nMols);

        // #pragma omp parallel for
        for(i = 0; i < nBonds; i++) {

            const int ii = topology.bonds.atoms[2*i];
            const int jj = topology.bonds.atoms[(2*i) + 1];

            const int mol = topology.molecules.atomsToMolecules[ii];
            if(mol != topology.molecules.atomsToMolecules[jj]) {
                cerr.printf("Error. We do not support SHAKE between different molecules\n");
                dd.abort();
            }

            this->bondsPerMolecule[mol].push_back(i);
        }
    }

    VHR_INLINE
    bool apply(Configuration &conf, const Topology &topology, const SimParameters &simParam, Box &simBox)
    {
        const int NSHAKEITERATIONS = 50000;

        // When we apply SHAKE we apply for all bonds!
        int i; //, nbonds = topology.bonds.nBonds;
        int mol, nMols = topology.nMolecules;


        // number tolerance = simParam.tolerance;
        const number tolerance = 2.0 * simParam.shaketol;

        const number dt2 = simParam.dt * simParam.dt;

        bool notError = true;


        // This is a workaround to omp reduce
        number tensorXX, tensorXY, tensorXZ,
               tensorYX, tensorYY, tensorYZ,
               tensorZX, tensorZY, tensorZZ;

        // This is ugly!!!
        tensorXX=tensorXY=tensorXZ=tensorYX=tensorYY=tensorYZ=tensorZX=tensorZY=tensorZZ = 0.0;

        #pragma omp parallel for reduction(+: tensorXX, tensorXY, tensorXZ, \
               tensorYX, tensorYY, tensorYZ, \
               tensorZX, tensorZY, tensorZZ) schedule(guided, 1)
        for(mol = 0; mol < nMols; mol++) {

            Tensor virialTensor;
            virialTensor = 0.0;

            bool iterateSHAKE = true;
            int iterations = 0;

            while(iterateSHAKE && notError) {

                iterateSHAKE = false;

                const int nBonds = this->bondsPerMolecule[mol].size();

                // When we apply SHAKE we apply for all bonds
                for(i = 0; i < nBonds; i++) {

                    const int bond = this->bondsPerMolecule[mol][i];

                    const int ii = topology.bonds.atoms[2*bond];
                    const int jj = topology.bonds.atoms[(2*bond) + 1];
                    const int type = topology.bonds.type[bond];

                    const number b0 = topology.bonds.b0[bond];

                    // TODO
                    // correct this! This is a huge cache miss!!!
                    // const number b0 = topology.bondStretchTypes.b0[type];
                    const number b0Sqr = b0 * b0;

                    const RVec rij = simBox.pbc(conf.current->x[ii] - conf.current->x[jj]);

                    const number rijSqr = (trans(rij) * rij);

                    const number diffSqr = b0Sqr - rijSqr;

                    if(fabs(diffSqr) >=  tolerance * b0Sqr) {

                        iterations++;

                        iterateSHAKE = true;
                        // cerr.printf("where are shaking!!!\n");

                        // we are constraining!!!!
                        const RVec rijOld = simBox.pbc(conf.old->x[ii] - conf.old->x[jj]);

                        // inner product
                        const number dotProduct = (trans(rijOld) * rij);

                        // checking shake error
                        if(dotProduct < b0Sqr * EPSILONLIMIT) {
                            cerr.printf("\n\nSHAKE ERROR between atoms %d and %d\n", ii+1, jj+1);
                            #pragma omp atomic write
                            notError = false;
                        }


                        if(iterations > NSHAKEITERATIONS) {
                            cerr.printf("\n\nSHAKE ERROR on molecule %d\n", mol + 1);
                            cerr.printf("Number of SHAKE iterations exceeded %d\n", NSHAKEITERATIONS);
                            #pragma omp atomic write
                            notError = false;
                        }

                        const number invMassii = topology.invMasses[ii];
                        const number invMassjj = topology.invMasses[jj];

                        const number lambda = diffSqr / (2.0 * dotProduct * (invMassii + invMassjj));

                        const RVec constraintForce = rijOld * lambda;

                        conf.current->x[ii] += constraintForce * invMassii;
                        conf.current->x[jj] -= constraintForce * invMassjj;

                        // conf.old->virialTensor += (rijOld * trans(rijOld)) * (lambda / dt2);
                        virialTensor = (rijOld * trans(rijOld)) * (lambda / dt2);
                        tensorXX += virialTensor(0,0);
                        tensorXY += virialTensor(0,1);
                        tensorXZ += virialTensor(0,2);
                        tensorYX += virialTensor(1,0);
                        tensorYY += virialTensor(1,1);
                        tensorYZ += virialTensor(1,2);
                        tensorZX += virialTensor(2,0);
                        tensorZY += virialTensor(2,1);
                        tensorZZ += virialTensor(2,2);
                    }
                }
            }
        }


        if(!notError)
            return false;

        Tensor virialTensor;
        virialTensor(0,0) = tensorXX;
        virialTensor(0,1) = tensorXY;
        virialTensor(0,2) = tensorXZ;
        virialTensor(1,0) = tensorYX;
        virialTensor(1,1) = tensorYY;
        virialTensor(1,2) = tensorYZ;
        virialTensor(2,0) = tensorZX;
        virialTensor(2,1) = tensorZY;
        virialTensor(2,2) = tensorZZ;

        conf.old->virialTensor += virialTensor;


        // cerr.printf("number of SHAKE iterations: %d\n", iterations);

        // TODO
        // correct the constraint velocities!!!
        const number invTimeStep = 1.0 / simParam.dt;
        #pragma omp parallel for
        for(i = 0; i < topology.constrainedAtoms.size(); i++) {

            int ii = topology.constrainedAtoms[i];

            conf.current->v[ii] = simBox.pbc(conf.current->x[ii] - conf.old->x[ii]) * invTimeStep;
        }

        return true;
    }

};






#endif //SHAKE_HPP
