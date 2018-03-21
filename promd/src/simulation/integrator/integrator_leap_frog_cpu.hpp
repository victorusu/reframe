#ifndef INTEGRATOR_LEAP_FROG_CPU_HPP
#define INTEGRATOR_LEAP_FROG_CPU_HPP

#include "../../main.hpp"

#include "../../simulation/box.hpp"
#include "../../interaction/pairlist/pairlist.hpp"
#include "../../interaction/nonbonded/nonbonded.hpp"
#include "../../simulation/simparam.hpp"
#include "../../topology/topology.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif

#define mvsq2e 1;

struct Integrator {

    int step;

    // TODO remove the "molecular virial preparation" from this function
    VHR_INLINE
    void computeKineticEnergyAndTemperature(Configuration &conf, Topology &topology, const SimParameters &simParam, const number invBoltz)
    {
        startTimer(kineticEnergyTimer);

        int i, nAtoms = conf.nAtoms;
        number ekin = 0.0;
        number halfMass = 0.0;

        Tensor kineticTensor = Tensor(0.0);
        Tensor tempTensor = Tensor(0.0);



        // TODO
        // parallelise this loop
        // #pragma omp parallel for reduction(+:ekin)

        // This is a workaround to omp reduce
        number tensorXX, tensorXY, tensorXZ,
               tensorYX, tensorYY, tensorYZ,
               tensorZX, tensorZY, tensorZZ;

        // This is ugly!!!
        tensorXX=tensorXY=tensorXZ=tensorYX=tensorYY=tensorYZ=tensorZX=tensorZY=tensorZZ = 0.0;

        #pragma omp parallel for reduction(+: tensorXX, tensorXY, tensorXZ, \
               tensorYX, tensorYY, tensorYZ, \
               tensorZX, tensorZY, tensorZZ) schedule(guided, 1)
        for (i = 0; i < nAtoms; i++) {

            halfMass = 0.5 * topology.masses[i];

            tempTensor = (conf.current->v[i] * trans(conf.current->v[i])) * halfMass;

            tensorXX += tempTensor(0,0);
            tensorXY += tempTensor(0,1);
            tensorXZ += tempTensor(0,2);
            tensorYX += tempTensor(1,0);
            tensorYY += tempTensor(1,1);
            tensorYZ += tempTensor(1,2);
            tensorZX += tempTensor(2,0);
            tensorZY += tempTensor(2,1);
            tensorZZ += tempTensor(2,2);

            // kineticTensor += tempTensor;
            // ekin += trace(tempTensor);

            // inner product (v * v) time 1/2 mass
            // ekin += trans(conf.current->v[i]) * conf.current->v[i] * 0.5 * topology.masses[i];
            // ekin += (trans(conf.current->v[i]) * conf.current->v[i]) * halfMass;
            // ekin += trace(tempTensor);
        }

        kineticTensor(0,0) = tensorXX;
        kineticTensor(0,1) = tensorXY;
        kineticTensor(0,2) = tensorXZ;
        kineticTensor(1,0) = tensorYX;
        kineticTensor(1,1) = tensorYY;
        kineticTensor(1,2) = tensorYZ;
        kineticTensor(2,0) = tensorZX;
        kineticTensor(2,1) = tensorZY;
        kineticTensor(2,2) = tensorZZ;

        ekin = trace(kineticTensor);

        conf.current->temperature  = 2.0 * ekin / ((3.0*nAtoms)-simParam.boxdof) * invBoltz;
        conf.current->kineticEnergy = ekin;
        conf.current->kineticEnergyTensor = kineticTensor;

        // number kineticEnergy = 0.0;

        // conf.current->kineticEnergyTensor = 0.0;

        // // looping over molecules
        // for(mol = 0; mol < topology.nMolecules; mol++) {

        //     int first = topology.molecules.molsBegin[mol];
        //     int last = topology.molecules.molsEnd[mol];

        //     RVec comX = RVec(0.0);
        //     RVec comV = RVec(0.0);

        //     number totalKin = 0.0;
        //     number totalMass = 0.0;

        //     // looping over atoms i in molecule mol
        //     for(i = first; i < last; i++) {

        //         number mass = topology.masses[i];
        //         totalMass += mass;

        //         // finished gathering atom i in molecule mol
        //         comX += conf.current->x[i] * mass;
        //         comV += conf.current->v[i] * mass;

        //         // inner product
        //         totalKin += 0.5 * mass * trans(conf.current->v[i]) * conf.current->v[i];
        //     }
        //     comX /= totalMass;
        //     comV /= totalMass;

        //     const number halfTotalMass = 0.5 * totalMass;

        //     // saving the data

        //     // inner product
        //     const number transKin = halfTotalMass * trans(comV) * comV;

        //     topology.molecules.totalKin[mol] = totalKin;
        //     topology.molecules.transKin[mol] = transKin;
        //     topology.molecules.rotIntKin[mol] = totalKin - transKin;

        //     topology.molecules.comX[mol] = comX;
        //     topology.molecules.comV[mol] = comV;

        //     // outer product
        //     topology.molecules.kineticTensor[mol] = halfTotalMass * comV * trans(comV);

        //     conf.current->kineticEnergyTensor += topology.molecules.kineticTensor[mol];


        //     kineticEnergy += totalKin;
        // }
        // conf.current->temperature  = 2.0 * kineticEnergy / ((3.0*nAtoms)-simParam.boxdof) * invBoltz;
        // conf.current->kineticEnergy = kineticEnergy;

        addToTime(kineticEnergyTimer, kineticEnergyTime);
    }

    VHR_INLINE
    void propagateVelocities(Configuration &conf, const Topology &topology, const SimParameters &simParam)
    {
        startTimer(integrationTimer);

        int i, nAtoms = conf.nAtoms;
        // number ekin = 0.0;
        number dt = simParam.dt;
        // int boxdof = simParam.boxdof;
        // number invBoltz = 1.0 / simParam.boltz;

        // swap the positions and velocities
        conf.swap();
        // TODO
        // missing box swap ??? and molecule state swap???

        // conf.current->v = conf.old->v + (conf.f * dt * topology.invMasses);

        // #pragma omp parallel for reduction(+:ekin)
        #pragma omp parallel for
        for (i = 0; i < nAtoms; i++) {

            conf.current->v[i] = conf.old->v[i] + (conf.f[i] * (dt * topology.invMasses[i]));

            // inner product (v * v) time 1/2 mass
            // ekin += (trans(conf.current->v[i]) * conf.current->v[i]) * 0.5 * topology.masses[i];
        }
        // conf.current->temperature  = 2.0 * ekin / ((3.0 * nAtoms) - simParam.boxdof) * invBoltz;

        // conf.current->kineticEnergy = ekin;

        addToTime(integrationTimer, integrationTime);
    }

    // TODO
    // the maximum displacement should be done after shake, not after propagete coordinates
    VHR_INLINE
    void propagatePositions(Configuration & conf, PairList &pairlist, /*const*/ Box &simBox,
                            const number dt)
    {
        startTimer(integrationTimer);

        int i;
        const int nAtoms = conf.nAtoms;

        RVec box = simBox.len;
        RVec halfBox = 0.5 * box;
        RVec invBox;
        invBox = 1.0;
        invBox /= box;

        number displacement = 0.0;
        #pragma omp parallel for reduction(max: displacement)
        for (i = 0; i < nAtoms; i++) {

            conf.current->x[i] = conf.old->x[i] + (dt * conf.current->v[i]);

            // // put atom back into box (if needed)
            RVec shift = (conf.current->x[i] - halfBox) * invBox;
            conf.current->x[i] -= rvecNINT(shift) * box;

            if(pairlist.doMaxDisplacement) {
                RVec dx = simBox.pbc(pairlist.refPosition[i] - conf.current->x[i]);

                // inner product of dx * dx
                const number disp = std::sqrt(trans(dx) * dx);
                if(disp > displacement)
                    displacement = disp;
            }
        }
        pairlist.maxDisplacement = displacement;

        addToTime(integrationTimer, integrationTime);
    }


    // template <typename T>
    // struct summer
    // {
    //     T operator()(const T &lhs, const T &rhs) const
    //     {
    //         return lhs + rhs;
    //     }
    // };

    VHR_INLINE
    void removeCOMMotion(Configuration &conf, const Topology &topology, const Box &simBox)
    {
        number totalMass = 0.0;
        number comVx = 0.0, comVy = 0.0, comVz = 0.0;
        number comXx = 0.0, comXy = 0.0, comXz = 0.0;

        int i, j, nAtoms = conf.nAtoms;
        RVec comV, comX;
        comV = 0.0;
        comX = 0.0;

        // int first = topology.molecules.molsBegin[0];
        // int last = topology.molecules.molsEnd[0];

        // RVec box = simBox.len;
        // RVec halfBox = 0.5 * box;
        // RVec invBox;
        // invBox = 1.0;
        // invBox /= box;


        // // looping over molecules
        // for(i = first; i < last; i++) {

        //     number mass = topology.masses[i];
        //     totalMass += mass;

        //     comX += conf.current->x[i] * mass;

        // }
        // comX /= totalMass;
        // // bringing the com to the center of the box
        // comX -= simBox.len * 0.5;


        // totalMass = 0.0;;

        // // computing
        // #pragma omp declare reduction( + : RVec : std::transform(omp_in.begin( ),  omp_in.end( ), omp_out.begin( ), omp_out.begin( ), summer<RVec>() ) ) initializer (omp_priv=RVec())
        // #pragma omp declare reduction( + : RVec : std::transform(omp_in.begin( ),  omp_in.end( ), omp_out.begin( ), omp_out.begin( ), summer( )) ) initializer (omp_priv(omp_orig))


        #pragma omp parallel for reduction(+: totalMass, comVx, comVy, comVz, comXx, comXy, comXz)
        // #pragma omp parallel for reduction(+: totalMass, comVx, comVy, comVz)
        // #pragma omp parallel for reduction(+: totalMass) reduction(+: comV, comX)
        for (i = 0; i < nAtoms; i++) {
            number mass = topology.masses[i];

            totalMass += mass;

            // comV += conf.current->v[i] * mass;
            // comX += conf.current->x[i] * mass;

            comVx += conf.current->v[i][0] * mass;
            comVy += conf.current->v[i][1] * mass;
            comVz += conf.current->v[i][2] * mass;

            comXx += conf.current->x[i][0] * mass;
            comXy += conf.current->x[i][1] * mass;
            comXz += conf.current->x[i][2] * mass;

        }
        // // workaround for openmp reduction
        comX[0] = comXx;
        comX[1] = comXy;
        comX[2] = comXz;
        comX /= totalMass;

        // bringing the com to the center of the box
        comX -= simBox.len * 0.5;

        comV[0] = comVx;
        comV[1] = comVy;
        comV[2] = comVz;
        comV /= totalMass;

        // removing
        #pragma omp parallel for
        for (i = 0; i < nAtoms; i++) {
            conf.current->x[i] -= comX;
            conf.current->v[i] -= comV;

            // // // put atom back into box (if needed)
            // RVec shift = (conf.current->x[i] - halfBox) * invBox;
            // conf.current->x[i] -= rvecNINT(shift) * box;

        }

    }

};

#endif //INTEGRATOR_LEAP_FROG_CPU_HPP