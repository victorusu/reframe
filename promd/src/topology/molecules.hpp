#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include "../domain/domain.hpp"            // For function prototypes and run-time constants.
#include "topologyfields.hpp"

// struct MolecularInfo
// {
//     // center of mass velocities
//     RVector comX;
//     RVector comKinEnergy;

//     std::vector<number> totalKin, transKin, rotIntKin;
// };

struct Molecules {

    // MolecularInfo *molInfoAtTMinusDtoverTwo;
    // MolecularInfo *molInfoAtTPlusDtoverTwo;

    std::vector<int> molsBegin;
    std::vector<int> molsEnd;

    RVector comX;
    RVector comV;

    RVector totalKin, transKin, rotIntKin;

    TVector kineticTensor;
    TVector virialTensor;

    std::vector<int> atomsToMolecules;


    // number totalKineticEnergyOfBath;

    // Molecules()
    // {
    //     molInfoAtTMinusDtoverTwo = &infoA;
    //     molInfoAtTPlusDtoverTwo = &infoB;
    // }

    // VHR_INLINE
    // void gatherConfIntoBoxAndRemoveCOMMotion(Configuration &conf, const Topology &topology, const Box &simBox)
    // {
    //     int i;

    //     // loop over molecules
    //     #pragma omp parallel for
    //     for (i = 0; i < nMols; i++) {

    //         int j;
    //         // number totalMass = 0.0;
    //         number totalMass = masses[i];

    //         comX[i] = 0.0;
    //         comY[i] = 0.0;
    //         comZ[i] = 0.0;

    //         tMinusDtoverTwo->comVx[i] = 0.0;
    //         tMinusDtoverTwo->comVy[i] = 0.0;
    //         tMinusDtoverTwo->comVz[i] = 0.0;

    //         tMinusDtoverTwo->ekinX[i] = 0.0;
    //         tMinusDtoverTwo->ekinY[i] = 0.0;
    //         tMinusDtoverTwo->ekinZ[i] = 0.0;

    //         for (j=groupsBegin[i]; j < groupsEnd[j]; j++) {

    //             number mass = topology.masses[i];
    //             number halfMass = 0.5 * mass;
    //             // totalMass += mass;

    //             comX[i] += conf.x[j] * mass;
    //             comY[i] += conf.y[j] * mass;
    //             comZ[i] += conf.z[j] * mass;

    //             number vx = conf.tMinusDtoverTwo->vx[j];
    //             number vy = conf.tMinusDtoverTwo->vy[j];
    //             number vz = conf.tMinusDtoverTwo->vz[j];

    //             tMinusDtoverTwo->comVx[i] += vx * mass;
    //             tMinusDtoverTwo->comVy[i] += vy * mass;
    //             tMinusDtoverTwo->comVz[i] += vz * mass;

    //             // kinetic energy tensor of molecule i at time t-dt/2
    //             tMinusDtoverTwo->ekinX[i] += vx * vx * halfMass;
    //             tMinusDtoverTwo->ekinY[i] += vy * vy * halfMass;
    //             tMinusDtoverTwo->ekinZ[i] += vz * vz * halfMass;

    //         }
    //         // COM of position molecule i
    //         comX[i] /= totalMass;
    //         comY[i] /= totalMass;
    //         comZ[i] /= totalMass;

    //         // COM of velocities of molecule i
    //         tMinusDtoverTwo->comVx[i] /= totalMass;
    //         tMinusDtoverTwo->comVy[i] /= totalMass;
    //         tMinusDtoverTwo->comVz[i] /= totalMass;

    //         number halfTotalMass = totalMass * 0.5;
    //         // COM kinetic energy of molecule i
    //         tMinusDtoverTwo->comEkinX[i] = tMinusDtoverTwo->comVx[i] * tMinusDtoverTwo->comVx[i] * halfTotalMass;
    //         tMinusDtoverTwo->comEkinY[i] = tMinusDtoverTwo->comVy[i] * tMinusDtoverTwo->comVy[i] * halfTotalMass;
    //         tMinusDtoverTwo->comEkinZ[i] = tMinusDtoverTwo->comVz[i] * tMinusDtoverTwo->comVz[i] * halfTotalMass;

    //         // total kinetic energy
    //         tMinusDtoverTwo->ekin[i] += tMinusDtoverTwo->ekinX[i] + tMinusDtoverTwo->ekinY[i] + tMinusDtoverTwo->ekinZ[i];

    //         // centre of mass kinetic energy
    //         tMinusDtoverTwo->comEkin[i] += tMinusDtoverTwo->comEkinX[i] + tMinusDtoverTwo->comEkinY[i] + tMinusDtoverTwo->comEkinZ[i];

    //         // rotational + internal kinetic energy
    //         tMinusDtoverTwo->rotIntKin[i] = tMinusDtoverTwo->ekin[i] - tMinusDtoverTwo->comEkin[i];

    //         // Now gatherring based on COM
    //         for (j=groupsBegin[i]; j < groupsEnd[j]; j++) {

    //         }
    //     }

    // }


  //   VHR_INLINE
  //   void gatherConfIntoBoxAndRemoveCOMMotion(Configuration &conf, const Topology &topology, const Box &simBox)
  //   {
  //    int i;

        // number ekin = 0.0;

  //    // loop over molecules
  //       #pragma omp parallel for
  //       for (i = 0; i < nMols; i++) {

  //        int j;
  //           // number totalMass = 0.0;
  //           number totalMass = masses[i];

  //           comX[i] = 0.0;
  //           comY[i] = 0.0;
  //           comZ[i] = 0.0;

        //  tMinusDtoverTwo->comVx[i] = 0.0;
        //  tMinusDtoverTwo->comVy[i] = 0.0;
        //  tMinusDtoverTwo->comVz[i] = 0.0;

        //  tMinusDtoverTwo->ekinX[i] = 0.0;
        //  tMinusDtoverTwo->ekinY[i] = 0.0;
        //  tMinusDtoverTwo->ekinZ[i] = 0.0;

     //        for (j=groupsBegin[i]; j < groupsEnd[j]; j++) {

  //            number mass = topology.masses[i];
  //            number halfMass = 0.5 * mass;
  //            // totalMass += mass;

  //            comX[i] += conf.x[j] * mass;
  //            comY[i] += conf.y[j] * mass;
  //            comZ[i] += conf.z[j] * mass;

  //            number vx = conf.tMinusDtoverTwo->vx[j];
  //            number vy = conf.tMinusDtoverTwo->vy[j];
  //            number vz = conf.tMinusDtoverTwo->vz[j];

  //            tMinusDtoverTwo->comVx[i] += vx * mass;
  //            tMinusDtoverTwo->comVy[i] += vy * mass;
  //            tMinusDtoverTwo->comVz[i] += vz * mass;

        //      // kinetic energy tensor of molecule i at time t-dt/2
  //            tMinusDtoverTwo->ekinX[i] += vx * vx * halfMass;
  //            tMinusDtoverTwo->ekinY[i] += vy * vy * halfMass;
  //            tMinusDtoverTwo->ekinZ[i] += vz * vz * halfMass;

  //        }
  //        // COM of position molecule i
  //        comX[i] /= totalMass;
  //        comY[i] /= totalMass;
  //        comZ[i] /= totalMass;

        //  // COM of velocities of molecule i
        //  tMinusDtoverTwo->comVx[i] /= totalMass;
        //  tMinusDtoverTwo->comVy[i] /= totalMass;
        //  tMinusDtoverTwo->comVz[i] /= totalMass;

  //           number halfTotalMass = totalMass * 0.5;
        //  // COM kinetic energy of molecule i
        //  tMinusDtoverTwo->comEkinX[i] = tMinusDtoverTwo->comVx[i] * tMinusDtoverTwo->comVx[i] * halfTotalMass;
        //  tMinusDtoverTwo->comEkinY[i] = tMinusDtoverTwo->comVy[i] * tMinusDtoverTwo->comVy[i] * halfTotalMass;
        //  tMinusDtoverTwo->comEkinZ[i] = tMinusDtoverTwo->comVz[i] * tMinusDtoverTwo->comVz[i] * halfTotalMass;

        //  // total kinetic energy
        //  tMinusDtoverTwo->ekin[i] += tMinusDtoverTwo->ekinX[i] + tMinusDtoverTwo->ekinY[i] + tMinusDtoverTwo->ekinZ[i];

        //  // centre of mass kinetic energy
        //  tMinusDtoverTwo->comEkin[i] += tMinusDtoverTwo->comEkinX[i] + tMinusDtoverTwo->comEkinY[i] + tMinusDtoverTwo->comEkinZ[i];

        //  // rotational + internal kinetic energy
        //  tMinusDtoverTwo->rotIntKin[i] = tMinusDtoverTwo->ekin[i] - tMinusDtoverTwo->comEkin[i];

        //  // Now gatherring based on COM
        //  for (j=groupsBegin[i]; j < groupsEnd[j]; j++) {

        //  }
  //       }

  //   }

  //   void computeCOMandKineticEnergies(const Configuration &conf, const Topology &topology)
  //   {
  //    int i;

        // number ekin = 0.0;
  //    // loop over molecules
  //       #pragma omp parallel for reduction(+:ekin) reduction(+:totalMass)
  //       for (i=0; i < nMols; i++) {

  //        int j;
  //           number totalMass = 0.0;

  //           current->comX = 0.0;
  //           current->comY = 0.0;
  //           current->comZ = 0.0;

        //  current->comVx[i] = 0.0;
        //  current->comVy[i] = 0.0;
        //  current->comVz[i] = 0.0;

        //  current->ekinX[i] = 0.0;
        //  current->ekinY[i] = 0.0;
        //  current->ekinZ[i] = 0.0;

     //        for (j=groupsBegin[i]; j < groupsEnd[j]; j++) {

  //            number mass = topology.masses[i];
  //            number halfMass = 0.5 * mass;
  //            totalMass += mass;

  //            // These are all new values
  //            current->comX[i] += conf.x[j] * mass;
  //            current->comY[i] += conf.y[j] * mass;
  //            current->comZ[i] += conf.z[j] * mass;

  //            current->comVx[i] += conf.current->Vx[j] * mass;
  //            current->comVy[i] += conf.current->Vy[j] * mass;
  //            current->comVz[i] += conf.current->Vz[j] * mass;

        //      // kinetic energy tensor of molecule i at time T+dT/2
  //            current->ekinX[i] += conf.current->Vx[j] * conf.current->Vx[j] * halfMass;
  //            current->ekinY[i] += conf.current->Vy[j] * conf.current->Vy[j] * halfMass;
  //            current->ekinZ[i] += conf.current->Vz[j] * conf.current->Vz[j] * halfMass;

  //            // Average between T-dT/2 and T+dT/2
  //            // infoAtT.comVx[i] += (conf.current->Vx[j] + conf.old->Vx[j]) * halfMass;
  //            // infoAtT.comVy[i] += (conf.current->Vy[j] + conf.old->Vy[j]) * halfMass;
  //            // infoAtT.comVz[i] += (conf.current->Vz[j] + conf.old->Vz[j]) * halfMass;

  //            // infoAtT.>ekinX[i] += conf.current->Vx[j] * conf.current->Vx[j] * halfMass;
  //            // infoAtT.>ekinY[i] += conf.current->Vy[j] * conf.current->Vy[j] * halfMass;
  //            // infoAtT.>ekinZ[i] += conf.current->Vz[j] * conf.current->Vz[j] * halfMass;


  //        }
  //        // COM of position molecule i
  //        current->comX[i] /= totalMass;
  //        current->comY[i] /= totalMass;
  //        current->comZ[i] /= totalMass;

        //  // COM of velocities of molecule i
        //  current->comVx[i] /= totalMass;
        //  current->comVy[i] /= totalMass;
        //  current->comVz[i] /= totalMass;

  //           number halfTotalMass = totalMass * 0.5;
        //  // COM kinetic energy of molecule i
        //  current->comEkinX[i] = current->comVx[i] * current->comVx[i] * halfTotalMass;
        //  current->comEkinY[i] = current->comVy[i] * current->comVy[i] * halfTotalMass;
        //  current->comEkinZ[i] = current->comVz[i] * current->comVz[i] * halfTotalMass;

        //  current->ekin[i] += current->ekinX[i] + current->ekinY[i] + current->ekinZ[i];
        //  current->comEkin[i] += current->comEkinX[i] + current->comEkinY[i] + current->comEkinZ[i];

        //  current->rotIntKin[i] = current->ekin[i] - current->comEkin[i];

  //        // Average between T-dT/2 and T+dT/2
  //        // infoAtT.comVx[i] /= totalMass;
  //        // infoAtT.comVy[i] /= totalMass;
  //        // infoAtT.comVz[i] /= totalMass;

        //  // COM of velocities of molecule i at time T
  //        infoAtT.comVx[i] = 0.5 * (current->comVx[i] + old->comVx[i]);
        //  infoAtT.comVy[i] = 0.5 * (current->comVy[i] + old->comVy[i]);
        //  infoAtT.comVz[i] = 0.5 * (current->comVz[i] + old->comVz[i]);

        //  // COM kinetic energy of molecule i at time T
  //        infoAtT.comEkinX[i] = 0.5 * (current->comEkinX[i] + old->comEkinX[i]);
        //  infoAtT.comEkinY[i] = 0.5 * (current->comEkinY[i] + old->comEkinY[i]);
        //  infoAtT.comEkinZ[i] = 0.5 * (current->comEkinZ[i] + old->comEkinZ[i]);

        //  // kinetic energy tensor of molecule i
        //  infoAtT.>ekinX[i] = 0.5 * (current->ekinX[i] + old->ekinX[i]);
        //  infoAtT.>ekinY[i] = 0.5 * (current->ekinY[i] + old->ekinY[i]);
        //  infoAtT.>ekinZ[i] = 0.5 * (current->ekinZ[i] + old->ekinZ[i]);

        //  current->ekin[i] += current->ekinX[i] + current->ekinY[i] + current->ekinZ[i];
        //  current->comEkin[i] += current->comEkinX[i] + current->comEkinY[i] + current->comEkinZ[i];

        //  current->rotIntKin[i] = current->ekin[i] - current->comEkin[i];


        //  infoAtT.rotIntKin[i] = 0.5 * (current->rotIntKin[i] + current->rotIntKin[i]);
  //       }

  //   }


};


#endif //MOLECULE_HPP