#ifndef BONDED_CPU_HPP
#define BONDED_CPU_HPP

#include "../../main.hpp"
#include "../../domain/domain.hpp"

#include "../../configuration/conf.hpp"
#include "../../simulation/box.hpp"
#include "../../simulation/simparam.hpp"
#include "../../topology/topology.hpp"

#include "../constraints/shake.hpp"

VHR_ALWAYS_INLINE
number g96BondEnergyForce(const RVec &rij, const number kq, const number b0, RVec &force)
{
    number dx2 = (trans(rij) * rij) - (b0 * b0);

    force = (- kq * dx2) * rij;

    return 0.25 * kq * dx2 * dx2;
};

VHR_ALWAYS_INLINE
number g96BondAngleEnergyForce(const RVec &rij, const RVec &rkj, const number kq, const number cos0, RVec &fi, RVec &fk)
{
    const number invNormrij = 1.0 / (std::sqrt(trans(rij) * rij));
    const number invNormrkj = 1.0 / (std::sqrt(trans(rkj) * rkj));

    const number costheta = (trans(rij) * rkj) * invNormrij * invNormrkj;

    const number dcos = costheta - cos0;


    fi = (- kq * dcos * invNormrij) * ( (rkj * invNormrkj) - (rij * invNormrij * costheta) );
    fk = (- kq * dcos * invNormrkj) * ( (rij * invNormrij) - (rkj * invNormrkj * costheta) );

    return 0.5 * kq * dcos * dcos;
};

// VHR_ALWAYS_INLINE
// number g96ImproperDihedralEnergyForce(const RVec &rij, const RVec &rkj, const number kq, const number cos0, RVec &fi, RVec &fk)
// {
//     // bondAngleEnergy += g96BondAngleEnergyForce(rij, rkj, kq, cos0, fi, fk);
// // return (trans(u) * v) / (std::sqrt((trans(u) * u)) * std::sqrt((trans(v) * v)));

//     const number invNormrij = 1.0 / (std::sqrt(trans(rij) * rij));
//     const number invNormrkj = 1.0 / (std::sqrt(trans(rkj) * rkj));

//     const number costheta = (trans(rij) * rkj) * invNormrij * invNormrkj;

//     const number dcos = costheta - cos0;


//     fi = (- kq * dcos * invNormrij) * ( (rkj * invNormrkj) - (rij * invNormrij * costheta) );
//     fk = (- kq * dcos * invNormrkj) * ( (rij * invNormrij) - (rkj * invNormrkj * costheta) );

//     return 0.5 * kq * dcos * dcos;
// };


VHR_ALWAYS_INLINE
number g96ProperDihedralEnergyForce(const RVec &rij, const RVec &rkj, const RVec &rkl, const number cp, const number pd, const int np, RVec &fi, RVec &fj, RVec &fk, RVec &fl)
{
    const RVec rmj = rij % rkj;
    const RVec rnk = rkj % rkl;

    const number invRkjSqr = 1.0 / (trans(rkj) * rkj);

    number frim = trans(rij) * rkj * invRkjSqr;
    number frln = trans(rkl) * rkj * invRkjSqr;

    const RVec rim = (rij - RVec(frim)) * rkj;
    const RVec rln = (frln * rkj) - rkl;

    const number invDim = std::sqrt(trans(rim) * rim);
    const number invDln = std::sqrt(trans(rln) * rln);

    const number cosphi = (trans(rim) * rln) * invDim * invDln;

    const number cosphi2 = cosphi * cosphi;
    const number cosphi3 = cosphi2 * cosphi;
    const number cosphi4 = cosphi3 * cosphi;

    number cosmphi = 0.0;
    number dcosmphi = 0.0;

    switch (np) {
      case 0:
        cosmphi = 0.0;
        dcosmphi = 0.0;
        break;
      case 1:
        cosmphi = cosphi;
        dcosmphi = 1.0;
        break;
      case 2:
        cosmphi = 2.0 * cosphi2 - 1.0;
        dcosmphi = 4.0 * cosphi;
        break;
      case 3:
        cosmphi = 4.0 * cosphi3 - 3.0 * cosphi;
        dcosmphi = 12.0 * cosphi2 - 3.0;
        break;
      case 4:
        cosmphi = 8.0 * cosphi4 - 8.0 * cosphi2 + 1.0;
        dcosmphi = 32.0 * cosphi3 - 16.0 * cosphi;
        break;
      case 5:
        cosmphi = 16.0 * cosphi4 * cosphi - 20.0 * cosphi3 + 5.0 * cosphi;
        dcosmphi = 80.0 * cosphi4 - 60.0 * cosphi2 + 5.0;
        break;
      case 6:
        cosmphi = 32.0 * cosphi4 * cosphi2 - 48.0 * cosphi4 + 18.0 * cosphi2 - 1.0;
        dcosmphi = 192.0 * cosphi4 * cosphi - 192.0 * cosphi3 + 36.0 * cosphi;
        break;
    }

    number ki = -cp * pd * dcosmphi * invDim;
    number kl = -cp * pd * dcosmphi * invDln;
    number kj1 = frim - 1.0;
    number kj2 = frln;

    fi = ki * (rln * invDln - rim * invDim * cosphi);
    fl = kl * (rim * invDim - rln * invDln * cosphi);
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0 * (fi + fj + fl);

    return cp * (1 + pd * cosmphi);;
};

struct Bonded {

    number energy;
    number bondEnergy;
    number bondAngleEnergy;
    number properDihedralEnergy;

    // TODO
    // parallelize this function
    void compute(Configuration &conf, const Topology &topology, const SimParameters &simParam, Box &simBox, const bool computeVirial)
    {
        startTimer(bondTimer);

        const int nAtoms = conf.nAtoms;
        const int nBonds = topology.bonds.nBonds;
        const int nBondAngles = topology.bondAngles.nBondAngles;
        const int nProperDihedrals = topology.properDihedrals.nProperDihedrals;

        number bondEnergy = 0.0;
        number bondAngleEnergy = 0.0;
        number properDihedralEnergy = 0.0;

        TVector virialTensors;

        virialTensors.resize(dd.nprocs);

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:bondEnergy, bondAngleEnergy, properDihedralEnergy)
    #endif
        {
            int i, tid, fromidx, toidx;


    #if defined(_OPENMP)
            tid = omp_get_thread_num();
    #else
            tid = 0;
    #endif

            fromidx = tid * nAtoms;
            toidx = fromidx + nAtoms;

            RVec zero(0.0);

            // zeroing forces for all but the first tid
            // this is due to the fact that the force vector actually
            // contains the nonbonded forces
            if(tid > 0) {
                for(i = fromidx; i < toidx; i++)
                    conf.f[i] = zero;
            }

            // zeroing the virial tensor
            virialTensors[tid] = 0.0;

            // Bonds
            if(!simParam.shake) {

                // computing the bonds
                // #pragma omp parallel for reduction(+: bondEnergy)
                // for(i = 0; i < nBonds; i++) {
                for(i = 0; i < nBonds; i += dd.nprocs) {

                    // the actual loop index is j
                    int j = i + tid;
                    if (j >= nBonds)
                        break;

                    RVec force;

                    const int ii = topology.bonds.atoms[2*j];
                    const int jj = topology.bonds.atoms[(2*j) + 1];

                    const number kq = topology.bonds.kq[j];
                    const number b0 = topology.bonds.b0[j];

                    // TODO
                    // Are we going to keep the pbc here??
                    const RVec rij = simBox.pbc(conf.current->x[ii] - conf.current->x[jj]);

                    bondEnergy += g96BondEnergyForce(rij, kq, b0, force);

                    conf.f[ii + fromidx] += force;
                    conf.f[jj + fromidx] -= force;

                    // virial contribution
                    if(computeVirial) {
                        virialTensors[tid] += rij * trans(force);
                    }
                }
            }

            // BondAngles
            // #pragma omp parallel for reduction(+: bondAngleEnergy)
            for(i = 0; i < nBondAngles; i += dd.nprocs) {

                // the actual loop index is j
                int j = i + tid;
                if (j >= nBondAngles)
                    break;

                RVec fi, fk;

                const int ii = topology.bondAngles.atoms[3*j];
                const int jj = topology.bondAngles.atoms[(3*j) + 1];
                const int kk = topology.bondAngles.atoms[(3*j) + 2];

                const number kq = topology.bondAngles.kq[j];
                const number cos0 = topology.bondAngles.cos0[j];

                const RVec rij = simBox.pbc(conf.current->x[ii] - conf.current->x[jj]);
                const RVec rkj = simBox.pbc(conf.current->x[kk] - conf.current->x[jj]);


                bondAngleEnergy += g96BondAngleEnergyForce(rij, rkj, kq, cos0, fi, fk);

                conf.f[ii + fromidx] += fi;
                conf.f[jj + fromidx] -= fi + fk;
                conf.f[kk + fromidx] += fk;

                // virial contribution
                if(computeVirial) {
                    virialTensors[tid] += (rij * trans(fi)) + (rkj * trans(fk));
                }

            }

            // Proper dihedrals
            // #pragma omp parallel for reduction(+: properDihedralEnergy)
            for(i = 0; i < nProperDihedrals; i += dd.nprocs) {

                // the actual loop index is j
                int j = i + tid;
                if (j >= nProperDihedrals)
                    break;

                RVec fi, fj, fk, fl;

                const int ii = topology.properDihedrals.atoms[4*j];
                const int jj = topology.properDihedrals.atoms[(4*j) + 1];
                const int kk = topology.properDihedrals.atoms[(4*j) + 2];
                const int ll = topology.properDihedrals.atoms[(4*j) + 3];

                const number cp = topology.properDihedrals.cp[j];
                const number pd = topology.properDihedrals.pd[j];
                const int np = topology.properDihedrals.np[j];

                const RVec rij = simBox.pbc(conf.current->x[ii] - conf.current->x[jj]);
                const RVec rkj = simBox.pbc(conf.current->x[kk] - conf.current->x[jj]);
                const RVec rkl = simBox.pbc(conf.current->x[kk] - conf.current->x[ll]);
                const RVec rlj = simBox.pbc(conf.current->x[ll] - conf.current->x[jj]);

                properDihedralEnergy += g96ProperDihedralEnergyForce(rij, rkj, rkl, cp, pd, np, fi, fj, fk, fl);

                conf.f[ii + fromidx] += fi;
                conf.f[jj + fromidx] += fj;
                conf.f[kk + fromidx] += fk;
                conf.f[ll + fromidx] += fl;

                // virial contribution
                if(computeVirial) {
                    virialTensors[tid] += (rij * trans(fi)) + (rkj * trans(fk)) + (rlj * trans(fl));
                }
            }

            // before reducing the forces, we have to make sure
            // that all threads are done adding to them.
            #if defined (_OPENMP)
            #pragma omp barrier

            // set equal chunks of index ranges
            i = 1 + (nAtoms / dd.nprocs);
            fromidx = tid * i;
            toidx = fromidx + i;
            if (toidx > nAtoms)
                toidx = nAtoms;

            // reduce forces from threads with tid != 0 into
            // the storage of the first thread. since we have
            // threads already spawned, we do this in parallel.
            for (i = 1; i < dd.nprocs; ++i) {
                int offsets, j;

                offsets = i * nAtoms;

                for (j = fromidx; j < toidx; j++) {
                    conf.f[j] += conf.f[offsets + j];
                }
            }

            if(tid == 0) {
                // cerr.printf("reducing virialTensors\n");
                // cerr.flush();
                for (i = 1; i < dd.nprocs; ++i) {
                    virialTensors[0] += virialTensors[i];
                }
            }
            #endif


        }
        this->energy = bondEnergy + bondAngleEnergy + properDihedralEnergy;

        this->bondEnergy = bondEnergy;
        this->bondAngleEnergy = bondAngleEnergy;
        this->properDihedralEnergy = properDihedralEnergy;

        conf.current->virialTensor += virialTensors[0];


        addToTime(bondTimer, bondTime);
    }

};


#endif //BONDED_CPU_HPP
