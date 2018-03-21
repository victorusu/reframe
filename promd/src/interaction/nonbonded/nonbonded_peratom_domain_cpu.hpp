#ifndef NONBONDED_PERATOM_DOMAIN_CPU_HPP
#define NONBONDED_PERATOM_DOMAIN_CPU_HPP

#include "../../main.hpp"
#include "../../domain/domain.hpp"

#include "../pairlist/pairlist.hpp"

#include "../../configuration/conf.hpp"
#include "../../simulation/box.hpp"
#include "../../simulation/simparam.hpp"
#include "../../simulation/barostat/barostat.hpp"

#include "../../topology/topology.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif

// number pbc(number x, const number boxby2, const number box)
// {
//     while (x >  boxby2) x -= box;
//     while (x < -boxby2) x += box;
//     return x;
// }

struct Nonbonded {

    number energy;

    number vdwEnergy;
    number coulEnergy;
    number vdw14Energy;
    number coul14Energy;

    // VHR_ALWAYS_INLINE
    // static void azzero(number *d, const int n)
    // {
    //     int i;
    //     for (i=0; i<n; ++i) {
    //         d[i]=0.0;
    //     }
    // }

    void compute(Configuration &conf, const Topology &topology, const PairList &pairlist, const SimParameters &simParam, Box &simBox, const bool computeVirial)
    {
        startTimer(nonbondedComputeTimer);

        number eVdWTotal = 0.0;
        number e14VdwTotal = 0.0;

        number eCoulTotal = 0.0;
        number e14CoulTotal = 0.0;

        TVector virialTensors;

        virialTensors.resize(dd.nprocs);

        // specialized cut-offs
        const number rcutVdwSqr = simParam.rcutvdw * simParam.rcutvdw;
        const number rcutCoulSqr = simParam.rcutcoul * simParam.rcutcoul;
        const number rcutVdWAndCoulSqr = std::max(rcutVdwSqr, rcutCoulSqr);

        // reaction-field term
        // const crf = 2.0 * simParam.


        const int nAtoms = conf.nAtoms;

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:eVdWTotal) reduction(+:e14VdwTotal) reduction(+:eCoulTotal) reduction(+:e14CoulTotal)
    #endif
        {
            number c12, c6;

            int i, tid, fromidx, toidx;

            // TODO
            // precompute some constants

            // let each thread operate on a different
            // part of the enlarged force array
    #if defined(_OPENMP)
            tid = omp_get_thread_num();
    #else
            tid = 0;
    #endif

            fromidx = tid * nAtoms;
            toidx = fromidx + nAtoms;

            RVec zero(0.0);

            // zeroing forces
            for(i = fromidx; i < toidx; i++)
                conf.f[i] = zero;

            // zeroing the virial tensor
            virialTensors[tid] = 0.0;

            if(tid == 0) {
                startTimer(nonbondedFirstTimer);
            }

            // self interaction of atoms in cell
            for(i = 0; i < nAtoms; i += dd.nprocs) {

                int j = i + tid;
                if (j >= nAtoms)
                    break;

                int ii = j;

                const RVec rii = conf.current->x[ii];

                // vdW pair nonbonded calculation
                {
                    const NeighList *neighList = pairlist.vdWNeigh.data() + j;

                    for (j = 0; j < neighList->nAtoms; j++) {

                        int jj = neighList->listIds[j];

                        number rsq;

                        // get distance between particle i and j
                        const RVec r = simBox.pbc(rii - conf.current->x[jj]);

                        // rsq = rx*rx + ry*ry + rz*rz;
                        // inner product r * r;
                        rsq = trans(r) * r;

                        // compute force and energy if within cutoff
                        if (rsq < rcutVdwSqr) {

                            number r6inv, rinvSqr, ffac;
                            int lower, max;

                            rinvSqr = 1.0 / rsq;
                            r6inv = rinvSqr * rinvSqr * rinvSqr;

                            lower = topology.atomTypes[ii];
                            max = topology.atomTypes[jj];

                            if(lower > max) {
                                std::swap(lower,max);
                                // fprintf(stderr, "%d (lower) > %d (max) from atoms %d and %d! this should not happend\n", lower, max, ii, jj);
                            }

                            int iac = ((max+1)*max)/2 + lower;
                            LJParameters lj = topology.ljParameters[iac];
                            c12 = lj.c12;
                            c6 = lj.c6;

                            ffac = (12.0 * c12 * r6inv - 6.0 * c6) *r6inv * rinvSqr;
                            eVdWTotal += r6inv * (c12 * r6inv - c6);

                            RVec force = r * ffac;

                            conf.f[ii + fromidx] += force;
                            conf.f[jj + fromidx] -= force;

                            if(computeVirial) {
                                // virial calculation   r  * f
                                // outer product of r and f
                                // cerr.printf("computing virialTensors[%d]\n", tid);
                                virialTensors[tid] += r * trans(force);
                            }
                        }
                    }
                }

                // coul pair nonbonded calculation
                {
                    // const NeighList *neighList = pairlist.coulNeigh.data() + j;

                    // for (j = 0; j < neighList->nAtoms; j++) {

                    //     int jj = neighList->listIds[j];

                    //     number rsq;

                    //     // get distance between particle i and j
                    //     const RVec r = simBox.pbc(rii - conf.current->x[jj]);

                    //     // rsq = rx*rx + ry*ry + rz*rz;
                    //     // inner product r * r;
                    //     rsq = trans(r) * r;

                    //     // compute force and energy if within cutoff
                    //     if (rsq < rcutCoulSqr) {

                    //         number r6inv, rinvSqr, ffac;
                    //         int lower, max;

                    //         rinvSqr = 1.0 / rsq;
                    //         r6inv = rinvSqr * rinvSqr * rinvSqr;

                    //         lower = topology.atomTypes[ii];
                    //         max = topology.atomTypes[jj];

                    //         if(lower > max) {
                    //             std::swap(lower,max);
                    //             // fprintf(stderr, "%d (lower) > %d (max) from atoms %d and %d! this should not happend\n", lower, max, ii, jj);
                    //         }

                    //         ffac = (12.0 * c12 * r6inv - 6.0 * c6) *r6inv * rinvSqr;
                    //         eCoulTotal += q_eps * (disti - m_crf_2cut3i[eps] * dist2 - m_crf_cut[eps]);

                    //         RVec force = r * ffac;

                    //         conf.f[ii + fromidx] += force;
                    //         conf.f[jj + fromidx] -= force;

                    //         if(computeVirial) {
                    //             // virial calculation   r  * f
                    //             // outer product of r and f
                    //             // cerr.printf("computing virialTensors[%d]\n", tid);
                    //             virialTensors[tid] += r * trans(force);
                    //         }
                    //     }
                    // }
                }

                // TODO
                // add the coulomb interaction
                // 1-4 interaction
                {
                    const int n14s = topology.onefours[ii].size();

                    for (j = 0; j < n14s; j++) {

                        int jj = topology.onefours[ii][j];

                        number rsq;

                        // get distance between particle i and j
                        const RVec r = simBox.pbc(rii - conf.current->x[jj]);

                        // rsq = rx*rx + ry*ry + rz*rz;
                        // inner product r * r;
                        rsq = trans(r) * r;

                        // computing the 1-4 interaction
                        number r6inv, rinvSqr, ffac;
                        int lower, max;

                        rinvSqr = 1.0 / rsq;
                        r6inv = rinvSqr * rinvSqr * rinvSqr;

                        lower = topology.atomTypes[ii];
                        max = topology.atomTypes[jj];

                        if(lower > max) {
                            std::swap(lower,max);
                            // fprintf(stderr, "%d (lower) > %d (max) from atoms %d and %d! this should not happend\n", lower, max, ii, jj);
                        }

                        int iac = ((max+1)*max)/2 + lower;
                        LJParameters lj = topology.ljParameters[iac];
                        c12 = lj.c12_14;
                        c6 = lj.c6_14;

                        ffac = (12.0 * c12 * r6inv - 6.0 * c6) *r6inv * rinvSqr;
                        e14VdwTotal += r6inv * (c12 * r6inv - c6);

                        RVec force = r * ffac;

                        conf.f[ii + fromidx] += force;
                        conf.f[jj + fromidx] -= force;

                        if(computeVirial) {
                            // virial calculation   r  * f
                            // outer product of r and f
                            virialTensors[tid] += r * trans(force);
                        }
                    }
                }

            } // for(i = 0; i < nAtoms; i += dd.nprocs) {

            // before reducing the forces, we have to make sure
            // that all threads are done adding to them.
            #if defined (_OPENMP)
            #pragma omp barrier


            if(tid == 0) {
                addToTime(nonbondedFirstTimer, nonbondedFirstTime);
                startTimer(nonbondedReductionTimer);
            }

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

            if(tid == 0) {
                addToTime(nonbondedReductionTimer, nonbondedReductionTime);
            }
            #endif

        }
        this->energy = eVdWTotal  + e14VdwTotal + eCoulTotal + e14CoulTotal;

        this->vdwEnergy = eVdWTotal;
        this->coulEnergy = eCoulTotal;
        this->vdw14Energy = e14VdwTotal;
        this->coul14Energy = e14CoulTotal;

        conf.current->virialTensor = virialTensors[0];

        addToTime(nonbondedComputeTimer, nonbondedComputeTime);
    }
};

#endif //NONBONDED_PERATOM_DOMAIN_CPU_HPP
