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
    void compute(Configuration &conf, const Topology &topology, const SimParameters &simParam, Box &simBox)
    {
        startTimer(bondTimer);

        int i;
        int nBonds = topology.bonds.nBonds;
        int nBondAngles = topology.bondAngles.nBondAngles;
        int nProperDihedrals = topology.properDihedrals.nProperDihedrals;

        number bondEnergy = 0.0;
        number bondAngleEnergy = 0.0;
        number properDihedralEnergy = 0.0;


        // Bonds
        if(!simParam.shake) {

            // computing the bonds
            // #pragma omp parallel for reduction(+: bondEnergy)
            for(i = 0; i < nBonds; i++) {

                RVec force;

                const int ii = topology.bonds.atoms[2*i];
                const int jj = topology.bonds.atoms[(2*i) + 1];

                const number kq = topology.bonds.kq[i];
                const number b0 = topology.bonds.b0[i];

                // TODO
                // Are we going to keep the pbc here??
                const RVec rij = simBox.pbc(conf.current->x[ii] - conf.current->x[jj]);

                bondEnergy += g96BondEnergyForce(rij, kq, b0, force);

                conf.f[ii] += force;
                conf.f[jj] -= force;

                // virial contribution
                conf.current->virialTensor += rij * trans(force);

            }
        }

        // BondAngles
        // #pragma omp parallel for reduction(+: bondAngleEnergy)
        for(i = 0; i < nBondAngles; i++) {

            RVec fi, fk;

            const int ii = topology.bondAngles.atoms[3*i];
            const int jj = topology.bondAngles.atoms[(3*i) + 1];
            const int kk = topology.bondAngles.atoms[(3*i) + 2];

            const RVec rij = simBox.pbc(conf.current->x[ii] - conf.current->x[jj]);
            const RVec rkj = simBox.pbc(conf.current->x[kk] - conf.current->x[jj]);

            const number kq = topology.bondAngles.kq[i];
            const number cos0 = topology.bondAngles.cos0[i];


            bondAngleEnergy += g96BondAngleEnergyForce(rij, rkj, kq, cos0, fi, fk);

            conf.f[ii] += fi;
            conf.f[jj] -= fi + fk;
            conf.f[kk] += fk;

            // virial contribution
            conf.current->virialTensor += (rij * trans(fi)) + (rkj * trans(fk));

        }

        // Proper dihedrals
        // #pragma omp parallel for reduction(+: properDihedralEnergy)
        for(i = 0; i < nProperDihedrals; i++) {

            RVec fi, fj, fk, fl;

            const int ii = topology.properDihedrals.atoms[3*i];
            const int jj = topology.properDihedrals.atoms[(3*i) + 1];
            const int kk = topology.properDihedrals.atoms[(3*i) + 2];
            const int ll = topology.properDihedrals.atoms[(3*i) + 3];

            const RVec rij = simBox.pbc(conf.current->x[ii] - conf.current->x[jj]);
            const RVec rkj = simBox.pbc(conf.current->x[kk] - conf.current->x[jj]);
            const RVec rkl = simBox.pbc(conf.current->x[kk] - conf.current->x[ll]);
            const RVec rlj = simBox.pbc(conf.current->x[ll] - conf.current->x[jj]);

            const number cp = topology.properDihedrals.cp[i];
            const number pd = topology.properDihedrals.pd[i];
            const int np = topology.properDihedrals.np[i];

            properDihedralEnergy += g96ProperDihedralEnergyForce(rij, rkj, rkl, cp, pd, np, fi, fj, fk, fl);

            conf.f[ii] += fi;
            conf.f[jj] += fj;
            conf.f[kk] += fk;
            conf.f[ll] += fl;

            conf.current->virialTensor += (rij * trans(fi)) + (rkj * trans(fk)) + (rlj * trans(fl));

        }

        this->energy = bondEnergy + bondAngleEnergy + properDihedralEnergy;

        this->bondEnergy = bondEnergy;
        this->bondAngleEnergy = bondAngleEnergy;
        this->properDihedralEnergy = properDihedralEnergy;

        addToTime(bondTimer, bondTime);
    }

};


#endif //BONDED_CPU_HPP
