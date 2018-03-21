#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include <set>
#include <map>

#include "../main.hpp"            // For function prototypes and run-time constants.
#include "topologyfields.hpp"  // For function prototypes and run-time constants.
#include "molecules.hpp"

// #define addexclusion(e,i,j)   (e)[((int) (i))] |= (1<<((int) (j)))
// #define remveexclusion(e,i,j) (e)[((int) (i))] &= (~(1<<((int) (j))))
// #define isexcluded(e,i,j)     (bool) (j-i-1 > 31 ? 0 : (e)[((int) (i))] & (1<<((int) (j))))
// #define notexcluded(e,i,j)   (!(isexcluded(e,i,j)))

#define addexclusion(e,i,j)   (e)[((int) (i))] |= (1<<((int) (j)))
#define removeexclusion(e,i,j) (e)[((int) (i))] &= (~(1<<((int) (j))))
#define isexcluded(e,i,j)     (bool) (j-i > 30 ? 0 : (e)[((int) (i))] & (1<<((int) (j))))
#define notexcluded(e,i,j)   !(isexcluded((e),(i),(j)))

// #define add14(e,i,j)   (e)[((int) (i))] |= (1<<((int) (j)))
// #define remove14(e,i,j) (e)[((int) (i))] &= (~(1<<((int) (j))))
// #define has14(e,i,j)     (bool) (j-i > 30 ? 0 : (e)[((int) (i))] & (1<<((int) (j))))
// #define doesnothave14(e,i,j)   !(isexcluded((e),(i),(j)))

struct Topology {

    std::string title;

    int nResidues;
    std::vector<std::string> resNames;

    int nAtoms;
    std::vector<std::string> atomNames;

    NVector masses;
    NVector invMasses;
    NVector charges;

    std::vector<int> resIds;
    std::vector<int> excl;
    std::vector<std::vector<int>> onefours;

    int nAtomTypes;
    std::vector<std::string> atomTypeNames;
    std::vector<int> atomTypes;

    int nLJParameters;
    std::vector<LJParameters> ljParameters;
    // std::map<int,LJParameters> usedLJParameters;

    BondStretchTypes bondStretchTypes;
    Bonds bonds;

    BondAngleTypes bondAngleTypes;
    BondAngles bondAngles;


    ImproperDihedralTypes improperDihedralTypes;

    ProperDihedralTypes properDihedralTypes;
    ProperDihedrals properDihedrals;

    // int nSoluteMolecules;
    // std::vector<int> soluteMolecules;
    int nMolecules;
    Molecules molecules;

    std::vector<int> constrainedAtoms;

    // vector charges;

    // int nBondStretchTypes;
    // std::vector<NonBondedParameters> nonBondedTypes;

    void reserve(int nAtoms, int nmasses, int nLJParameters)
    {
        this->atomNames.reserve(nAtoms);
        this->charges.reserve(nAtoms);

        this->masses.reserve(nmasses);
        this->invMasses.reserve(nmasses);
        this->atomTypes.reserve(nLJParameters);
        this->ljParameters.reserve(nLJParameters);
    }

    VHR_ALWAYS_INLINE
    void populateConstrainedAtoms()
    {
        int i, nbonds = this->bonds.nBonds;
        for(i = 0; i < nbonds; i++) {
            const int ii = this->bonds.atoms[2*i];
            const int jj = this->bonds.atoms[(2*i) + 1];

            constrainedAtoms.push_back(ii);
            constrainedAtoms.push_back(jj);
        }

        // Removing duplicated
        std::sort(constrainedAtoms.begin(), constrainedAtoms.end());
        constrainedAtoms.erase(std::unique(constrainedAtoms.begin(), constrainedAtoms.end()), constrainedAtoms.end());


    }

};


#endif //TOPOLOGY_HPP