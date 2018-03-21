#ifndef TOPOLOGYFIELDS_HPP
#define TOPOLOGYFIELDS_HPP

#include "../main.hpp"  // For function prototypes and run-time constants.


struct BondStretchTypes
{
    int ntypes;
    std::vector<number> kh; // harmonic force constant
    std::vector<number> kq; // quartic force constant
    std::vector<number> b0; // equilibrium distance
};

struct Bonds
{
    int nBonds;
    std::vector<int> atoms; // atoms involved in the bond
    std::vector<int> type; // bond type
    std::vector<number> kq; // quartic force constant
    std::vector<number> b0; // equilibrium distance
};

struct BondAngleTypes
{
    int ntypes;
    std::vector<number> kh; // harmonic force constant
    std::vector<number> kq; // cosine force constant
    std::vector<number> a0; // equilibrium angle
};

struct BondAngles
{
    int nBondAngles;
    std::vector<int> atoms; // atoms involved in the bond
    std::vector<int> type; // angle type
    std::vector<number> kq; // force constant
    std::vector<number> cos0; // equilibrium cosine!! We do not same the equilibrium angle!!!
};


struct ImproperDihedralTypes
{
    int ntypes;
    std::vector<number> kh; // harmonic force constant
    std::vector<number> a0; // equilibrium angle
};

struct ProperDihedralTypes
{
    int ntypes;
    std::vector<number> cp; // force constant
    std::vector<number> pd; // phase-shift
    std::vector<number> np; // multiplicity
};

struct ProperDihedrals
{
    int nProperDihedrals;
    std::vector<int> atoms; // atoms involved in the bond
    std::vector<int> type; // angle type
    std::vector<number> cp; // force constant
    std::vector<number> pd; // equilibrium cosine!! We do not same the equilibrium angle!!!
    std::vector<int> np; // equilibrium cosine!! We do not same the equilibrium angle!!!
};


struct LJParameters {

    int ntypes;
    int iac, jac;
    number c12, c6, c12_14, c6_14;

};

#endif //TOPOLOGYFIELDS_HPP