#ifndef TOPOLOGYFILEREADER_HPP
#define TOPOLOGYFILEREADER_HPP

#include "../main.hpp"
#include "../topology/topology.hpp"

#include "blocksfile.hpp"

#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <cstdlib>
#include <exception>


class TOPOLOGYFileReader : public BlocksFile
{

public:

    TOPOLOGYFileReader();

    virtual ~TOPOLOGYFileReader();

    bool readfile(Configuration &conf, Box &simBox, Topology &topology, SimParameters &simParam)
    {
        // number num;

        std::string line;
        int i;

        // **************************************************
        // BEGIN of TITLE block
        {
            if(!findBlock("TITLE")) {
                fprintf(stderr, "TITLE not found in file %s\n", getFileName().c_str());
                return false;
            }
            while(readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                topology.title.append(line + '\n');
            }
        }

        // if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END"))
        //     topology.title = "";
        // else {
        //     topology.title = line;

        //     if(!findBlock("END")) {
        //         fprintf(stderr, "TITLE block incomplete in file %s\n", getFileName().c_str());
        //         return false;
        //     }
        // }
        // END of TITLE block
        // **************************************************


        // **************************************************
        // BEGIN of PHYSICALCONSTANTS block
        {
            if(!findBlock("PHYSICALCONSTANTS")) {
                fprintf(stderr, "PHYSICALCONSTANTS not found in file %s\n", getFileName().c_str());
                return false;
            }

            // first vacuum permittivity
            line = readLineSkippingCommentedLines('#');
            simParam.vacuumpermittivity = atof(line.c_str());
            // second hbar
            line = readLineSkippingCommentedLines('#');
            simParam.hbar = atof(line.c_str());
            // third speed of light
            line = readLineSkippingCommentedLines('#');
            simParam.speedoflight = atof(line.c_str());
            // forth boltz constant (actually it is the ideal gas constant because the units are in kJ / mol)
            line = readLineSkippingCommentedLines('#');
            simParam.boltz = atof(line.c_str());

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE PHYSICALCONSTANTS incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }

        // END of PHYSICALCONSTANTS block
        // **************************************************

        // **************************************************
        // BEGIN of ATOMTYPENAME block
        {
            if(!findBlock("ATOMTYPENAME")) {
                fprintf(stderr, "ATOMTYPENAME not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            topology.nAtomTypes = atoi(line.c_str());
            // topology.reserve(topology.nAtomTypes);


            for(i = 0; i < topology.nAtomTypes; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the ATOMTYPENAME block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    topology.atomTypeNames.push_back(line.c_str());
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE ATOMTYPENAME incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of ATOMTYPENAME block
        // **************************************************


        // **************************************************
        // BEGIN of RESNAME block
        {
            if(!findBlock("RESNAME")) {
                fprintf(stderr, "RESNAME not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            topology.nResidues = atoi(line.c_str());
            topology.resNames.reserve(topology.nResidues);

            for(i = 0; i < topology.nResidues; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the RESNAME block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    topology.resNames.push_back(line.c_str());
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE RESNAME incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of RESNAME block
        // **************************************************


        // **************************************************
        // BEGIN of SOLUTEATOM block
        {
            if(!findBlock("SOLUTEATOM")) {
                fprintf(stderr, "SOLUTEATOM not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            topology.nAtoms = atoi(line.c_str());
            topology.atomNames.reserve(topology.nAtoms);
            topology.masses.reserve(topology.nAtoms);
            topology.invMasses.reserve(topology.nAtoms);
            topology.resIds.reserve(topology.nAtoms);

            int tmp, resnm, iac, chargegroup, nexclusions, n14interactions, j;
            number mass, charge;
            std::string atomName;

            std::stringstream sstr;
            while(readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                sstr << line << '\n';
            }

            topology.excl.resize(topology.nAtoms);
            // std::vector<int> cg;
            // cg.reserve(100);

            topology.onefours.resize(topology.nAtoms);
            for(i = 0; i < topology.nAtoms; i++) {
                sstr >> tmp;
                sstr >> resnm;
                sstr >> atomName;
                sstr >> iac;
                sstr >> mass;
                sstr >> charge;
                sstr >> chargegroup;
                sstr >> nexclusions;

                if(nexclusions > 31) {
                    fprintf(stderr, "Atom %d in topology %s cannot have more than 31 exclusion\n", i+1, getFileName().c_str());
                    return false;
                }

                // fprintf(stderr, "%d %d %s %d %f %f %d %d\n", tmp, resnm, resname, iac, mass, charge, chargegroup, nexclusions);
                topology.resIds.push_back(--resnm);
                topology.atomNames.push_back(atomName);
                topology.atomTypes.push_back(--iac);

                pushBack(topology.masses, mass);
                pushBack(topology.invMasses, 1.0 / mass);
                pushBack(topology.charges, charge);
                // topology.charges.push_back(charge);

                for(j = 0; j < nexclusions; j++) {
                    sstr >> tmp;
                    tmp--;
                    addexclusion(topology.excl, i, tmp-i-1);
                    if(tmp-i-1 < 0) {
                        fprintf(stderr, "Atom %d has an excluded atom with index (%d) smaller than itself in file %s\n", i+1, tmp, getFileName().c_str());
                        return false;
                    }
                    else if(tmp-i-1 > 31){
                        fprintf(stderr, "Atom %d cannot be excluded from atom %d because (%d - %d) - 1 > 31\n", i+1, tmp, tmp, i+1);
                        return false;
                    }
                }
                // if(nexclusions > 0) {
                //     for(j = 0; j < nexclusions; j++) {
                //         fprintf(stderr, "Atom %d is excluded from %d in topology\n", i, j);
                //         fprintf(stderr, "It shows as isexcluded(topology.excl, %d, %d): %x\n", i, j, isexcluded(topology.excl, i, j));
                //     }
                //     return false;

                // }
                // topology.excl.push_back(excl);

                // TODO
                // correct this!
                // 1-4 interactions should be added to the exclusions
                // just in case!!!

                ;
                sstr >> n14interactions;
                // fprintf(stderr, "                    %d\n", n14interactions);
                for(j = 0; j < n14interactions; j++) {
                    sstr >> tmp;
                    topology.onefours[i].push_back(tmp - 1);
                }

                if (sstr.fail()) {
                    fprintf(stderr, "Bad line at SOLUTEATOM in file %s, atom %d\n", getFileName().c_str(), i+1);
                    return false;
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE SOLUTEATOM incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of SOLUTEATOM block
        // **************************************************


        // **************************************************
        // BEGIN of BONDSTRETCHTYPE block
        {
            if(!findBlock("BONDSTRETCHTYPE")) {
                fprintf(stderr, "BONDSTRETCHTYPE not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int ntypes = atoi(line.c_str());
            topology.bondStretchTypes.ntypes = ntypes;
            topology.bondStretchTypes.kh.reserve(ntypes);
            topology.bondStretchTypes.kq.reserve(ntypes);
            topology.bondStretchTypes.b0.reserve(ntypes);

            for(i = 0; i < ntypes; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the BONDSTRETCHTYPE block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    number kh, kq, b0;
                    sstr >> kq;
                    sstr >> kh;
                    sstr >> b0;
                    topology.bondStretchTypes.kq.push_back(kq);
                    topology.bondStretchTypes.kh.push_back(kh);
                    topology.bondStretchTypes.b0.push_back(b0);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE BONDSTRETCHTYPE incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of BONDSTRETCHTYPE block
        // **************************************************


        // **************************************************
        // BEGIN of BONDH block
        {
            if(!findBlock("BONDH")) {
                fprintf(stderr, "BONDH not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int nBonds = atoi(line.c_str());
            topology.bonds.nBonds = nBonds;
            topology.bonds.atoms.reserve(nBonds * 2);
            topology.bonds.type.reserve(nBonds);

            for(i = 0; i < nBonds; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the BONDH block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    int  ai, aj, type;
                    sstr >> ai;
                    sstr >> aj;
                    sstr >> type;
                    topology.bonds.atoms.push_back(--ai);
                    topology.bonds.atoms.push_back(--aj);
                    topology.bonds.type.push_back(--type);

                    // HUGE cache miss.
                    // but it should only happend while we read the file
                    // not while simulating!!!
                    topology.bonds.kq.push_back(topology.bondStretchTypes.kq[type]);
                    topology.bonds.b0.push_back(topology.bondStretchTypes.b0[type]);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE BONDH incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of BONDH block
        // **************************************************

        // **************************************************
        // BEGIN of BOND block
        {
            if(!findBlock("BOND")) {
                fprintf(stderr, "BOND not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int nBonds = atoi(line.c_str());
            topology.bonds.nBonds += nBonds;
            topology.bonds.atoms.reserve(topology.bonds.nBonds * 2);
            topology.bonds.type.reserve(topology.bonds.nBonds);

            for(i = 0; i < nBonds; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the BOND block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    int  ai, aj, type;
                    sstr >> ai;
                    sstr >> aj;
                    sstr >> type;
                    topology.bonds.atoms.push_back(--ai);
                    topology.bonds.atoms.push_back(--aj);
                    topology.bonds.type.push_back(--type);

                    // HUGE cache miss.
                    // but it should only happend while we read the file
                    // not while simulating!!!
                    topology.bonds.kq.push_back(topology.bondStretchTypes.kq[type]);
                    topology.bonds.b0.push_back(topology.bondStretchTypes.b0[type]);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE BOND incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of BOND block
        // **************************************************



        // **************************************************
        // BEGIN of BONDANGLEBENDTYPE block
        {
            if(!findBlock("BONDANGLEBENDTYPE")) {
                fprintf(stderr, "BONDANGLEBENDTYPE not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int ntypes = atoi(line.c_str());
            topology.bondAngleTypes.ntypes = ntypes;
            topology.bondAngleTypes.kh.reserve(ntypes);
            topology.bondAngleTypes.kq.reserve(ntypes);
            topology.bondAngleTypes.a0.reserve(ntypes);

            for(i = 0; i < ntypes; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the BONDANGLEBENDTYPE block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    number kh, kq, a0;
                    sstr >> kq;
                    sstr >> kh;
                    sstr >> a0;
                    topology.bondAngleTypes.kq.push_back(kq);
                    topology.bondAngleTypes.kh.push_back(kh);
                    topology.bondAngleTypes.a0.push_back(a0);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE BONDANGLEBENDTYPE incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of BONDANGLEBENDTYPE block
        // **************************************************



        // **************************************************
        // BEGIN of BONDANGLEH block
        {
            if(!findBlock("BONDANGLEH")) {
                fprintf(stderr, "BONDANGLEH not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int nBondAngles = atoi(line.c_str());
            topology.bondAngles.nBondAngles = nBondAngles;
            topology.bondAngles.atoms.reserve(nBondAngles * 3);
            topology.bondAngles.type.reserve(nBondAngles);

            for(i = 0; i < nBondAngles; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the BONDANGLEH block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    int  ai, aj, ak, type;
                    sstr >> ai;
                    sstr >> aj;
                    sstr >> ak;
                    sstr >> type;
                    topology.bondAngles.atoms.push_back(--ai);
                    topology.bondAngles.atoms.push_back(--aj);
                    topology.bondAngles.atoms.push_back(--ak);
                    topology.bondAngles.type.push_back(--type);

                    // HUGE cache miss.
                    // but it should only happend while we read the file
                    // not while simulating!!!
                    topology.bondAngles.kq.push_back(topology.bondAngleTypes.kq[type]);
                    topology.bondAngles.cos0.push_back(std::cos(topology.bondAngleTypes.a0[type]));
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE BONDANGLEH incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of BONDANGLEH block
        // **************************************************



        // **************************************************
        // BEGIN of BONDANGLEH block
        {
            if(!findBlock("BONDANGLE")) {
                fprintf(stderr, "BONDANGLE not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int nBondAngles = atoi(line.c_str());
            topology.bondAngles.nBondAngles = nBondAngles;
            topology.bondAngles.atoms.reserve(nBondAngles * 3);
            topology.bondAngles.type.reserve(nBondAngles);

            for(i = 0; i < nBondAngles; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the BONDANGLE block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    int  ai, aj, ak, type;
                    sstr >> ai;
                    sstr >> aj;
                    sstr >> ak;
                    sstr >> type;
                    topology.bondAngles.atoms.push_back(--ai);
                    topology.bondAngles.atoms.push_back(--aj);
                    topology.bondAngles.atoms.push_back(--ak);
                    topology.bondAngles.type.push_back(--type);

                    // HUGE cache miss.
                    // but it should only happend while we read the file
                    // not while simulating!!!
                    topology.bondAngles.kq.push_back(topology.bondAngleTypes.kq[type]);
                    topology.bondAngles.cos0.push_back(std::cos(topology.bondAngleTypes.a0[type]));
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE BONDANGLE incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of BONDANGLE block
        // **************************************************




        // **************************************************
        // BEGIN of IMPDIHEDRALTYPE block
        {
            if(!findBlock("IMPDIHEDRALTYPE")) {
                fprintf(stderr, "IMPDIHEDRALTYPE not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int ntypes = atoi(line.c_str());
            topology.improperDihedralTypes.ntypes = ntypes;
            topology.improperDihedralTypes.kh.reserve(ntypes);
            topology.improperDihedralTypes.a0.reserve(ntypes);

            for(i = 0; i < ntypes; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the IMPDIHEDRALTYPE block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    number kh, a0;
                    sstr >> kh;
                    sstr >> a0;
                    topology.improperDihedralTypes.kh.push_back(kh);
                    topology.improperDihedralTypes.a0.push_back(a0);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE IMPDIHEDRALTYPE incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of IMPDIHEDRALTYPE block
        // **************************************************


        // TODO
        // IMPLEMENT THE IMPROPER DIHEDRAL CODE!!!!


        // **************************************************
        // BEGIN of TORSDIHEDRALTYPE block
        {
            if(!findBlock("TORSDIHEDRALTYPE")) {
                fprintf(stderr, "TORSDIHEDRALTYPE not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int ntypes = atoi(line.c_str());
            topology.properDihedralTypes.ntypes = ntypes;
            topology.properDihedralTypes.cp.reserve(ntypes);
            topology.properDihedralTypes.pd.reserve(ntypes);
            topology.properDihedralTypes.np.reserve(ntypes);

            for(i = 0; i < ntypes; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the TORSDIHEDRALTYPE block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    number cp, pd;
                    int np;

                    sstr >> cp;
                    sstr >> pd;
                    sstr >> np;
                    topology.properDihedralTypes.cp.push_back(cp);
                    topology.properDihedralTypes.pd.push_back(pd);
                    topology.properDihedralTypes.np.push_back(np);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE TORSDIHEDRALTYPE incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of TORSDIHEDRALTYPE block
        // **************************************************




        // **************************************************
        // BEGIN of DIHEDRALH block
        {
            if(!findBlock("DIHEDRALH")) {
                fprintf(stderr, "DIHEDRALH not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int nProperDihedrals = atoi(line.c_str());
            topology.properDihedrals.nProperDihedrals = nProperDihedrals;
            topology.properDihedrals.atoms.reserve(nProperDihedrals * 4);
            topology.properDihedrals.type.reserve(nProperDihedrals);

            for(i = 0; i < nProperDihedrals; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the DIHEDRALH block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    int  ai, aj, ak, al, type;
                    sstr >> ai;
                    sstr >> aj;
                    sstr >> ak;
                    sstr >> al;
                    sstr >> type;
                    topology.properDihedrals.atoms.push_back(--ai);
                    topology.properDihedrals.atoms.push_back(--aj);
                    topology.properDihedrals.atoms.push_back(--ak);
                    topology.properDihedrals.atoms.push_back(--al);
                    topology.properDihedrals.type.push_back(--type);

                    // HUGE cache miss.
                    // but it should only happend while we read the file
                    // not while simulating!!!
                    topology.properDihedrals.cp.push_back(topology.properDihedralTypes.cp[type]);
                    topology.properDihedrals.pd.push_back(topology.properDihedralTypes.pd[type]);
                    topology.properDihedrals.np.push_back(topology.properDihedralTypes.np[type]);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE DIHEDRALH incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of DIHEDRALH block
        // **************************************************


        // **************************************************
        // BEGIN of DIHEDRAL block
        {
            if(!findBlock("DIHEDRAL")) {
                fprintf(stderr, "DIHEDRAL not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            int nProperDihedrals = atoi(line.c_str());
            topology.properDihedrals.nProperDihedrals = nProperDihedrals;
            topology.properDihedrals.atoms.reserve(nProperDihedrals * 4);
            topology.properDihedrals.type.reserve(nProperDihedrals);

            for(i = 0; i < nProperDihedrals; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the DIHEDRAL block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    int  ai, aj, ak, al, type;
                    sstr >> ai;
                    sstr >> aj;
                    sstr >> ak;
                    sstr >> al;
                    sstr >> type;
                    topology.properDihedrals.atoms.push_back(--ai);
                    topology.properDihedrals.atoms.push_back(--aj);
                    topology.properDihedrals.atoms.push_back(--ak);
                    topology.properDihedrals.atoms.push_back(--al);
                    topology.properDihedrals.type.push_back(--type);

                    // HUGE cache miss.
                    // but it should only happend while we read the file
                    // not while simulating!!!
                    topology.properDihedrals.cp.push_back(topology.properDihedralTypes.cp[type]);
                    topology.properDihedrals.pd.push_back(topology.properDihedralTypes.pd[type]);
                    topology.properDihedrals.np.push_back(topology.properDihedralTypes.np[type]);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE DIHEDRAL incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of DIHEDRAL block
        // **************************************************





        // **************************************************
        // BEGIN of LJPARAMETERS block
        {
            if(!findBlock("LJPARAMETERS")) {
                fprintf(stderr, "LJPARAMETERS not found in file %s\n", getFileName().c_str());
                return false;
            }

            line = readLineSkippingCommentedLines('#');
            topology.nLJParameters = atoi(line.c_str());
            topology.ljParameters.reserve(topology.nLJParameters);

            for(i = 0; i < topology.nLJParameters; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the LJPARAMETERS block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {
                    std::stringstream sstr(line);
                    LJParameters lj;
                    sstr >> lj.iac;
                    sstr >> lj.jac;
                    sstr >> lj.c12;
                    sstr >> lj.c6;
                    sstr >> lj.c12_14;
                    sstr >> lj.c6_14;
                    topology.ljParameters.push_back(lj);
                }
            }

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE LJPARAMETERS incomplete in file %s\n", getFileName().c_str());
                return false;
            }

            // for (std::set<int>::iterator it=topology.usedAtomTypes.begin(); it!=topology.usedAtomTypes.end(); ++it) {
            //     for (std::set<int>::iterator jt=topology.usedAtomTypes.begin(); jt!=topology.usedAtomTypes.end(); ++jt) {
            //         int lower;
            //         int max;

            //         if(*it < *jt) {
            //             lower = *it;
            //             max = *jt;
            //         }
            //         else {
            //             lower = *jt;
            //             max = *it;
            //         }

            //         int iac = ((max+1)*max)/2 + lower;
            //         // fprintf(stderr, "it, jt, iac: %d, %d, %d\n", *it, *jt, iac);
            //         topology.usedLJParameters.insert(std::pair<int,LJParameters>(iac, topology.ljParameters[iac]));
            //     }
            // }
        }
        // END of LJPARAMETERS block
        // **************************************************



        // **************************************************
        // BEGIN of SOLUTEMOLECULES block
        {
            if(!findBlock("SOLUTEMOLECULES")) {
                fprintf(stderr, "SOLUTEMOLECULES not found in file %s\n", getFileName().c_str());
                return false;
            }

            int it;
            line = readLineSkippingCommentedLines('#');
            int nMolecules = atoi(line.c_str());
            topology.nMolecules = nMolecules;
            // topology.soluteMolecules.reserve(topology.nSoluteMolecules);
            topology.molecules.molsBegin.reserve(nMolecules);
            topology.molecules.molsEnd.reserve(nMolecules);

            int nLinesToRead = nMolecules / 10;
            if(nMolecules % 10)
                nLinesToRead++;

            int molCounter = 0, lineCounter;
            for(i = 0; i < nLinesToRead; i++) {
                if(!readLineSkippingCommentedLinesAndCheckEndOfBlock(line, '#', "END")) {
                    fprintf(stderr, "Could not read the SOLUTEMOLECULES block from file %s\n", getFileName().c_str());
                    return false;
                }
                else {

                    lineCounter = 0;
                    std::stringstream sstr(line);
                    while(molCounter < nMolecules && lineCounter < 10) {
                        sstr >> it;
                        topology.molecules.molsEnd.push_back(it);
                        molCounter++;
                        lineCounter++;
                    }
                }
            }
            topology.molecules.molsBegin.push_back(0);
            std::copy(topology.molecules.molsEnd.begin(), topology.molecules.molsEnd.end() - 1, topology.molecules.molsBegin.begin() + 1);

            int nAtoms = topology.nAtoms;
            int j;

            topology.molecules.atomsToMolecules.resize(nAtoms);
            for(i = 0; i < nAtoms; i++) {

                for(j = 0; j < nMolecules; j++) {

                    const int first = topology.molecules.molsBegin[j];
                    const int last = topology.molecules.molsEnd[j];

                    if(i >= first && i < last) {
                        topology.molecules.atomsToMolecules[i] = j;
                        break;
                    }
                }
            }

            topology.molecules.comX.resize(nMolecules);
            topology.molecules.comV.resize(nMolecules);
            topology.molecules.totalKin.resize(nMolecules);
            topology.molecules.transKin.resize(nMolecules);
            topology.molecules.rotIntKin.resize(nMolecules);
            topology.molecules.kineticTensor.resize(nMolecules);
            topology.molecules.virialTensor.resize(nMolecules);

            if(!findBlock("END")) {
                fprintf(stderr, "TITLE SOLUTEMOLECULES incomplete in file %s\n", getFileName().c_str());
                return false;
            }
        }
        // END of SOLUTEMOLECULES block
        // **************************************************

        return true;
    }

};

#endif //TOPOLOGYFILEREADER_HPP