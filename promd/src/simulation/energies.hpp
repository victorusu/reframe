#ifndef ENERGIES_HPP
#define ENERGIES_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <string.h>

#include "../main.hpp"

struct Energies
{
	number avgTotalEnergy;
	number avgKineticEnergy;
	number avgPotentialEnergy;
	number avgNonbondedEnergy;
	number avgNonbondedVdwEnergy;
	number avgNonbondedCoulEnergy;
	number avgNonbondedVdw14Energy;
	number avgNonbondedCoul14Energy;
	number avgBondedEnergy;
	number avgBondedBondEnergy;
	number avgBondedBondAngleEnergy;
	number avgBondedImproperDihedralEnergy;
	number avgBondedProperDihedralEnergy;
	number avgTemperature;
	number avgPressure;

	number stdDevTotalEnergy;
	number stdDevKineticEnergy;
	number stdDevPotentialEnergy;
	number stdDevNonbondedEnergy;
	number stdDevNonbondedVdwEnergy;
	number stdDevNonbondedCoulEnergy;
	number stdDevNonbondedVdw14Energy;
	number stdDevNonbondedCoul14Energy;
	number stdDevBondedEnergy;
	number stdDevBondedBondEnergy;
	number stdDevBondedBondAngleEnergy;
	number stdDevBondedImproperDihedralEnergy;
	number stdDevBondedProperDihedralEnergy;
	number stdDevTemperature;
	number stdDevPressure;

    std::vector<number> totalEnergy;
    std::vector<number> kineticEnergy;
    std::vector<number> potentialEnergy;
    std::vector<number> nonbondedEnergy;
    std::vector<number> nonbondedVdwEnergy;
    std::vector<number> nonbondedCoulEnergy;
    std::vector<number> nonbondedVdw14Energy;
    std::vector<number> nonbondedCoul14Energy;
    std::vector<number> bondedEnergy;
    std::vector<number> bondedBondEnergy;
    std::vector<number> bondedBondAngleEnergy;
	std::vector<number> bondedImproperDihedralEnergy;
	std::vector<number> bondedProperDihedralEnergy;
    std::vector<number> temperature;
    std::vector<number> pressure;

    void extend(const int nstlim)
    {
	    totalEnergy.resize(nstlim);
	    kineticEnergy.resize(nstlim);
	    potentialEnergy.resize(nstlim);
	    nonbondedEnergy.resize(nstlim);
	    nonbondedVdwEnergy.resize(nstlim);
	    nonbondedCoulEnergy.resize(nstlim);
	    nonbondedVdw14Energy.resize(nstlim);
	    nonbondedCoul14Energy.resize(nstlim);
	    bondedEnergy.resize(nstlim);
	    bondedBondEnergy.resize(nstlim);
	    bondedBondAngleEnergy.resize(nstlim);
		bondedImproperDihedralEnergy.resize(nstlim);
		bondedProperDihedralEnergy.resize(nstlim);
	    temperature.resize(nstlim);
	    pressure.resize(nstlim);
    }

    void computeAveragesAndStdDevs()
    {
    	int i, nSteps = totalEnergy.size();
    	// int nprocs = 1;

	    avgTotalEnergy = 0.0;
	    avgKineticEnergy = 0.0;
	    avgPotentialEnergy = 0.0;
	    avgNonbondedEnergy = 0.0;
	    avgNonbondedVdwEnergy = 0.0;
	    avgNonbondedCoulEnergy = 0.0;
	    avgNonbondedVdw14Energy = 0.0;
	    avgNonbondedCoul14Energy = 0.0;
	    avgBondedEnergy = 0.0;
	    avgBondedBondEnergy = 0.0;
	    avgBondedBondAngleEnergy = 0.0;
		avgBondedImproperDihedralEnergy = 0.0;
		avgBondedProperDihedralEnergy = 0.0;
	    avgTemperature = 0.0;
	    avgPressure = 0.0;

	    stdDevTotalEnergy = 0.0;
	    stdDevKineticEnergy = 0.0;
	    stdDevPotentialEnergy = 0.0;
	    stdDevNonbondedEnergy = 0.0;
	    stdDevNonbondedVdwEnergy = 0.0;
	    stdDevNonbondedCoulEnergy = 0.0;
	    stdDevNonbondedVdw14Energy = 0.0;
	    stdDevNonbondedCoul14Energy = 0.0;
	    stdDevBondedEnergy = 0.0;
	    stdDevBondedBondEnergy = 0.0;
	    stdDevBondedBondAngleEnergy = 0.0;
		stdDevBondedImproperDihedralEnergy = 0.0;
		stdDevBondedProperDihedralEnergy = 0.0;
	    stdDevTemperature = 0.0;
	    stdDevPressure = 0.0;

    	#pragma parallel for \
					    reduction(+: avgTotalEnergy) \
					    reduction(+: avgKineticEnergy) \
					    reduction(+: avgPotentialEnergy) \
					    reduction(+: avgNonbondedEnergy) \
					    reduction(+: avgNonbondedVdwEnergy) \
					    reduction(+: avgNonbondedCoulEnergy) \
					    reduction(+: avgNonbondedVdw14Energy) \
					    reduction(+: avgNonbondedCoul14Energy) \
					    reduction(+: avgBondedEnergy) \
					    reduction(+: avgBondedBondEnergy) \
					    reduction(+: avgBondedBondAngleEnergy) \
						reduction(+: avgBondedImproperDihedralEnergy) \
						reduction(+: avgBondedProperDihedralEnergy) \
					    reduction(+: avgTemperature) \
					    reduction(+: avgPressure)
	    for(i = 0; i < nSteps; i++) {
			avgTotalEnergy += totalEnergy[i];
			avgKineticEnergy += kineticEnergy[i];
			avgPotentialEnergy += potentialEnergy[i];
			avgNonbondedEnergy += nonbondedEnergy[i];
			avgNonbondedVdwEnergy += nonbondedVdwEnergy[i];
			avgNonbondedCoulEnergy += nonbondedCoulEnergy[i];
			avgNonbondedVdw14Energy += nonbondedVdw14Energy[i];
			avgNonbondedCoul14Energy += nonbondedCoul14Energy[i];
			avgBondedEnergy += bondedEnergy[i];
			avgBondedBondEnergy += bondedBondEnergy[i];
			avgBondedBondAngleEnergy += bondedBondAngleEnergy[i];
			avgBondedImproperDihedralEnergy += bondedImproperDihedralEnergy[i];
			avgBondedProperDihedralEnergy += bondedProperDihedralEnergy[i];
			avgTemperature += temperature[i];
			avgPressure += pressure[i];
	    }

		avgTotalEnergy /= nSteps;
		avgKineticEnergy /= nSteps;
		avgPotentialEnergy /= nSteps;
		avgNonbondedEnergy /= nSteps;
		avgNonbondedVdwEnergy /= nSteps;
		avgNonbondedCoulEnergy /= nSteps;
		avgNonbondedVdw14Energy /= nSteps;
		avgNonbondedCoul14Energy /= nSteps;
		avgBondedEnergy /= nSteps;
		avgBondedBondEnergy /= nSteps;
		avgBondedBondAngleEnergy /= nSteps;
		avgBondedImproperDihedralEnergy /= nSteps;
		avgBondedProperDihedralEnergy /= nSteps;
		avgTemperature /= nSteps;
		avgPressure /= nSteps;

    	#pragma parallel for \
					    reduction(+: stdDevTotalEnergy) \
					    reduction(+: stdDevKineticEnergy) \
					    reduction(+: stdDevPotentialEnergy) \
					    reduction(+: stdDevNonbondedEnergy) \
					    reduction(+: stdDevNonbondedVdwEnergy) \
					    reduction(+: stdDevNonbondedCoulEnergy) \
					    reduction(+: stdDevNonbondedVdw14Energy) \
					    reduction(+: stdDevNonbondedCoul14Energy) \
					    reduction(+: stdDevBondedEnergy) \
					    reduction(+: stdDevBondedBondEnergy) \
					    reduction(+: stdDevBondedBondAngleEnergy) \
						reduction(+: stdDevBondedImproperDihedralEnergy) \
						reduction(+: stdDevBondedProperDihedralEnergy) \
					    reduction(+: stdDevTemperature) \
					    reduction(+: stdDevPressure)
	    for(i = 0; i < nSteps; i++) {
			stdDevTotalEnergy += (totalEnergy[i] - avgTotalEnergy) * (totalEnergy[i] - avgTotalEnergy);
			stdDevKineticEnergy += (kineticEnergy[i] - avgKineticEnergy) * (kineticEnergy[i] - avgKineticEnergy);
			stdDevPotentialEnergy += (potentialEnergy[i] - avgPotentialEnergy) * (potentialEnergy[i] - avgPotentialEnergy);
			stdDevNonbondedEnergy += (nonbondedEnergy[i] - avgNonbondedEnergy) * (nonbondedEnergy[i] - avgNonbondedEnergy);
			stdDevNonbondedVdwEnergy += (nonbondedVdwEnergy[i] - avgNonbondedVdwEnergy) * (nonbondedVdwEnergy[i] - avgNonbondedVdwEnergy);
			stdDevNonbondedCoulEnergy += (nonbondedCoulEnergy[i] - avgNonbondedCoulEnergy) * (nonbondedCoulEnergy[i] - avgNonbondedCoulEnergy);
			stdDevNonbondedVdw14Energy += (nonbondedVdw14Energy[i] - avgNonbondedVdw14Energy) * (nonbondedVdw14Energy[i] - avgNonbondedVdw14Energy);
			stdDevNonbondedCoul14Energy += (nonbondedCoul14Energy[i] - avgNonbondedCoul14Energy) * (nonbondedCoul14Energy[i] - avgNonbondedCoul14Energy);
			stdDevBondedEnergy += (bondedEnergy[i] - avgBondedEnergy) * (bondedEnergy[i] - avgBondedEnergy);
			stdDevBondedBondEnergy += (bondedBondEnergy[i] - avgBondedBondEnergy) * (bondedBondEnergy[i] - avgBondedBondEnergy);
			stdDevBondedBondAngleEnergy += (bondedBondAngleEnergy[i] - avgBondedBondAngleEnergy) * (bondedBondAngleEnergy[i] - avgBondedBondAngleEnergy);
			stdDevBondedImproperDihedralEnergy += (bondedImproperDihedralEnergy[i] - avgBondedImproperDihedralEnergy) * (bondedImproperDihedralEnergy[i] - avgBondedImproperDihedralEnergy);
			stdDevBondedProperDihedralEnergy += (bondedProperDihedralEnergy[i] - avgBondedProperDihedralEnergy) * (bondedProperDihedralEnergy[i] - avgBondedProperDihedralEnergy);
			stdDevTemperature += (temperature[i] - avgTemperature) * (temperature[i] - avgTemperature);
			stdDevPressure += (pressure[i] - avgPressure) * (pressure[i] - avgPressure);
	    }

		stdDevTotalEnergy = std::sqrt(stdDevTotalEnergy / (nSteps - 1));
		stdDevKineticEnergy = std::sqrt(stdDevKineticEnergy / (nSteps - 1));
		stdDevPotentialEnergy = std::sqrt(stdDevPotentialEnergy / (nSteps - 1));
		stdDevNonbondedEnergy = std::sqrt(stdDevNonbondedEnergy / (nSteps - 1));
		stdDevNonbondedVdwEnergy = std::sqrt(stdDevNonbondedVdwEnergy / (nSteps - 1));
		stdDevNonbondedCoulEnergy = std::sqrt(stdDevNonbondedCoulEnergy / (nSteps - 1));
		stdDevNonbondedVdw14Energy = std::sqrt(stdDevNonbondedVdw14Energy / (nSteps - 1));
		stdDevNonbondedCoul14Energy = std::sqrt(stdDevNonbondedCoul14Energy / (nSteps - 1));
		stdDevBondedEnergy = std::sqrt(stdDevBondedEnergy / (nSteps - 1));
		stdDevBondedBondEnergy = std::sqrt(stdDevBondedBondEnergy / (nSteps - 1));
		stdDevBondedBondAngleEnergy = std::sqrt(stdDevBondedBondAngleEnergy / (nSteps - 1));
		stdDevBondedImproperDihedralEnergy = std::sqrt(stdDevBondedImproperDihedralEnergy / (nSteps - 1));
		stdDevBondedProperDihedralEnergy = std::sqrt(stdDevBondedProperDihedralEnergy / (nSteps - 1));
		stdDevTemperature = std::sqrt(stdDevTemperature / (nSteps - 1));
		stdDevPressure = std::sqrt(stdDevPressure / (nSteps - 1));
    }

};


#endif //ENERGIES_HPP
