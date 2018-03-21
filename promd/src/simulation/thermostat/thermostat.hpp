#ifndef THERMOSTAT_HPP
#define THERMOSTAT_HPP

#include "../../main.hpp"
#include "../../domain/domain.hpp"

#include "../../configuration/conf.hpp"
#include "../../simulation/box.hpp"
#include "../../simulation/simparam.hpp"
#include "../../topology/topology.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif

class Thermostat
{

public:
	// TODO
	// each bath should
	number refTemperature;
	number degreesOfFreedom;
	number invBoltz;
	number invTauT;
	number timestepSize;

	number scaleFactor;
	number temperature;

public:

	virtual void computeScale(const number kineticEnergy) {};

	virtual void scaleVelocities(Configuration &conf){};

};

class AndersenThermostat : public Thermostat
{

};

class LangevinThermostat : public Thermostat
{

};

class BerendsenThermostat : public Thermostat
{

	void computeScale(const number kineticEnergy)
	{
		if(this->temperature < EPSILONLIMIT)
			this->scaleFactor = 1.0;
		else {
			number temperature = 2.0 * kineticEnergy * this->invBoltz / this->degreesOfFreedom;

			this->scaleFactor = sqrt(1.0 + ((timestepSize * invTauT) * ((this->refTemperature / temperature) - 1.0)));
		}
	}

	void scaleVelocities(Configuration &conf)
	{
		int i, nAtoms = conf.nAtoms;

		// implementation without temperature groups
        #pragma omp parallel for
        for (i = 0; i < nAtoms; i++)
            conf.current->v[i] *= scaleFactor;

        conf.current->kineticEnergy *= scaleFactor * scaleFactor;
	}

};

VHR_INLINE
bool allocateThermostat(Thermostat **thermostat, const SimParameters &simParam, const number nAtoms)
{
	switch(simParam.ntt) {

		case 1:
			*thermostat = new BerendsenThermostat();
			(*thermostat)[0].refTemperature = simParam.temp0;

			(*thermostat)[0].degreesOfFreedom = (3 * nAtoms) - simParam.boxdof;
			(*thermostat)[0].timestepSize = simParam.dt;

			(*thermostat)[0].invTauT = 1.0 / simParam.taut;
			(*thermostat)[0].invBoltz = 1.0 / simParam.boltz;

		break;

		case 0:
			*thermostat = new Thermostat();
			(*thermostat)[0].degreesOfFreedom = (3 * nAtoms) - simParam.boxdof;
			(*thermostat)[0].invBoltz = 1.0 / simParam.boltz;
		break;

		default:
			return false;
		break;
	}

    return true;
}


#endif //THERMOSTAT_HPP
