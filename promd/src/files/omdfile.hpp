/*
 * omdfile.h
 *
 *  Created on: Oct 6, 2010
 *      Author: victor
 */

#ifndef OMDFILE_H_
#define OMDFILE_H_

#include "outputfile.hpp"

#include "../configuration/conf.hpp"

#include "../simulation/simparam.hpp"
#include "../simulation/box.hpp"
// #include "../simulation/energies.hpp"
#include "../simulation/integrator/integrator.hpp"
#include "../simulation/thermostat/thermostat.hpp"
#include "../simulation/barostat/barostat.hpp"

#include "../interaction/nonbonded/nonbonded.hpp"
#include "../interaction/bonded/bonded.hpp"

#include <stdarg.h>

class OMDFile : public OutputFile
{

public:

    OMDFile();

    OMDFile(std::fstream::openmode openMode);

    virtual ~OMDFile();

    int printf (const char *format, ...)
    {
        if(dd.master()) {
            va_list arg;
            int done;

            va_start (arg, format);
            done = vfprintf(this->file, format, arg);
            va_end (arg);

            return done;
        }
        return 0;
}

    void printStep(Configuration &conf, SimParameters &simParam, Nonbonded &nonbonded, Bonded &bonded, Box &simBox, Integrator &integrator, Barostat &barostat)
    {
        if(dd.master()) {
            // fprintf(file, "# Step   Kinetic Energy / (kJ / mol)   Temperature / K   Nonbonded Energy / (kJ / mol) Bonded Energy / (kJ / mol)  Pressure / (kJ / mol / nm^3)\n");
            // fprintf(file, "%6d %16.8e %26.6f %22.8e     %22.8e %22.8e\n", integrator.step, (conf.current->kineticEnergy+conf.old->kineticEnergy) / 2.0,  (conf.current->temperature+conf.old->temperature) / 2.0, nonbonded.energy, bonded.energy, barostat.pressure());
            const number kineticEnergy = (conf.current->kineticEnergy+conf.old->kineticEnergy) * 0.5;
            const number potentialEnergy = nonbonded.energy + bonded.energy;

            fprintf(file, "# Step                           %d\n", integrator.step);
            fprintf(file, "# Timestep                       %f\n", integrator.step * simParam.dt);
            fprintf(file, "# Total Energy                  %+e\n", potentialEnergy + kineticEnergy);
            fprintf(file, "#    Kinetic Energy                %+e\n", kineticEnergy);
            fprintf(file, "#    Potential Energy              %+e\n", potentialEnergy);
            fprintf(file, "#       Nonbonded Energy              %+e\n", nonbonded.energy);
            fprintf(file, "#          Vdw Energy                    %+e\n", nonbonded.vdwEnergy);
            fprintf(file, "#          Coulomb Energy                %+e\n", nonbonded.coulEnergy);
            fprintf(file, "#          Vdw_14 Energy                 %+e\n", nonbonded.vdw14Energy);
            fprintf(file, "#          Coulomb_14 Energy             %+e\n", nonbonded.coul14Energy);
            fprintf(file, "#       Bonded Energy                    %+e\n", bonded.energy);
            fprintf(file, "#          Bond Energy                   %+e\n", bonded.bondEnergy);
            fprintf(file, "#          BondAngle Energy              %+e\n", bonded.bondAngleEnergy);
            fprintf(file, "#          Improper dihedral Energy      %+e\n", 0.0);
            fprintf(file, "#          Proper dihedral Energy        %+e\n", bonded.properDihedralEnergy);
            fprintf(file, "# Temperature                   %+e\n", (conf.current->temperature+conf.old->temperature) * 0.5);
            fprintf(file, "# Pressure                      %+e\n", barostat.pressure());

            fprintf(file, "\n");
        }
    }

    void printSimParameters(SimParameters &simParam)
    {
        if(dd.master()) {
            fprintf(file, "\n#    NSTLIM         DT        COM\n");
            if(simParam.dt >= 1e-3)
                fprintf(file, "%11d %10.3f %10d\n", simParam.nstlim, simParam.dt, simParam.com);
            else
                fprintf(file, "%11d %10e %10d\n", simParam.nstlim, simParam.dt, simParam.com);

            fprintf(file, "\n#    GENVEL         IG      TEMPI\n");
            fprintf(file, "%11d %10d %10.1f\n", simParam.genvel, simParam.ig, simParam.tempi);

            fprintf(file, "\n#     SHAKE        TOL\n");
            fprintf(file, "%11d %10.2e\n", simParam.shake, simParam.shaketol);


            fprintf(file, "\n#       NTT      TEMP0       TAUT\n");
            fprintf(file, "%11d %10.1f %10.2f\n", simParam.ntt, simParam.temp0, simParam.taut);


            fprintf(file, "\n#       NTP    KAPPA_T       TAUP   PRESSURE\n");
            fprintf(file, "%11d %10.4e %10.2f %10.6f\n", simParam.ntp, simParam.compressibility, simParam.taup, simParam.refPressureTensor(XX, XX));

            fprintf(file, "\n#  RCUT-VDW  RCUT-COUL      RLIST    NSTLIST\n");
            fprintf(file, "%11.2f %10.2f %10.2f %10d\n", simParam.rcutvdw, simParam.rcutcoul, simParam.rlist, simParam.nstlist);

            fprintf(file, "\n#    EPS-CS     EPS-RF\n");
            fprintf(file, "%11.4f %10.4f\n", simParam.epscs, simParam.epsrf);

            fprintf(file, "\n#      NTPR       NTWX       NTWF\n");
            fprintf(file, "%11d %10d %10d\n", simParam.ntpr, simParam.ntwx, simParam.ntwf);

            fprintf(file, "\n#  PROPERTY\n");
            fprintf(file, "%11d\n", simParam.property);

            fflush(file);
        }
    }
};

#endif /* OMDFILE_H_ */
