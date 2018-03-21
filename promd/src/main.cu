#include <iostream>  // For output to terminal
#include <fstream>   // For file I/O
#include <sstream>   // For file I/O
#include <iomanip>   // For output format
#include <cmath>     // For atan(), sqrt() etc.
#include <ctime>     // For timing functions
#include <string>    // For string manipulation
#include <exception> // For array manipulation

#include "main.hpp"  // For function prototypes and run-time constants.

// #include "configuration/conf.hpp"
#include "domain/domain.hpp"
#include "utilities/getopt.hpp"

#include "interaction/pairlist/pairlist.hpp"
#include "interaction/nonbonded/nonbonded.hpp"

#include "simulation/energies.hpp"

#include "simulation/integrator/integrator.hpp"
#include "simulation/thermostat/thermostat.hpp"
#include "simulation/barostat/barostat.hpp"
#include "simulation/gather.hpp"

#include "interaction/bonded/bonded.hpp"

#include "topology/molecules.hpp"

#include "files/imdfilereader.hpp"
#include "files/topologyfilereader.hpp"
#include "files/trajfiles.hpp"
#include "files/omdfile.hpp"

int main(int argc, char *argv[])
{
    startTimer(totalProgramTimer);

    dd.init(argc, argv);
    // cerr.setDomainDecomposition(dd);
    // cout.setDomainDecomposition(dd);

    GetOpt getOpt(argc, argv);

    getOpt.help = "Help intruction";
    getOpt  << ArgumentOptions("-nt",  "number of processors",                          ArgumentOptions::INT,     "max available", "help", false, true, false)
            << ArgumentOptions("-imd", "simulation parameters file",                    ArgumentOptions::STRING,     "system.imd", "help", false, true, false)
            << ArgumentOptions("-cnf", "coordinate file. Formats: gro, cnf",            ArgumentOptions::STRING,     "system.cnf", "help", false, true, false)
            << ArgumentOptions("-top", "topology file",                                 ArgumentOptions::STRING,     "system.top", "help", false, true, false)
            << ArgumentOptions("-out", "output coordinate file Formats: gro, cnf, g96", ArgumentOptions::STRING,     "system.g96", "help", false, true, false)
            << ArgumentOptions("-trc", "trajectory file Formats: gro, cnf, g96, trc",   ArgumentOptions::STRING,     "system.trc", "help", false, true, false)
            << ArgumentOptions("-trf", "trajectory file Formats: trf",                  ArgumentOptions::STRING,     "system.trf", "help", false, true, false)
            << ArgumentOptions("-omd", "log file",                                      ArgumentOptions::STRING,     "system.omd", "help", false, true, false);

    // Checking for help before anything else
    std::string help = getOpt.parse("-h");
    if(help != "false") {
        getOpt.printHelp();
        dd.finalize();
        return 0;
    }

    // Setting the number of processors per mpi rank
    // each mpi rank should take the maximum available per node
    // if the nprocs was not defined by the user
    {
        int nprocs = 1;

        std::string opt = getOpt.parse("-nt");
        if(!opt.empty())
            nprocs = atoi(opt.c_str());

#if defined(_OPENMP)
        if(nprocs < 1) {
            dd.setNumberProcs();
            nprocs = dd.nprocs;
        }
        else {
            dd.setNumberProcs(nprocs);
        }
        omp_set_num_threads(dd.nprocs);
#else
        if(nprocs != 0)
            cerr.printf("Cannot set the number of threads to %d for you are not running an OpenMP binary\n", nprocs);
        dd.nprocs = 1;
#endif
    }
    // blaze::setNumThreads(dd.nprocs);

    cerr.printf("Number of processors: %d\n"
                "          Running on: %s\n\n", dd.nprocs, dd.hostname);

    const std::string trcfilename = getOpt.parse("-trc");
    const std::string trffilename = getOpt.parse("-trf");
    const std::string omdfilename = getOpt.parse("-omd");
    const std::string cnffilename = getOpt.parse("-cnf");
    const std::string imdfilename = getOpt.parse("-imd");
    const std::string topfilename = getOpt.parse("-top");

    IMDFileReader       inputimdFile;             // simulation parameters
    OMDFile             omdfile;                  // omd file
    TOPOLOGYFileReader  topfile;                  // topo file
    TRFFileWriter       trffile;

    InputTraj  *inputCONFFile = NULL;    // input coordinate
    OutputTraj *outputTRAJFile = NULL;   // trajectory file


    SimParameters simParam;
    Configuration conf;
    Box simBox;
    PairList pairlist;
    // pairlist.nprocs = dd.nprocs;

    Nonbonded nonbonded;
    Bonded bonded;
    Integrator integrator;
    Topology topology;

    Thermostat *thermostat = NULL;
    Barostat *barostat = NULL;
    SHAKE shake;

    Energies energies;

    std::stringstream ss;

    G96FileWriter cellDebugger;

    //
    // Opening the OMD file
    //
    // This should allow us to print the evolution of the MD
    //
    if(dd.master()) {
        if(!omdfile.open(omdfilename)) {
            cerr.printf("Unable to open file: %s\n", omdfile.getFileName().c_str());
            dd.abort();
            return 1;
        }
    }

    omdfile.printf("Running %s with %d thread(s)\n\n", argv[0], dd.nprocs);
    omdfile.printf("Precision: %s\n", prec);
    omdfile.printf("Time Unit: %s\n\n", timeUnit);

    cout.printf("Running %s with %d thread%s\n\n", argv[0], dd.nprocs, dd.nprocs > 1 ? "(s)" : "");
    cout.printf("Precision: %s\n", prec);
    cout.printf("Time Unit: %s\n\n", timeUnit);

    omdfile.printf("\n PROGRAM SENSITIVITY ANALISES PERFORMS A MD-RUN\n\n");

    //
    // Reading the IMD file
    //
    // This should read the simulation conditions
    //
    if(!inputimdFile.open(imdfilename)) {
        cerr.printf("Unable to open file: %s\n", inputimdFile.getFileName().c_str());
        dd.abort();
        return 1;
    }
    if(!inputimdFile.readfile(simParam)) {
        cerr.printf("Unable to read imd file: %s\n", inputimdFile.getFileName().c_str());
        dd.abort();
        return 1;
    }
    inputimdFile.close();

    //
    // Printing the imd file
    //
    //
    omdfile.printf("=============================================\n");
    omdfile.printf(" R E A D I N G   T H E   I N P U T   D A T A\n");
    omdfile.printf("=============================================\n\n");

    omdfile.printSimParameters(simParam);

    //
    // Reading the topology file
    //
    // This should read the topology and the physical constants
    // The latter should be save at simParam
    //
    if(!topfile.open(topfilename)) {
        cerr.printf("Unable to open file: %s\n", topfilename.c_str());
        dd.abort();
        return 1;
    }
    if(!topfile.readfile(conf, simBox, topology, simParam)) {
        cerr.printf("Unable to read top file: %s\n", topfilename.c_str());
        dd.abort();
        return 1;
    }
    topfile.close();


    // cerr.printf("topology.bondAngleTypes: %d\n", topology.bondAngleTypes.ntypes);
    // int i;
    // for(i = 0; i < topology.bondAngleTypes.ntypes; i++ ) {
    //     cerr.printf("topology.bondAngleTypes[%d] %f %f %f\n", i, topology.bondAngleTypes.kq[i], topology.bondAngleTypes.kh[i], topology.bondAngleTypes.a0[i]);
    // }

    // cerr.printf("topology.bondAngles: %d\n", topology.bondAngles.nBondAngles);
    // for(i = 0; i < topology.bondAngles.nBondAngles; i++ ) {
    //     cerr.printf("topology.bondAngle[%d] has atoms %d %d %d type %d force constant %f and cos0 %f\n", i, topology.bondAngles.atoms[3*i], topology.bondAngles.atoms[3*i+1], topology.bondAngles.atoms[3*i+2], topology.bondAngles.type[i], topology.bondAngles.kq[i], topology.bondAngles.cos0[i]);
    // }

    // // return 0;


    //
    // Allocating the thermostat
    //
    // This must be done after the topology reading
    //
    if(!allocateThermostat(&thermostat, simParam, topology.nAtoms)) {
        cerr.printf("Could not identify thermostat\n");
        dd.abort();
        return 1;
    }

    //
    // Allocating the barostat
    //
    // We just read if we have thermostat or not
    //
    if(!allocateBarostat(&barostat, simParam)) {
        cerr.printf("Could not identify barostat\n");
        dd.abort();
        return 1;
    }

    //
    // Reading the coordinates and velocities file
    //
    // We should check the size of the vectors. Specially the force vectors
    //
    if(!allocateInputTRAJFile(&inputCONFFile, cnffilename))
    {
        cerr.printf("Unable to allocate file %s of type: %s\n", cnffilename.c_str(), FileHandler::getFileNameExt(cnffilename).c_str());
        dd.abort();
        return 1;
    }
    if(!inputCONFFile->open()) {
        cerr.printf("Unable to open file: %s\n", cnffilename.c_str());
        dd.abort();
        return 1;
    }

    // Reserving memory for the configuration before reading
    conf.nAtoms = topology.nAtoms;
    // conf.init();
    conf.reserve(topology.nAtoms, dd.nprocs);

    // This reading should already place the atoms inside the box
    if(!inputCONFFile->readfile(conf, simBox, topology, simParam)) {
        cerr.printf("Number of atoms is different in TOP and CNF files\n");
        dd.abort();
        return 1;
    }
    inputCONFFile->close();

    //
    // Correcting the force vector size
    //
    // In this implementation we allocate it to be: dd.nprocs * nAtoms
    //
    // conf.correctForceAllocation(dd.nprocs);
    conf.correctAllocation(dd.nprocs);


    // updating the constrained atoms list
    if(simParam.shake)
        topology.populateConstrainedAtoms();

    //
    // Finally we check the simulation parameters combined with the simulation box
    //
    //
    if(!checkSimParam(simParam, simBox)) {
        dd.abort();
        return 1;
    }
    // omdfile.printf("Done!\n");



    //
    // All the reads are done so we open the trajectory file
    //
    if(dd.master()) {
        if(!allocateOutputTRAJFile(&outputTRAJFile, trcfilename))
        {
            cerr.printf("Unable to allocate output file type: %s\n", FileHandler::getFileNameExt(trcfilename).c_str());
            dd.abort();
            return 1;
        }
        if(simParam.ntwx) {
            if(!outputTRAJFile->open()) {
                cerr.printf("Unable to open file: %s\n", trcfilename.c_str());
                dd.abort();
                return 1;
            }
            else {
                outputTRAJFile->writeTitle(conf.title);
            }
        }
        // if(!outputTRAJFile->writeFile(conf, simBox)) {
        //     cerr.printf("Unable to create file: %s\n", trcfilename.c_str());
        //     dd.abort();
        //     return 1;
        // }
    }

    if(dd.master()) {
        if(simParam.ntwf) {
            trffile.setFileName(trffilename);
            if(!trffile.open()) {
                cerr.printf("Unable to open file: %s\n", trffilename.c_str());
                dd.abort();
                return 1;
            } else {
                trffile.writeTitle(conf.title);
            }
        }

    }


    // preparing the SHAKE
    if(simParam.shake) {
        shake.prepare(topology);
    }

    // preparing the energies
    energies.extend(simParam.nstlim);


    //
    // We should generate Maxwell-Boltzmann distribution of velocities if requested
    //
    // But I haven't implemented it yet
    //
    omdfile.printf("\n===================================================================\n");
    omdfile.printf(" A T O M I C   C O O R D I N A T E S   A N D   V E L O C I T I E S\n");
    omdfile.printf("===================================================================\n\n");

    omdfile.printf("\nShould print here the coordinates and possibly generated velocities!\n");
    omdfile.printf("But I haven't implemented it yet\n\n");
    omdfile.printf("Done!\n");




    omdfile.printf("\n===================================\n");
    omdfile.printf(" P A I R L I S T   C R E A T I O N\n");
    omdfile.printf("===================================\n\n");

    //
    // Splitting box into "domains"
    //
    // Well, it was implemented like domain, now it is just a grid! :(
    //
    //
    const bool gridOK = dd.splitBoxIntoCells(conf.nAtoms, simParam.rlist, simBox);

    if(!gridOK && dd.master()) {
        cerr.printf("\nUnable to split the box into cells.\n"
            "Please decrease the number of processors or change the box size.\n\n");

        cerr.printf("current Cell decomposition: %dx%dx%d\n", dd.nCells[0], dd.nCells[1], dd.nCells[2]);
        omdfile.printf("current Cell decomposition: %dx%dx%d\n", dd.nCells[0], dd.nCells[1], dd.nCells[2]);

        omdfile.close();

        dd.abort();
        return 1;
    }

    // Printing the grid info
    omdfile.printf("     shifts: %dx%dx%d\n", dd.shifts[0], dd.shifts[1], dd.shifts[2]);
    omdfile.printf("nCellsTotal: %dx%dx%d\n", dd.nTotalCellsX(), dd.nTotalCellsY(), dd.nTotalCellsZ());
    omdfile.printf("nTotalCells: %d\n", dd.nTotalCells());

    cerr.printf("     shifts: %dx%dx%d\n", dd.shifts[0], dd.shifts[1], dd.shifts[2]);
    cerr.printf("nCellsTotal: %dx%dx%d\n", dd.nTotalCellsX(), dd.nTotalCellsY(), dd.nTotalCellsZ());
    cerr.printf("nTotalCells: %d\n", dd.nTotalCells());


    // Statistics on the pairlist update
    int stepsWithoutPairListUpdate = 0;
    int totalStepsWithoutPairListUpdate = 0;
    int numberPairListUpdate = 0;

    // Checking whether we are updating based on a fix frequency time or based on molecular motion
    const double halfRlistDiff = (simParam.rlist - std::max(simParam.rcutvdw, simParam.rcutcoul)) * 0.5;
    if(simParam.nstlist < 0) {
        pairlist.doMaxDisplacement = true;
        pairlist.maxDisplacement = 0.0;
        if(halfRlistDiff < 0.0) {
            cerr.printf("RLIST must be greater than RCUTF in order to set NSTLIST < 1\n");
            return 1;
        }
        omdfile.printf("Updating pairlist if displacements are greater than %f\n", halfRlistDiff);
        cerr.printf("Updating pairlist if displacements are greater than %f\n", halfRlistDiff);
    }
    else {
        pairlist.doMaxDisplacement = false;
        omdfile.printf("\nUpdating pairlist every: %d steps\n", simParam.nstlist);
        cerr.printf("\nUpdating pairlist every: %d steps\n", simParam.nstlist);
    }

    // Creating the pairlist
    pairlist.create(conf, simBox, topology, true);
    // pairlist.print(omdfile, 0);

    if(pairlist.doMaxDisplacement) {
        pairlist.update(conf, simBox, topology);
        numberPairListUpdate++;
    }


    // {
    //     // int i = ii;
    //     // int j = jj;
    //     // if(i > j) {
    //     //     std::swap(i,j);
    //     // }
    //     int i, j;
    //     for(i = 0; i < topology.excl.size(); i++) {
    //         for(j = i+1; j < topology.excl.size(); j++) {
    //             if((j-i < 32) && isexcluded(topology.excl, i, j-i-1)) {
    //             // if((j-i > 31) || notexcluded(topology.excl, i, j-i-1)) {
    //                 omdfile.printf("atom %d is excluded from %d\n", i, j);
    //             }
    //         }
    //     }

    // }

    // return 0;

    cerr.printf("\nPairList creation time: %4.3f %s\n\n", pairListCreateTime, timeUnit);
    omdfile.printf("\nPairList creation time: %4.3f %s\n\n", pairListCreateTime, timeUnit);

    omdfile.printf("\n===========================\n");
    omdfile.printf(" M D   S I M U L A T I O N\n");
    omdfile.printf("===========================\n\n");


    cerr.printf("Simulation time: %f ps\n\n", simParam.nstlim * simParam.dt);
    omdfile.printf("Simulation time: %f ps\n\n", simParam.nstlim * simParam.dt);

    // calc number of dof

    // compute total mass of the system
    //CALL CLCMAS(NPM,NSM,NSPM,NSP,TOTMAS,TMASS,SUBMAS,SUBMIN)

    // calc the bath
    // CALL CLCBTH(NDOF,NBATH,NBNUM,NBNDX,EKREF,DTBATH,TFACBT)

    // int i, j;
    // for(i = 0; i < topology.excl.size(); i++) {
    //     for(j = 0; j < 32; j++)
    //         if(isexcluded(topology.excl, i, j))
    //             cerr.printf("atom %d is excluded from: %d\n", i, i+j+1);
    // }


    // int i = 0;
    // int lower = topology.atomTypes[i];
    // int max = topology.atomTypes[i];

    // // topology.atomTypes[ii];

    // int iac = ((max+1)*max)/2 + lower;
    // LJParameters lj = topology.ljParameters[iac];
    // number c12 = lj.c12;
    // number c6 = lj.c6;

    // cerr.printf("atom %d, type: %d, lj pos: %d with c12: %13.6e and c6: %13.6e\n", i, lower, iac, c12, c6);
    // return 0;

    // cerr.printf("computeKineticEnergyAndTemperature\n");
    // cerr.flush();
    integrator.computeKineticEnergyAndTemperature(conf, topology, simParam, thermostat->invBoltz);
    // cerr.printf("computeKineticEnergyAndTemperature\n");
    // cerr.flush();




    // Remove initial COM motion
    // integrator.removeCOMMotion(conf, topology, simBox);

    // **************************************************
    // main MD loop
    startTimer(mdLoopTimer);
    for(integrator.step = 0; integrator.step < simParam.nstlim; integrator.step++) {


        // if(integrator.step > 0)
        //     return 0;

        // cerr.printf("integrator.step: %d\n", integrator.step);
        // cerr.flush();

        // TODO
        // review this comment
        // Place particles back into box if necessary!
        // We can do it here, or we can do it after propagating the positions
        // We chose to do it when we read the conf file, before we create the pairlist (requirement)
        // and to replace atoms inside the box after propagating the positions

        // cerr.printf("gather\n");
        // cerr.flush();
        gather(conf, simBox, topology);
        // cerr.printf("gather\n");
        // cerr.flush();

        // cerr.printf("begin of gatherAndComputePerMoleculeCOMAndCOMV\n");
        // cerr.flush();
        // cerr.printf("gatherAndComputePerMoleculeCOMAndCOMV\n");
        // cerr.flush();
        // gatherAndComputePerMoleculeCOMAndCOMV(conf, simBox, topology, simParam, thermostat->invBoltz);
        // cerr.printf("gatherAndComputePerMoleculeCOMAndCOMV\n");
        // cerr.flush();

        // cerr.printf("end of gatherAndComputePerMoleculeCOMAndCOMV\n");
        // cerr.flush();

        // compute initial kinetic energy
        // integrator.computeKineticEnergyAndTemperature(conf.vx, conf.vy, conf.vz, topology.masses, simParam.boxdof, thermostat->invBoltz, conf.nAtoms);

        // cout.printf("kineticEnergy: %f and temperature: %f\n", conf.current->kineticEnergy, conf.current->temperature);



        // cerr.printf("outputTRAJFile->writeFile\n");
        // cerr.flush();

        // Prepare virial calculation if necessary
        // IF (LDOVIR) THEN
        //     CALL PRPVIR(NATTOT,NPM,NSM,X,V,XR,TMASS,NSPM,NSP,SUBMAS,EKCM,EKCMTO,LEVERY)
        // ENDIF

        // update pair list
        // cerr.printf("pairlist\n");
        // cerr.flush();
        // outputTRAJFile->writeFile(conf, topology, simBox, integrator.step * simParam.dt, integrator.step);
        // outputTRAJFile->flush();

        // cellDebugger.writeCells(conf, topology, simBox, pairlist, integrator.step * simParam.dt, integrator.step);

        if(pairlist.doMaxDisplacement) {

            if(pairlist.maxDisplacement > halfRlistDiff) {

                // Update statistics on the pairlist update
                totalStepsWithoutPairListUpdate += stepsWithoutPairListUpdate;
                stepsWithoutPairListUpdate=0;
                numberPairListUpdate++;

                pairlist.update(conf, simBox, topology);
            }
            else
                stepsWithoutPairListUpdate++;

        } else if (integrator.step % simParam.nstlist == 0) {

            // Update statistics on the pairlist update
            totalStepsWithoutPairListUpdate += stepsWithoutPairListUpdate;
            stepsWithoutPairListUpdate=0;
            numberPairListUpdate++;

            pairlist.update(conf, simBox, topology);

        }
        else
        {
            stepsWithoutPairListUpdate++;
        }

        // cerr.printf("pairlist\n");
        // cerr.flush();

        // COMPUTE FORCES
        // computed bonded and nonbonded interactions

        // cerr.printf("nonbonded\n");
        // cerr.flush();
        nonbonded.compute(conf, topology, pairlist, simParam, simBox, barostat->compute);
        // cerr.printf("nonbonded\n");
        // cerr.flush();

        // cerr.printf("bonded\n");
        // cerr.flush();
        bonded.compute(conf, topology, simParam, simBox, barostat->compute);
        // cerr.printf("bonded\n");
        // cerr.flush();

        // write pos and box at time t to trajectory and velocities at time t - dt/2
        // cerr.printf("outputTRAJFile->writeFile\n");
        // cerr.flush();
        if (simParam.ntwx && ((integrator.step % simParam.ntwx) == 0)) {
            outputTRAJFile->writeFile(conf, topology, simBox, integrator.step * simParam.dt, integrator.step);
            outputTRAJFile->flush();
        }

        if (simParam.ntwf && ((integrator.step % simParam.ntwf) == 0)) {
            trffile.writeFile(conf, topology, simBox, integrator.step * simParam.dt, integrator.step);
            trffile.flush();
        }

        // ss << "time step: " << integrator.step;
        // conf.title = ss.str();
        // ss.str(std::string());

        // outputTRAJFile->writeFile(conf, topology, simBox);

        // outputTRAJFile->flush();
        // return 0;


    // {
    //     cerr.printf("pressute tensor\n");
    //     int i, j;
    //     for(i = 0; i < 3; i++) {
    //         for(j = 0; j < 3; j++) {
    //             cerr.printf("%f  ", barostat->pressureTensor(i, j));
    //         }
    //         cerr.printf("\n");
    //     }
    //     cerr.printf("kineticEnergy tensor\n");
    //     for(i = 0; i < 3; i++) {
    //         for(j = 0; j < 3; j++) {
    //             cerr.printf("%f  ", conf.current->kineticEnergyTensor(i, j));
    //         }
    //         cerr.printf("\n");
    //     }
    //     cerr.printf("virial tensor\n");
    //     for(i = 0; i < 3; i++) {
    //         for(j = 0; j < 3; j++) {
    //             cerr.printf("%f  ", conf.current->virialTensor(i, j));
    //         }
    //         cerr.printf("\n");
    //     }
    // }

        // Compute virial
        // IF (LDOVIR) THEN
        //    CALL CLCVIR(EKCM,VIR,PRES,EKCMTO,VIRTOT,PRESTO)
        // ENDIF
        // cerr.printf("pressureCalculation\n");
        // cerr.flush();
        barostat->pressureCalculation(conf, simBox, topology);
        // cerr.printf("pressureCalculation\n");
        // cerr.flush();


        // computing the scaling factor of Berendsen's thermostat
        // cerr.printf("thermostat->computeScale\n");
        // cerr.flush();
        thermostat->computeScale(conf.current->kineticEnergy);
        // cerr.printf("thermostat->computeScale\n");
        // cerr.flush();

        // *************************************************************************************** \\
        // START OF LEAP FROG STEP
        // propagate velocities (unconstrained) computing the center of mass velocity
        // cerr.printf("propagateVelocities\n");
        // cerr.flush();
        integrator.propagateVelocities(conf, topology, simParam);
        // cerr.printf("propagateVelocities\n");
        // cerr.flush();

        // do the following calculations
        // compute the virial, pressure, box scaling factors at the same time
        // compute the temperature scaling factor, correct the velocities at t+dt/2
        // compute the temperature and at t-dt/2, t and t+dt/2

        // scale the velocities if coupled to a bath
        // TODO integrate the computeKineticEnergyAndTemperature with the prograpate velocties
        // integrator.computeKineticEnergyAndTemperature(conf.vx, conf.vy, conf.vz, topology.masses, simParam.boxdof, thermostat->invBoltz, conf.nAtoms);
        // cerr.printf("scaleVelocities\n");
        // cerr.flush();
        thermostat->scaleVelocities(conf);
        // cerr.printf("scaleVelocities\n");
        // cerr.flush();

        // propagate coordinates (unconstrained)
        // cerr.printf("propagatePositions\n");
        // cerr.flush();
        integrator.propagatePositions(conf, pairlist, simBox, simParam.dt);
        // cerr.printf("propagatePositions\n");
        // cerr.flush();

        // END OF LEAP FROG STEP
        // *************************************************************************************** //

        // apply shake
        // cerr.printf("apply shake\n");
        // cerr.flush();
        if(simParam.shake) {
            if(!shake.apply(conf, topology, simParam, simBox)) {
                cerr.printf("Step: %d\n", integrator.step);
                cerr.flush();

                ss << "time step: " << integrator.step;
                conf.title = ss.str();
                ss.str(std::string());

                outputTRAJFile->writeFile(conf, topology, simBox, integrator.step * simParam.dt, integrator.step);
                outputTRAJFile->flush();

                omdfile.printf("SHAKE error at step: %d\n", integrator.step);
                omdfile.flush();
                dd.abort();
            }
        }
        // cerr.printf("apply shake\n");
        // cerr.flush();


        // calculate the constrained velocities

        // calculate the kinetic energies
        // C calc temperatures
        //          DO 1110 II=1,NFTMAX
        //             IF (TFACPR(II) .GE. EPS) THEN
        //                TEMP(II) = EKNOW(II)/TFACPR(II)
        //             ELSE
        //                TEMP(II) = 0.0
        //             ENDIF
        //  1110    CONTINUE

        // cerr.printf("computeKineticEnergyAndTemperature\n");
        // cerr.flush();
        integrator.computeKineticEnergyAndTemperature(conf, topology, simParam, thermostat->invBoltz);
        // cerr.printf("computeKineticEnergyAndTemperature\n");
        // cerr.flush();

        // cerr.printf("STEP %d\n", integrator.step);
        // cerr.printf("    conf.old->kineticEnergy: %f\n", conf.old->kineticEnergy);
        // cerr.printf("conf.current->kineticEnergy: %f\n", conf.current->kineticEnergy);
        // cerr.printf("      average kineticEnergy: %f\n", 0.5 * (conf.current->kineticEnergy + conf.old->kineticEnergy));

        // cerr.printf("    conf.old->temperature: %f\n", conf.old->temperature);
        // cerr.printf("conf.current->temperature: %f\n", conf.current->temperature);
        // cerr.printf("      average temperature: %f\n", 0.5 * (conf.current->temperature + conf.old->temperature));

        // omdfile.printf("STEP %d\n", integrator.step);
        // omdfile.printf("    conf.old->kineticEnergy: %f\n", conf.old->kineticEnergy);
        // omdfile.printf("conf.current->kineticEnergy: %f\n", conf.current->kineticEnergy);
        // omdfile.printf("      average kineticEnergy: %f\n", 0.5 * (conf.current->kineticEnergy + conf.old->kineticEnergy));

        // omdfile.printf("    conf.old->temperature: %f\n", conf.old->temperature);
        // omdfile.printf("conf.current->temperature: %f\n", conf.current->temperature);
        // omdfile.printf("      average temperature: %f\n", 0.5 * (conf.current->temperature + conf.old->temperature));

        //
        // integrator.putMoleculesBackIntoBox(conf, pairlist, simBox, simParam.dt);


        // Now the conf is update to the new constrained velocities
        // integrator.computeKineticEnergyAndTemperature(conf.vx, conf.vy, conf.vz, topology.masses, simParam.boxdof, thermostat->invBoltz, conf.nAtoms);
        // integrator.computeKineticEnergyAndTemperature(conf, topology, simParam, thermostat->invBoltz);
        // cerr.printf("temperature: %f, invBoltz: %f, kineticEnergy: %f\n", integrator.temperature, thermostat->invBoltz, integrator.kineticEnergy);

        // C rescale coords if we have pressure coupling
        //          IF (NTP .NE. NTPOFF) THEN
        //             CALL SCLCRD(NATTOT,X,XC,PRES,PRESTO)
        //          ENDIF
        barostat->scaleBoxAndCoordinates(conf, simBox, simParam, pairlist);

        // write energy blocks, volume, pressure and pressure scaling
        // C writing of energies and volume,pressure etc. to energy trajectory
        //          IF (NTWE .NE. 0) THEN
        //             IF (MOD(NSTEP,NTWE) .EQ. 0) THEN
        //                CALL WRTIME(IUTRJE,LFORM,NSTEP,TIME)
        //                CALL WRNRG(IUTRJE,LFORM,
        //      $              MXEWRT,ENER,
        //      $              MXCTBL,ENERES,
        //      $              NUSNRE,EPLJ,EPEL,EPRF,EPRC)
        //                CALL WRVPRT(IUTRJE,LFORM,MXVWRT,VOLPRT)
        //             ENDIF
        //          ENDIF
        // write output, if requested
        if (simParam.ntpr && ((integrator.step % simParam.ntpr) == 0)) {
            omdfile.printStep(conf, simParam, nonbonded, bonded, simBox, integrator, *barostat);
        }

        // C add energies, volprt and temperatures to averages and average square
        //          DO 210 II=1,MXETBL
        //             DTMP = ENER(II)
        //             EPSUM(II) = EPSUM(II) + DTMP
        //             EPSQ(II)  = EPSQ(II)  + DTMP**2
        //  210     CONTINUE


        // C centre of mass printing (and removal if necessary)
        //          LREMCM = (NSCM .NE. 0)
        //          IF (LREMCM) THEN
        //             LREMCM = (MOD(NSTEP+1,NSCM) .EQ. 0)
        //          ENDIF

        //          LPRLSQ = (NTPL .NE. 0)
        //          IF (LPRLSQ) THEN
        //             LPRLSQ = (MOD(NSTEP+1,NTPL) .EQ. 0)
        //          ENDIF

        //          IF (LREMCM .OR. LPRLSQ) THEN
        //             DO 80 I3 = 1,NATTO3
        //                F(I3) = X(I3) - V(I3)*DTHALF
        //  80         CONTINUE
        //          ENDIF

        //          IF (LREMCM) THEN
        //             CALL CENMAS(NATTOT,NPM,NRP,NSM,NRAM,0,NDIM,NDRMAX,F,V,
        //      $           TOTMAS,0,WMAS,WMASS,
        //      $           EKCMTO,XCM,VCM,ACM,EKROT,OCM,ICMROT)

        //             CALL STOPCM(NATTOT,NDIM,F,V,XCM,VCM,OCM,ISCROT)
        //          ENDIF


        // removing the center of mass motion
        if(simParam.com && ((integrator.step % simParam.com) == 0)) {
            // cerr.printf("removeCOMMotion at step: %d\n", integrator.step);
            // cerr.flush();
            integrator.removeCOMMotion(conf, topology, simBox);
            // cerr.printf("removeCOMMotion\n");
            // cerr.flush();
        }

        {
            const number kineticEnergy = (conf.current->kineticEnergy+conf.old->kineticEnergy) * 0.5;
            const number potentialEnergy = nonbonded.energy + bonded.energy;

            // saving the energetic data
            energies.totalEnergy[integrator.step] = potentialEnergy + kineticEnergy;
            energies.kineticEnergy[integrator.step] = kineticEnergy;
            energies.potentialEnergy[integrator.step] = potentialEnergy;
            energies.nonbondedEnergy[integrator.step] = nonbonded.energy;
            energies.nonbondedVdwEnergy[integrator.step] = nonbonded.vdwEnergy;
            energies.nonbondedCoulEnergy[integrator.step] = nonbonded.coulEnergy;
            energies.nonbondedVdw14Energy[integrator.step] = nonbonded.vdw14Energy;
            energies.nonbondedCoul14Energy[integrator.step] = nonbonded.coul14Energy;
            energies.bondedEnergy[integrator.step] = bonded.energy;
            energies.bondedBondEnergy[integrator.step] = bonded.bondEnergy;
            energies.bondedBondAngleEnergy[integrator.step] = bonded.bondAngleEnergy;
            energies.bondedImproperDihedralEnergy[integrator.step] = 0.0;
            energies.bondedProperDihedralEnergy[integrator.step] = bonded.properDihedralEnergy;
            energies.temperature[integrator.step] = (conf.current->temperature+conf.old->temperature) * 0.5;
            energies.pressure[integrator.step] = barostat->pressure();
        }

        if(integrator.step && ((integrator.step % 100) == 0)) {
            const double tmpTime = getElapsed(totalProgramTimer);
            const double etaTime = tmpTime / integrator.step * simParam.nstlim - tmpTime;
            cerr.printf("\r%6.2f%% spent: %16.3f %s. ETA til finish: %16.3f %s. ETA total time: %16.3f %s. %8.3f ns/day", 100.0 * integrator.step/simParam.nstlim, tmpTime, timeUnit, etaTime, timeUnit, etaTime+tmpTime, timeUnit, (integrator.step * simParam.dt * nsdayconstant) / tmpTime);
            // cout.printf("Avg. # steps without pairlist update: %d\n", totalStepsWithoutPairListUpdate / numberPairListUpdate);
            // cout.printf("Max displacement: %f and halfRlistDiff: %f\n", pairlist.maxDisplacement, halfRlistDiff);
        }
        // cerr.printf("end of step %d\n", integrator.step);
        // cerr.flush();
        // omdfile.flush();

    }
    cerr.printf("\n");
    // **************************************************
    addToTime(mdLoopTimer, mdLoopTime);

    energies.computeAveragesAndStdDevs();

    omdfile.printf("############################################\n");
    omdfile.printf("#              A V E R A G E S             #\n");
    omdfile.printf("############################################\n");
    omdfile.printf("# Number of Steps                %d\n", simParam.nstlim);
    omdfile.printf("# Total Energy                  %+e\n", energies.avgTotalEnergy);
    omdfile.printf("#    Kinetic Energy                %+e\n", energies.avgKineticEnergy);
    omdfile.printf("#    Potential Energy              %+e\n", energies.avgPotentialEnergy);
    omdfile.printf("#       Nonbonded Energy              %+e\n", energies.avgNonbondedEnergy);
    omdfile.printf("#          Vdw Energy                    %+e\n", energies.avgNonbondedVdwEnergy);
    omdfile.printf("#          Coulomb Energy                %+e\n", energies.avgNonbondedCoulEnergy);
    omdfile.printf("#          Vdw_14 Energy                 %+e\n", energies.avgNonbondedVdw14Energy);
    omdfile.printf("#          Coulomb_14 Energy             %+e\n", energies.avgNonbondedCoul14Energy);
    omdfile.printf("#       Bonded Energy                    %+e\n", energies.avgBondedEnergy);
    omdfile.printf("#          Bond Energy                   %+e\n", energies.avgBondedBondEnergy);
    omdfile.printf("#          BondAngle Energy              %+e\n", energies.avgBondedBondAngleEnergy);
    omdfile.printf("#          Improper dihedral Energy      %+e\n", energies.avgBondedImproperDihedralEnergy);
    omdfile.printf("#          Proper dihedral Energy        %+e\n", energies.avgBondedProperDihedralEnergy);
    omdfile.printf("# Temperature                   %+e\n", energies.avgTemperature);
    omdfile.printf("# Pressure                      %+e\n", energies.avgPressure);
    omdfile.printf("\n\n");

    omdfile.printf("############################################\n");
    omdfile.printf("#          F L U C T U A T I O N S         #\n");
    omdfile.printf("############################################\n");
    omdfile.printf("# Number of Steps                %d\n", simParam.nstlim);
    omdfile.printf("# Total Energy                  %+e\n", energies.stdDevTotalEnergy);
    omdfile.printf("#    Kinetic Energy                %+e\n", energies.stdDevKineticEnergy);
    omdfile.printf("#    Potential Energy              %+e\n", energies.stdDevPotentialEnergy);
    omdfile.printf("#       Nonbonded Energy              %+e\n", energies.stdDevNonbondedEnergy);
    omdfile.printf("#          Vdw Energy                    %+e\n", energies.stdDevNonbondedVdwEnergy);
    omdfile.printf("#          Coulomb Energy                %+e\n", energies.stdDevNonbondedCoulEnergy);
    omdfile.printf("#          Vdw_14 Energy                 %+e\n", energies.stdDevNonbondedVdw14Energy);
    omdfile.printf("#          Coulomb_14 Energy             %+e\n", energies.stdDevNonbondedCoul14Energy);
    omdfile.printf("#       Bonded Energy                    %+e\n", energies.stdDevBondedEnergy);
    omdfile.printf("#          Bond Energy                   %+e\n", energies.stdDevBondedBondEnergy);
    omdfile.printf("#          BondAngle Energy              %+e\n", energies.stdDevBondedBondAngleEnergy);
    omdfile.printf("#          Improper dihedral Energy      %+e\n", energies.stdDevBondedImproperDihedralEnergy);
    omdfile.printf("#          Proper dihedral Energy        %+e\n", energies.stdDevBondedProperDihedralEnergy);
    omdfile.printf("# Temperature                   %+e\n", energies.stdDevTemperature);
    omdfile.printf("# Pressure                      %+e\n", energies.stdDevPressure);
    omdfile.printf("\n\n");

    addToTime(totalProgramTimer, totalProgramTime);

    omdfile.printf("-------------------------------------\n");
    omdfile.printf("                         # thread(s): %d\n", dd.nprocs);
    omdfile.printf("-------------------------------------\n");
    omdfile.printf("Avg. # steps without pairlist update: %d\n", totalStepsWithoutPairListUpdate / numberPairListUpdate);
    omdfile.printf("                  # pairlist updates: %d\n", numberPairListUpdate);
    omdfile.printf("-------------------------------------\n");
    omdfile.printf("               Print trajectory time: %16.3f %s %16.3f%%\n", printTRJTime,         timeUnit, 100.0 * printTRJTime/totalProgramTime);
    omdfile.printf("              PairList creation time: %16.3f %s %16.3f%%\n", pairListCreateTime,   timeUnit, 100.0 * pairListCreateTime/totalProgramTime);
    omdfile.printf("                PairList update time: %16.3f %s %16.3f%%\n", pairListUpdateTime,   timeUnit, 100.0 * pairListUpdateTime/totalProgramTime);

    omdfile.printf("           atomic neighbor list time: %16.3f %s %16.3f%%\n", updateAtomNeighListOfLocalCellsPairlistTime,   timeUnit, 100.0 * updateAtomNeighListOfLocalCellsPairlistTime/totalProgramTime);
    omdfile.printf("            clean neighbor list time: %16.3f %s %16.3f%%\n", resetLocalCellsAndTheAtomNeighListPairlistTime,   timeUnit, 100.0 * resetLocalCellsAndTheAtomNeighListPairlistTime/totalProgramTime);
    omdfile.printf("         place atoms into cells time: %16.3f %s %16.3f%%\n", placeAtomsIntoAllCellsPairlistTime,   timeUnit, 100.0 * placeAtomsIntoAllCellsPairlistTime/totalProgramTime);

    omdfile.printf("                 Kinetic Energy time: %16.3f %s %16.3f%%\n", kineticEnergyTime,    timeUnit, 100.0 * kineticEnergyTime/totalProgramTime);
    omdfile.printf("              Nonbonded compute time: %16.3f %s %16.3f%%\n", nonbondedComputeTime, timeUnit, 100.0 * nonbondedComputeTime/totalProgramTime);
    omdfile.printf("                Nonbonded first time: %16.3f %s %16.3f%%\n", nonbondedFirstTime,   timeUnit, 100.0 * nonbondedFirstTime/totalProgramTime);
    omdfile.printf("               Nonbonded second time: %16.3f %s %16.3f%%\n", nonbondedSecondTime,  timeUnit, 100.0 * nonbondedSecondTime/totalProgramTime);
    omdfile.printf("            Nonbonded reduction time: %16.3f %s %16.3f%%\n", nonbondedReductionTime,  timeUnit, 100.0 * nonbondedReductionTime/totalProgramTime);
    omdfile.printf("                           Bond time: %16.3f %s %16.3f%%\n", bondTime,   timeUnit, 100.0 * bondTime/totalProgramTime);
    omdfile.printf("------------------------------------- %35.3f%%\n", 100.0 * (printTRJTime + pairListCreateTime + pairListUpdateTime + kineticEnergyTime + nonbondedComputeTime + bondTime)/totalProgramTime);
    omdfile.printf("                      Integrate time: %16.3f %s %16.3f%%\n", integrationTime,      timeUnit, 100.0 * integrationTime/totalProgramTime);
    omdfile.printf("                             MD time: %16.3f %s %16.3f%%\n", mdLoopTime,           timeUnit, 100.0 * mdLoopTime/totalProgramTime);
    omdfile.printf("-------------------------------------\n");
    omdfile.printf("                          Total time: %16.3f %s %16.3f%%\n", totalProgramTime, timeUnit, 100.0);
    omdfile.printf("                              ns/day: %16.3f\n", (simParam.nstlim * simParam.dt * nsdayconstant) / (totalProgramTime));
    omdfile.printf("-------------------------------------\n");

    cout.printf("-------------------------------------\n");
    cout.printf("                         # thread(s): %d\n", dd.nprocs);
    cout.printf("-------------------------------------\n");
    cout.printf("Avg. # steps without pairlist update: %d\n", totalStepsWithoutPairListUpdate / numberPairListUpdate);
    cout.printf("                  # pairlist updates: %d\n", numberPairListUpdate);
    cout.printf("-------------------------------------\n");
    cout.printf("               Print trajectory time: %16.3f %s %16.3f%%\n", printTRJTime,         timeUnit, 100.0 * printTRJTime/totalProgramTime);
    cout.printf("              PairList creation time: %16.3f %s %16.3f%%\n", pairListCreateTime,   timeUnit, 100.0 * pairListCreateTime/totalProgramTime);
    cout.printf("                PairList update time: %16.3f %s %16.3f%%\n", pairListUpdateTime,   timeUnit, 100.0 * pairListUpdateTime/totalProgramTime);

    cout.printf("           atomic neighbor list time: %16.3f %s %16.3f%%\n", updateAtomNeighListOfLocalCellsPairlistTime,   timeUnit, 100.0 * updateAtomNeighListOfLocalCellsPairlistTime/totalProgramTime);
    cout.printf("            clean neighbor list time: %16.3f %s %16.3f%%\n", resetLocalCellsAndTheAtomNeighListPairlistTime,   timeUnit, 100.0 * resetLocalCellsAndTheAtomNeighListPairlistTime/totalProgramTime);
    cout.printf("         place atoms into cells time: %16.3f %s %16.3f%%\n", placeAtomsIntoAllCellsPairlistTime,   timeUnit, 100.0 * placeAtomsIntoAllCellsPairlistTime/totalProgramTime);

    cout.printf("                 Kinetic Energy time: %16.3f %s %16.3f%%\n", kineticEnergyTime,    timeUnit, 100.0 * kineticEnergyTime/totalProgramTime);
    cout.printf("              Nonbonded compute time: %16.3f %s %16.3f%%\n", nonbondedComputeTime, timeUnit, 100.0 * nonbondedComputeTime/totalProgramTime);
    cout.printf("                Nonbonded first time: %16.3f %s %16.3f%%\n", nonbondedFirstTime,   timeUnit, 100.0 * nonbondedFirstTime/totalProgramTime);
    cout.printf("               Nonbonded second time: %16.3f %s %16.3f%%\n", nonbondedSecondTime,  timeUnit, 100.0 * nonbondedSecondTime/totalProgramTime);
    cout.printf("            Nonbonded reduction time: %16.3f %s %16.3f%%\n", nonbondedReductionTime,  timeUnit, 100.0 * nonbondedReductionTime/totalProgramTime);
    cout.printf("                           Bond time: %16.3f %s %16.3f%%\n", bondTime,   timeUnit, 100.0 * bondTime/totalProgramTime);
    cout.printf("------------------------------------- %35.3f%%\n", 100.0 * (printTRJTime + pairListCreateTime + pairListUpdateTime + kineticEnergyTime + nonbondedComputeTime + bondTime)/totalProgramTime);
    cout.printf("                      Integrate time: %16.3f %s %16.3f%%\n", integrationTime,      timeUnit, 100.0 * integrationTime/totalProgramTime);
    cout.printf("                             MD time: %16.3f %s %16.3f%%\n", mdLoopTime,           timeUnit, 100.0 * mdLoopTime/totalProgramTime);
    cout.printf("-------------------------------------\n");
    cout.printf("                          Total time: %16.3f %s %16.3f%%\n", totalProgramTime, timeUnit, 100.0);
    cout.printf("                              ns/day: %16.3f\n", (simParam.nstlim * simParam.dt * nsdayconstant) / (totalProgramTime));
    cout.printf("-------------------------------------\n");

    outputTRAJFile->close();

    omdfile.printf("Normal Termination\n");
    // omdfile << "Normal Termination\n";
    std::cout << "Normal Termination\n";

    dd.finalize();

    return 0;
}

/* VELOCIT VERLET LOOP

    // **************************************************
    // main MD loop
    startTimer(mdLoopTimer);
    for(integrator.step = 1; integrator.step <= simParam.nstlim; integrator.step++) {


        // write output, if requested
        if (simParam.ntpr && ((integrator.step % simParam.ntpr) == 0)) {
            omdfile.printStep(simParam, nonbonded, simBox, integrator);

        // pairlist.print(omdfile, integrator.step);
        }

        if (simParam.ntwx && ((integrator.step % simParam.ntwx) == 0)) {
            ss << "time step: " << integrator.step;
            conf.title = ss.str();
            ss.str(std::string());

            // printGRO(trjfile, conf, simBox);
            outputTRAJFile->writeFile(conf, simBox);
        }

        // propagate system and recompute energies
        integrator.integrate(conf, topology, pairlist, simParam, nonbonded, simBox);
        // compute the kinetic energy
        integrator.computeKineticEnergyAndTemperature(conf, topology, simParam);

        // thermostating if requested
        if(thermostat != NULL)
            thermostat->apply(integrator.kineticEnergy, conf);


        // update cell list
        if(pairlist.doMaxDisplacement) {

            if(pairlist.maxDisplacement > halfRlistDiff) {

                // Update statistics on the pairlist update
                totalStepsWithoutPairListUpdate += stepsWithoutPairListUpdate;
                stepsWithoutPairListUpdate=1;
                numberPairListUpdate++;

                // omdfile.printf("\npairlist.maxDisplacement: %f, halfRlistDiff: %f\n", pairlist.maxDisplacement, halfRlistDiff);
                // omdfile.printf("Avg. # steps without pairlist update: %d\n", totalStepsWithoutPairListUpdate / numberPairListUpdate);

                // cerr.printf("\npairlist.maxDisplacement: %f, halfRlistDiff: %f\n", pairlist.maxDisplacement, halfRlistDiff);
                // cerr.printf("Avg. # steps without pairlist update: %d\n", totalStepsWithoutPairListUpdate / numberPairListUpdate);

                pairlist.update(conf, simBox, topology);
            }
            else
                stepsWithoutPairListUpdate++;

        } else if (integrator.step % simParam.nstlist == 0) {

            // Update statistics on the pairlist update
            totalStepsWithoutPairListUpdate += stepsWithoutPairListUpdate;
            stepsWithoutPairListUpdate=1;
            numberPairListUpdate++;
            // cerr.printf("pairlist.maxDisplacement: %f, rlistDiff: %f\n", pairlist.maxDisplacement, 0.5 * rlistDiff);
            // cerr.printf("Avg. # steps without pairlist update: %d\n", totalStepsWithoutPairListUpdate / numberPairListUpdate);
            pairlist.update(conf, simBox, topology);

        }
        else
        {
            stepsWithoutPairListUpdate++;
        }

        // removing the center of mass motion
        if(simParam.com && ((integrator.step % simParam.com) == 0))
            integrator.removeCOMMotion(conf, topology);

        if((integrator.step % 100) == 0) {
            const double tmpTime = getElapsed(totalProgramTimer);
            const double etaTime = tmpTime / integrator.step * simParam.nstlim - tmpTime;
            cerr.printf("\r%6.2f%% spent: %16.3f %s. ETA til finish: %16.3f %s. ETA total time: %16.3f %s. %8.3f ns/day", 100.0 * integrator.step/simParam.nstlim, tmpTime, timeUnit, etaTime, timeUnit, etaTime+tmpTime, timeUnit, (integrator.step * simParam.dt * nsdayconstant) / tmpTime);
            // cout.printf("Avg. # steps without pairlist update: %d\n", totalStepsWithoutPairListUpdate / numberPairListUpdate);
            // cout.printf("Max displacement: %f and halfRlistDiff: %f\n", pairlist.maxDisplacement, halfRlistDiff);
        }
    }
    */
