#!/bin/bash

make distclean && make -j12 omp promd && ./promd_omp -nt 4 -imd ../senseanalysis/examples/salome/butane/butane.imd -cnf ../senseanalysis/examples/salome/butane/butane.cnf -top ../senseanalysis/examples/salome/butane/butane.top -omd butane-run.omd -trc butane-run.g96 -trf butane-run.trf

./analyses.bash butane-run.omd
