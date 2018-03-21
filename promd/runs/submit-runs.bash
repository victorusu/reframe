#!/bin/bash

ntps=( 0 1 2 )
ntts=( 0 1 )

rcutfs=( 1.4 )
rlists=( 1.4 1.5 )
nstlists=( 1 5 10 -1 )

nprocs=( 0 1 8 )

cwd=$PWD
for ntp in ${ntps[@]}; do

    for ntt in ${ntts[@]}; do

        for rcutf in ${rcutfs[@]}; do

            for rlist in ${rlists[@]}; do

                for nstlist in ${nstlists[@]}; do

                    for nproc in ${nprocs[@]}; do

                        nwd=nproc_${nproc}_npt_${ntp}_ntt_${ntt}_rcutf_${rcutf}_rlist_${rlist}_nstlist_${nstlist}
                        mkdir -p ${nwd}
                        cd ${nwd}
                        cp ${cwd}/basedir/* .
                        sed -i "s/\\\$NTT/${ntt}/g;s/\\\$NTP/${ntp}/g;s/\\\$RCUTF/${rcutf}/g;s/\\\$RLIST/${rlist}/g;s/\\\$NSTLIST/${nstlist}/g;" system.imd

                        prog=${cwd}/promd_omp
                        if [ ${nproc} -eq 0 ]; then
                            prog=${cwd}/promd
                            nproc=1
                        fi

                        createjobpromd system.imd ${nproc} ${prog}
                        qsub system.job
                        #${prog} -nt ${proc}
                        ${cwd}/../analyses.bash

                        cd ${cwd}
                    done

                done

            done

        done

    done

done
