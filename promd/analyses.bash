#!/bin/bash

getfilename () {

    args=`echo $1 | awk -F"." '{print NF}'`

    getfilenameout=$1

    if [ ${args} -gt 1 ]; then

        getfilenameout=`echo $1 | awk -F'.' '{for(i=1;i<NF-1;i++){printf "%s.", $i}; printf "%s\n", $(NF-1)}'`

    fi

} # end of function getfilename () {

getscriptname () {

   args=`echo $1 | awk -F"/" '{print NF}'`

   scriptname=`echo $1 | awk -F"/" -v var=$args '{print $var}'`

} # end of function getscriptname ()

printError () {

    getscriptname $0
    echo "Error in usage."
    echo "Type: $scriptname <inputfile> [nprocs] [prog]"
    echo ""
    echo "where [nproc] is the number of processors. Default 4"
    echo "where [prog] is the program. Default promd_omp"
    echo ""

    exit 1
}

omdfile="system.omd"
if [ $# -gt 1 ]; then

    printError

elif [ $# -eq 1 ]; then

    omdfile=$1

fi

tmpfile=/tmp/jkajkdhjkdaskjhfasdfhjaskdhfasdkjf$RANDON
tmptime=${tmpfile}.time
tmpprop=${tmpfile}.prop

outfile=""
if [ -r ${omdfile} ]; then

    getfilename ${omdfile}
    outfile=${getfilenameout}

else

    echo "File ${omdfile}" does not exist
    exit 1

fi


printf "#%14s\n" "time" > ${tmptime}
grep "# Timestep" ${omdfile} | awk '{print $3}' >> ${tmptime}

     property=(         totene                totkin                  totpot                  nonbonded                  bonded                     bond                     bondangle                      improperdihedral                      properdihedral     temperature     pressure)
 propertyname=("# Total Energy" "#    Kinetic Energy" "#    Potential Energy" "#       Nonbonded Energy" "#       Bonded Energy" "#          Bond Energy" "#          BondAngle Energy" "#          Improper dihedral Energy" "#          Proper dihedral Energy" "# Temperature" "# Pressure")

for((i = 0; i < ${#property[@]}; i++)); do

    printf " %15s\n" ${property[$i]} > ${tmpprop}
    col=`echo ${propertyname[$i]} | awk '{print NF}'`
    grep "${propertyname[$i]}" ${omdfile} | awk -v col=${col} '{print $(col + 1)}' >> ${tmpprop}
    paste ${tmptime} ${tmpprop} > "${outfile}-${property[$i]}.dat"

done

rm -f ${tmpfile}
rm -f ${tmptime}
rm -f ${tmpprop}
