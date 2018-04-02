#!/bin/bash

inputName=$1
maxNumProcs=$2
numParts=$3
dirName=$4

#Check if files exist
if [ ! -f Makefile ]; then
	echo 'Makefile does not exist.'
	echo 'Exiting job...'
	exit
fi

if [ ! -f $inputName ]; then
	echo 'Template of the input file does not exist.'
	echo 'Exiting job...'
	exit
fi

#Check that numProcs does not exceed 350
if [ $maxNumProcs -gt 350 ]; then 
	echo 'Maximum number of requested  processors exceeded, 350.'
	echo 'Exiting job...'
fi

echo $numParts
sed -i 's/WWWW/'$numParts'/g' $inputName
mv $inputName 'input.dat'
numProcs=4
source ${MODULESHOME}/init/bash
module purge 
module load openmpi
module load gcc/3.3
echo $(module list)

while [ $numProcs -le $maxNumProcs ]
do

sed -i 's/XXXX/'$numProcs'/g' sub_bsc.sh
sed -i 's/YYYY/'$dirName'_'$numProcs'/g' sub_bsc.sh
make run-main-bsc num_procs=$numProcs dir_name=$dirName'_'$numProcs

sed -i 's/'$dirName'_'$numProcs'/YYYY/g' sub_bsc.sh
sed -i 's/'$numProcs'/XXXX/g' sub_bsc.sh
numProcs=$(( 2*$numProcs ))

done


 
sed -i 's/'$numParts'/WWWW/g' input.dat
cp 'input.dat' $inputName
