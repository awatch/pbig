#!/bin/bash


TOPDIR=`pwd`/..
EXEDIR=$TOPDIR/src/pbig
SRCDIR=$TOPDIR/src/pbig
CURDIR=`pwd`

append=

for setting in $@ 
do

echo $setting
OUTDIR=$CURDIR/$setting/$setting.set.dir
mkdir $OUTDIR
source ./$setting/$setting.set


function initSrc() {
	echo -e '\n\n\n'
	echo $algo
	echo -e '\n\n\n'

	echo "set(ALGO ${algo})" > config.cmake
	echo "set(NTPB ${NTPB})" >> config.cmake
	echo "set(COMPRESS ${compress})" >> config.cmake
	echo "set(NStreams ${nstream})" >> config.cmake
	echo "set(BALANCE ${balance})" >> config.cmake
	echo "set(MERGED ${merge})" >> config.cmake
	echo "set(NOuts ${nstream})" >> config.cmake
	echo "set(DIMS ${dim})" >> config.cmake
	echo "set(VERIFY 0)" >> config.cmake
	if [ $have_steps == "no" ]; then
		 echo "set(PROFILING 0)" >> config.cmake
	else 
		 echo "set(PROFILING 1)" >> config.cmake
	fi

}


for dim in ${DIMS[@]:0}
do

for NTPB in ${BlockArray[@]:0}
do

for nstream in ${NStreams[@]:0}
do

for compress in ${COMPRESS[@]:0}
do

for balance in ${BALANCE[@]:0}
do

for merge in ${MERGED[@]:0}
do

for algo in ${AlgosArray[@]:0} 
do

for have_steps in ${HAVE_STEPS[@]:0}
do
	if [ $have_steps != "no" ] && [ $algo != "PBIG" ]; then
		continue
	fi

	cd $SRCDIR
	initSrc;
	rm -f CMakeCache.txt
	cmake .
	make clean
	make 

	for rs in ${PossRegSizes[@]:0}
	do
	for obj in ${OBJarray[@]:0}
	do
	for device in ${DeviceArray[@]:0}
	do
	for cs in ${CSArray[@]:0}
	do
		
		sleep 1

		old_cs=$cs
		old_ntpb=$ntpb
		if [ $algo == "PRI" ] || [ $algo == "SEQ" ] || [ $algo == "CGAL_BOX" ] || [ $algo == "BULLET" ] ; then
			cs=0
		fi
		if [ $algo == "SEQ" ] || [ $algo == "CGAL_BOX" ] || [ $algo == "BULLET" ]; then
			NTPB=256
		fi


		echo // ${algo} ${dim} ${cs} ${ds} ${obj} $step $SimCase ${rs} ${NThreads} ${NTPB} ${device} ${compress} ${nstream} ${balance} ${have_steps}
		command="$EXEDIR/check.$algo -c ${cs} -g ${ds} -n ${obj} -p ${step} -t ${SimCase} -S ${rs} -r ${NThreads} -d ${device} -x ${have_steps} "
		outputfile="$OUTDIR/result_${algo}_${dim}_${cs}_${ds}_${obj}_${step}_${SimCase}_${rs}_${NThreads}_${NTPB}_${device}_${compress}_${nstream}_${balance}_${have_steps}.${append}"
		echo $command
		echo $outputfile
		$command > $outputfile 2>&1

#sleep 1

		if [ $algo == "SIMPLE" ] || [ $algo == "SEQ" ] || [ $algo == "CGAL_BOX" ] || [ $algo == "BULLET" ] ; then
			cs=$old_cs
			ntpb=$old_ntpb
			break
		fi

	done #cs
	done #device
	done #obj
	done #rs


done #have_steps
done # algo
done # merge
done # balance
done # compress
done # nstream

	if [ $algo == "SEQ" ] || [ $algo == "CGAL_BOX" ]; then
		break
	fi


done # NTPB

done # DIMS

done # settings

