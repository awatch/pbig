#!/bin/bash


CURDIR=`pwd`
TMP_FILE=proc.tmp

append=


for setting in $1
do

echo Process $setting

source ./$setting
OUTDIR=$CURDIR/$setting.dir


direc=${setting}
outF=$setting.csv
echo $outF
rm $outF
echo "" > $outF


ndataset=0

for device in ${DeviceArray[@]:0}
do
for cs in ${CSArray[@]:0}
do
for obj in ${OBJarray[@]:0}
do

	ndataset=`expr $ndataset + 1`
	printf "##device=%d_obj=%s_cs=%s\n" $device $obj $cs >> $outF
	echo "#$setting" | cut -d '.' -f 1 | awk '{printf "%10s,",$1}' >> $outF
	for algo in ${AlgosArray[@]:0} 
	do
		echo "$algo" | awk '{printf "%10s,",$1}' >> $outF
	done
	echo "" >> $outF

for rs in ${PossRegSizes[@]:0}
do
	echo $rs| awk '{printf "%10d,",$1}' >> $outF
for nstream in ${NStreams[@]:0}
do
for compress in ${COMPRESS[@]:0} 
do
for dim in ${DIMS[@]:0}
do
for NTPB in ${BlockArray[@]:0}
do
for balance in ${BALANCE[@]:0}
do
for merge in ${MERGED[@]:0}
do
for algo in ${AlgosArray[@]:0} 
do

	if [ $algo == "SEQ" ] || [ $algo == "CGAL_BOX" ] || [ $algo == "BULLET" ]; then
		file=$OUTDIR/result_${algo}_${dim}_no_${ds}_${obj}_${step}_${SimCase}_${rs}_${NThreads}_256_${device}_${compress}_${nstream}_${balance}_no.${append}
	elif [ $algo == "PRI" ]; then
		file=$OUTDIR/result_${algo}_${dim}_0_${ds}_${obj}_${step}_${SimCase}_${rs}_${NThreads}_256_${device}_${compress}_${nstream}_${balance}_no.${append}
	else
		file=$OUTDIR/result_${algo}_${dim}_${cs}_${ds}_${obj}_${step}_${SimCase}_${rs}_${NThreads}_${NTPB}_${device}_${compress}_${nstream}_${balance}_no.${append}
	fi

		echo $file
		del_line=`nl $file|grep 'step0'|cut -f 1`
		del_line=`echo 3+$del_line|perl -nle 'print eval $_'`
		cat $file|sed "1,${del_line}d" > filter_file

		nresults=`cat filter_file | grep -e 'step.*,' | wc -l ` 
		sum=`cat filter_file|grep -e 'step.*,'|awk 'BEGIN{FS=",";} {printf $2 "+"} END {print 0}'`
		avg=`echo "($sum)/$step" | perl -nle 'print eval $_'`


		if [ $sum == 0 ] || [ $avg == "" ] ; then
			echo $avg | awk '{printf "%10d,",$1}' >> $outF
		else
			echo $avg | awk '{printf "%10d,",$1}' >> $outF
		fi

done #algo
	echo "" >> $outF
done #merge
done #balance
done #NTPB
done #DIMS
done #compress
done #nstream
done #obj
done #rs
done #cs
done #device


echo "#ndataset=$ndataset" >> $outF
echo "#ncol=`echo ${AlgosArray[@]:0}|wc -w`" >> $outF

cp $outF $OUTDIR/$outF 
cat $outF


done
