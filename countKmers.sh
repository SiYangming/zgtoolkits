#!/bin/bash

set -e
set -u
set -o pipefail

#modified from NIKS script
#echo $PATH

if [ $# -lt 1 ] && [ $(wc -l $1) -lt 6 ]
then
    echo -e "You should provide a configure file including the follwing setting \n "
		echo "dir=/scratch/atif/1000_genomes/BEB		#directory for read files \n "
		echo "hawkDir=/scratch/atif/hawk			#directory where hawk is installed \n "
		echo "jellyfishDir=/home/atif/jellyfish-Hawk/bin		#directory where jellyfish is installed \n "
		echo "sortDir=/home/atif/coreutils/deps/bin		#directory where parallel sort is installed \n "
		echo "CORES=30 #number of cores to use for blast searches \n"
		echo "KMERSIZE=31 # RD:61 \n "
		exit 1
fi


config=$1
if [ ! -f ${config} ] && [ ! -r ${config} ]
then
	echo -e "$1 is not a file or is not readable \n "
fi

#  setting
dirs=$(awk 'BEGIN{FS="="}; $1~/dir/ { print $2}' $config)
hawkDir=$(awk 'BEGIN{FS="="}; $1~/hawkDir/ { print $2}' $config)
jellyfishDir=$(awk 'BEGIN{FS="="}; $1~/jellyfishDir/ { print $2}' $config)
sortDir=$(awk 'BEGIN{FS="="}; $1~/sortDir/ { print $2}' $config)
CORES=$(awk 'BEGIN{FS="="}; $1~/CORES/ { print $2}' $config)
KMERSIZE=$(awk 'BEGIN{FS="="}; $1~/KMERSIZE/ { print $2}' $config)

# be aware of newline
for dir in ${dirs}
do
    cd $dir
    OUTPREFIX=$(basename ${dir})
  	mkdir -p ${OUTPREFIX}_kmers
    echo $jellyfishDir
  	${jellyfishDir}/jellyfish count -C -o ${OUTPREFIX}_kmers/tmp -m ${KMERSIZE} -t ${CORES} -s 2G <( zcat *.fq.gz ) #change if not gzipped

  	COUNT=$(ls ${OUTPREFIX}_kmers/tmp* |wc -l)

  # merge k-mers results
  	if [ $COUNT -eq 1 ]
  	then
   		mv ${OUTPREFIX}_kmers/tmp_0 ${OUTPREFIX}_kmers_jellyfish
  	else
  		${jellyfishDir}/jellyfish merge -o ${OUTPREFIX}_kmers_jellyfish ${OUTPREFIX}_kmers/tmp*
  	fi
  	rm -rf ${OUTPREFIX}_kmers

  	COUNT=$(ls ${OUTPREFIX}_kmers_jellyfish |wc -l)

  	if [ $COUNT -eq 1 ]
  	then
  		${jellyfishDir}/jellyfish histo -f -o ${OUTPREFIX}.kmers.hist.csv -t ${CORES} ${OUTPREFIX}_kmers_jellyfish
  		awk '{print $2"\t"$1}' ${OUTPREFIX}.kmers.hist.csv > ${OUTPREFIX}_tmp
  		mv ${OUTPREFIX}_tmp ${OUTPREFIX}.kmers.hist.csv

  		awk -f ${hawkDir}/countTotalKmer.awk ${OUTPREFIX}.kmers.hist.csv >> ${dir}/total_kmers.txt

  		CUTOFF=1
  		echo $CUTOFF > ${OUTPREFIX}_cutoff.csv


  		${jellyfishDir}/jellyfish dump -c -L `expr $CUTOFF + 1` ${OUTPREFIX}_kmers_jellyfish > ${OUTPREFIX}_kmers.txt
  		${sortDir}/sort --parallel=${CORES} -n -k 1 ${OUTPREFIX}_kmers.txt > ${OUTPREFIX}_kmers_sorted.txt

  		rm ${OUTPREFIX}_kmers_jellyfish
  		rm ${OUTPREFIX}_kmers.txt

  		echo "${dir}/${OUTPREFIX}_kmers_sorted.txt" >> ${dir}/sorted_files.txt

  	fi

done
