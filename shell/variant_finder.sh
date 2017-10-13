#!/bin/bash
# version 1.1
# function:  provide different function to align and vairant calling
# author: xuzhougeng, xuzhougeng@163.com

set -e
set -u
set -o pipefail

# set environment and build index and databases
PATH=$HOME/biosoft/zgtoolkits:/usr/local/bin:/usr/bin:/bin
picard=~/biosoft/picard.jar
gatk=~/biosoft/GenomeAnalysisTK.jar
reference=/public1/wangjw/database/TAIR10/TAIR10.fa
dbsnp=/public1/wangjw/database/TAIR10/dbsnp.vcf
dict=$(echo $reference | sed s/fa/dict/)
test -f $dict || java -jar $picard CreateSequenceDictionary  R=$reference O=$dict
test -f $reference.fai || samtools faidx $reference

# process options
# set default aligner and variant_finder
aligner=bwa
vairant_finder=110

while getopts a:v:h opt
do
    case "$opt" in
	a) aligner="$OPTARG";;
	v) vairant_finder="$OPTARG";;
	h) printf "Usage: %s: [-a bwa/bowtie2] [-v 111 ]wd_path sample.txt \n" $0
	   printf "\t -a: align tool to use, bwa, bowtie2,  off means not using aligner \n"
	   printf "\t -v: gatk,bcftools or FreeBays, 1,on 0,off, default:110 means use gatk and bcftools to call vairant,\n"
	   printf "\t wd_path: your project root path \n"
	   printf "\t sample.txt: it contains the reads name and relative path. for exmaple \n"
	   printf "\t \t reads1 \t /data/seq/xxx_reads1.fq"
	   exit 0;;
	*) printf "Usage: %s: [-a value] [-h] args \n" $0
	    exit 2;;
	esac
done

shift $(($OPTIND - 1))

# align reads to reference
# $1: aligner ; $2-$#: samples
function aligner(){
    mkdir -p ${wd_dir}/analysis/align
    aligner=$1
	  shift
    for sample in $@
    do
    	result=$(basename $sample | cut -d '.' -f 1)
    	case $aligner in
    	"bwa")
    	    printf "use bwa to align reads \n"
    		test -f $reference.bwt || bwa index $reference
    		bwa mem -t 4 -M -R "@RG\tID:$(basename $sample)\tSM:$(basename $sample)\tPL:Illumina" $reference \
    		    ${wd_dir}${sample}1.fq ${wd_dir}${sample}2.fq > ${wd_dir}/analysis/align/${result}_bwa.sam;;
    	"bowtie2")
    	    printf " use bowtie2 to align reads \n"
    		bowtie2-build $reference ${reference%.*}
    		bowtie2 --rg-id $(basename $sample) --rg "PL:ILLUMINA" --rg "SM:$(basename $sample)" -x ${reference%.*} \
    		    -1 ${wd_dir}${sample}1.fq -2 ${wd_dir}${sample}2.fq -p 8 -S ${wd_dir}/analysis/align/${result}_bt.sam;;
    	*) exit 0
    	 esac
    done

}


# Vaiant Discovery
function gatk_caller(){
	for file in $@
	do
		java -Xmx16g -jar $gatk -T HaplotypeCaller \
		-R $reference  -I $file --genotyping_mode DISCOVERY \
		--standard_min_confidence_threshold_for_calling 10 \
		-A StrandAlleleCountsBySample \
		-o ../variant_gatk/${file%%.*}_raw_variants.vcf
	done

}


function bcftools_caller(){
    for file in $@
    do
        samtools mpileup  -vu -t AD,DP -f $reference $file | \
	    bcftools call -vm -Ov > ../variant_bcftools/${file%%.*}_raw_variants.vcf
    done

}


function FreeBays_caller(){
    for file in $@
	do
	    printf  "${file%%.*} \n"
        echo "freebays is under test  "
    done
}



wd_dir=$1
sample_path=$2

###########################
# use different function  #
#      to align           #
###########################

if [ $aligner != off ]
then
	samples=$(egrep -o '/.*reads' $sample_path | uniq)
	aligner $aligner $samples

	##  sort and convert

	# wd: PROJECT/analysis/align
	# output is sampleName.sorted.bam, such as A.sorted.bam B.sorted.bam
	cd ${wd_dir}/analysis/align
	samfiles=$(ls *.sam)
	echo $samfiles
	for file in $samfiles
	do
		output=${file%.*}
		samtools view -b -o ${output}.bam ${file}
		samtools sort -o ${output}.sorted.bam ${output}.bam
		samtools index ${output}.sorted.bam
	done

	## PCR-duplicate remove
	# wd: PROJECT/analysis/align_xxx
	# output is sampleName.sorted.redup.bam such as A.sorted.redup.bam B.sorted.redup.bam
	cd ${wd_dir}/analysis/align
	bamfiles=$(ls *sorted.bam)
	for file in $bamfiles
	do
		output=${file%.*}
		java -Xmx10G -jar ${picard} MarkDuplicates I=$file O=${output}.redup.bam \
			 METRICS_FILE=${output}.sort.metrics REMOVE_DUPLICATES=true \
			 ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
		samtools index ${output}.redup.bam

	done
fi

###########################
# use different function  #
#     to call variant     #
###########################


cd $wd_dir/analysis/align
pre_files=$(ls *.redup.bam)
tag1=$(echo $vairant_finder | cut -c 1)
tag2=$(echo $vairant_finder | cut -c 2)
tag3=$(echo $vairant_finder | cut -c 3)
if [ $tag1 = 0 -a $tag2 = 0 -a $tag3 = 0 ]
then
    printf "wanring: This program will not call varaint \n"
fi

if [ $tag1 = 1 ]
then
    mkdir -p ${wd_dir}/analysis/variant_gatk
	gatk_caller $pre_files
fi

if [ $tag2 = 1 ]
then
    mkdir -p ${wd_dir}/analysis/variant_bcftools
	bcftools_caller $pre_files
fi

if [ $tag3 = 1 ]
then
    FreeBays_caller $pre_files

fi
