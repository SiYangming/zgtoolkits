#!/bin/bash
# version 1.1
# function:  provide different function to align and vairant calling
# author: xuzhougeng, xuzhougeng@163.com
# parameter: samplePath, index, reference threads

set -e
set -u
set -o pipefail

PATH=/home/wangjw/miniconda3/bin:/usr/local/bin:/usr/bin:/bin
echo "you should run this program in the project root path "
if [ $# -lt 4 ]
then
  echo -e "$# less than 4, please enter samplePath, index, reference and threads"
  exit 1
fi

# /public1/wangjw/linzhi_20170715/Data/PD/FCH2LCHCCXY_L7_wHAXPI051602-43
samplePath=$1
index=$2
reference=$3
threads=$4
#index=/public1/wangjw/linzhi_20170715/chiIndex/hirsuta.fa
#reference=/public1/wangjw/linzhi_20170715/chiIndex/hirsuta.fa

mkdir -p alignment
alignDir=alignment
# alignment
for sample in `cat $samplePath`
do
    filename=${sample##*/}
    echo "processing $filename with bwa"
    if [ ! -f  $alignDir/${filename}.sam ]
    then
        bwa mem -t 8 -B 2 $index ${sample}_1.fq.gz ${sample}_2.fq.gz >\
        $alignDir/${filename}.sam 2>  $alignDir/${filename}.log
        echo "$filename done"
    else
        echo "$filename exists"
    fi
done

# convert sort and index
# warning total memoery > threads x memeory

for sample in `cat $samplePath`
do
  echo "processing $filename with samtools"
  output=${sample##*/}
  samtools view -@ $threads -b -o $alignDir/${output}.bam $alignDir/${filename}.sam
  samtools sort -@ $threads -m 2G -o $alignDir/${output}.sorted.bam $alignDir/${output}.bam
  samtools index -@ $threads $alignDir/${output}.sorted.bam
  echo "$filename done"
done


# vairant calling with bcftools
mkdir -p variant_bcftools
varDir=variant_bcftools
for sample in `cat $samplePath`
do
  output=${sample##*/}
  echo "processing $sample with bcftools"
  samtools mpileup  -vu -t AD,DP -f $reference $alignDir/${output}.sorted.bam | \
  bcftools call -vm -Ov > $varDir/${output%%.*}_raw_variants.vcf && echo "$sample done " &
done
