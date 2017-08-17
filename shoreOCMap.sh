#!/bin/bash
set -e
set -u
set -o pipefail
# SHOREmap need to be installed firstly
# use the vcf file called from parent and F2 population as input
# these files should store under the same directory
# input: pb, with the same phenotype as the mutant; pa, without the same phenotype as the mutant
Pa=$1
Pb=$2
F2=$3
chrfile=$4

# step1: convert the vcf format
declare -i  runid=1
mkdir -p OC
for file in $Pa $Pb $F2
do
  if [ ! -f OC/convert/${runid}_converted_variant.txt ]
  then
    if [ ${file##*.} = gz ]
    then
      SHOREmap convert --marker <( zcat $file ) --folder OC/convert -runid $runid 2&> OC/s1.log
    else
      SHOREmap convert --marker  $file  --folder OC/convert -runid $runid 2&> OC/s1.log
    fi
  fi
  runid=$runid+1
done

# step2: collect markers
mkdir -p OC/marker_creation
cat OC/convert/1_converted_variant.txt OC/convert/2_converted_variant.txt OC/convert/3_converted_variant.txt \
> OC/marker_creation/combined_quality_variants.txt

# step3: extract the information of the consensus base cells of the mapping population for all the candidate markers
if [ ! -f OC/marker_creation/extracted_consensus_0.txt ]
then
  SHOREmap extract --chrsizes $chrfile --folder OC/marker_creation --marker OC/marker_creation/combined_quality_variants.txt \
  --consen OC/convert/3_converted_consen.txt -verbose 2> OC/s3_log.txt
fi

echo -e "\n step3 done"
# step4: extract quality-reference bases of one parent respective to quality-variants that have been called in the other parental genome
# output bg-ref-base-pa, which  with   the same phenotype as the mutant
if [ ! -f OC/marker_creation/extracted_quality_ref_base_21.txt ]
then
  SHOREmap extract --chrsizes $chrfile --folder OC/marker_creation --marker OC/convert/1_converted_variant.txt --extract-bg-ref \
  --consen OC/convert/1_converted_reference.txt --row-first 21  -verbose 2> OC/s4_1.log
fi
# output bg-ref-base-pb, which  without  the same phenotype as the mutant
if [ ! -f OC/marker_creation/extracted_quality_ref_base_12.txt ]
then
  SHOREmap extract --chrsizes $chrfile --folder OC/marker_creation --marker OC/convert/2_converted_variant.txt --extract-bg-ref \
  --consen OC/convert/1_converted_reference.txt --row-first 12 -verbose  2> OC/s4_2.log
fi
echo -e "\nstep4 done"


# step5 : create markers with resequencing information of the parental lines
if [ ! -f OC/marker_creation/SHOREmap_created_F2Pab_specific.txt ]
then
  SHOREmap create --chrsizes $chrfile --folder OC/marker_creation --marker OC/convert/3_converted_variant.txt \
  --marker-pa OC/convert/1_converted_variant.txt --marker-pb OC/convert/2_converted_variant.txt  \
  --bg-ref-base-pa OC/marker_creation/extracted_quality_ref_base_21.txt \
  --bg-ref-base-pb OC/marker_creation/extracted_quality_ref_base_12.txt \
  --pmarker-score 20 --pmarker-min-cov 10 --pmarker-max-cov 80 --pmarker-min-freq 0.5 \
  --bg-ref-score 20 --bg-ref-cov 10 --bg-ref-cov-max 80 --bg-ref-freq 0.5 -verbose 2> s5.log
fi

# step6:  analyze AFs of the OCF2 population
# modify Sliding-window - --window-size, --window-step,--interval-min-mean -, based on specific project
SHOREmap outcross --chrsizes $chrfile --folder OC/SHOREmap_analysis --marker  OC/marker_creation/SHOREmap_created_F2Pab_specific.txt \
--consen OC/marker_creation/extracted_consensus_0.txt --min-marker 5 -plot-boost -plot-scale --window-step 1000  \
--window-size 5000 --interval-min-mean 0.7 --interval-max-cvar 0.04 --min-coverage 20 --max-coverage 80 \
--marker-score 25 --fg-N-cov 4 -plot-win --cluster 1 -rab -background2 -verbose
