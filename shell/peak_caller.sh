
index=
fastq_file=

bowtie2 -p 6 -3 5 --local -x ~/reference/index/bowtie/mm10 -U ../data/SRR620205.fastq| samtools sort -O bam -o cbx7.bam
bowtie2 -p 6 -3 5 --local -x ~/reference/index/bowtie/mm10 -U ../data/SRR620206.fastq| samtools sort -O bam -o suz12.bam
bowtie2 -p 6 -3 5 --local -x ~/reference/index/bowtie/mm10 -U ../data/SRR620207.fastq| samtools sort -O bam -o RYBP.bam
bowtie2 -p 6 -3 5 --local -x ~/reference/index/bowtie/mm10 -U ../data/SRR620208.fastq| samtools sort -O bam -o IgGold.bam
bowtie2 -p 6 -3 5 --local -x ~/reference/index/bowtie/mm10 -U ../data/SRR620209.fastq| samtools sort -O bam -o IgG.bam

source activate py27
macs2 callpeak -c IgGold.bam -t suz12.bam -m 10 30 -p 1e-5 -f BAM -g mm -n suz12 2>suz12.masc2.log
macs2 callpeak -c IgGold.bam -t ring1B.bam -m 10 30 -p 1e-5 -f BAM -g mm -n ring1B 2>ring1B.masc2.log
macs2 callpeak -c IgG.bam -t cbx7.bam -m 10 30 -p 1e-5 -f BAM -g mm -n cbx7 2>cbx7.masc2.log
macs2 callpeak -c IgG.bam -t RYBP.bam -m 10 30 -p 1e-5 -f BAM -g mm -n RYBP 2>RYBP.masc2.log

ls *.bam |while read id
 do
 samtools index $id $id.bai
 done
 ls *bam |while read id
 do
 file=$(basename $id )
 sample=${file%%.*}
 echo $sample
 bamCoverage -b $id -o $sample.bw -p 10 --normalizeUsingRPKM
 computeMatrix reference-point --referencePoint TSS -b 10000 -a 10000 -R ~/annotation/CHIPseq/mm10/ucsc.refseq.bed -S $sample.bw --skipZeros -o matrix1_${sample}_TSS.gz --outFileSortedRegions regions1_${sample}_genes.bed
 plotHeatmap -m matrix1_${sample}_TSS.gz -out ${sample}.png
 done

 computeMatrix reference-point -p 10 --referencePoint TSS -b 2000 -a 2000 -S *bw -R ~/annotation/CHIPseq/mm10/ucsc.refseq.bed --skipZeros -o tmp4.mat.gz
plotHeatmap -m tmp4.mat.gz -out tmp4.merge.png
plotProfile --dpi 720 -m tmp4.mat.gz -out tmp4.profile.pdf --plotFileFormat pdf --perGroup
plotHeatmap --dpi 720 -m tmp4.mat.gz -out tmp4.merge.pdf --plotFileFormat pdf
