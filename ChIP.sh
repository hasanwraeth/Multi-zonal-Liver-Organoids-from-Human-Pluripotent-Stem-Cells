#!/bin/bash

#activate right env
#conda activate msa

group=(HARK1 HARK2 HARK3 HARK4 HARK5 HARK6 HARK7 HARK8)
tag=(HARK1.rep1 HARK1.rep2 HARK3.rep1 HARK3.rep2 HARK5.rep1 HARK5.rep2 HARK7.rep1 HARK7.rep2)
species=hs38
prefix=bwa

for i in {0..7}
do
## fastQC
mkdir -p fastqc/${group[$i]}
fastqc -o fastqc/${group[$i]} -t 2 \
       fastq/${tag[$i]}.fastq.gz
## trim adapter, need trim_galore be installed, here we do not do this step
mkdir -p fastq.trimmed
trim_galore -q 15 --fastqc -o fastq.trimmed/${group[$i]} fastq/${tag[$i]}.fastq.gz

## mapping by bwa
mkdir -p sam
## -t: number of threads
## -M: mark shorter split hits as secondary, this is optional for Picard compatibility.
## >: save alignment to a SAM file
## 2>: save standard error to log file
bwa mem -M -t 4 GRCh38_Ensembl_106 \
           fastq/${tag[$i]}.fastq.gz \
           > sam/$prefix.$species.${tag[$i]}.sam \
           2> bwa.$prefix.$species.${tag[$i]}.log.txt

## convert sam file to bam and clean-up
mkdir -p bam
## -q: skip alignments with MAPQ samller than 30.
samtools view -bhS -q 30 sam/$prefix.$species.${tag[$i]}.sam > bam/$prefix.$species.${tag[$i]}.bam
## sort and index the bam file for quick access.
samtools sort bam/$prefix.$species.${tag[$i]}.bam -o bam/$prefix.$species.${tag[$i]}.srt.bam
samtools index bam/$prefix.$species.${tag[$i]}.srt.bam
## remove un-sorted bam file.
rm bam/$prefix.$species.${tag[$i]}.bam

## we remove the duplicated by picard::MarkDuplicates. 
mkdir -p bam/picard
picard MarkDuplicates \
       INPUT=bam/$prefix.$species.${tag[$i]}.srt.bam \
       OUTPUT=bam/$prefix.$species.${tag[$i]}.srt.markDup.bam \
       METRICS_FILE=bam/picard/$prefix.$species.${tag[$i]}.srt.fil.picard_info.txt \
       REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
samtools index bam/$prefix.$species.${tag[$i]}.srt.markDup.bam

## use deeptools::bamCoverage to generate bigwig files
## the bw file can be viewed in IGV
mkdir -p bw
bamCoverage -b bam/$prefix.$species.${tag[$i]}.srt.markDup.bam -o bw/$prefix.$species.${tag[$i]}.bw --normalizeUsing CPM

## call peaks by macs2
mkdir -p macs2/${tag[$i]}
## -g: mappable genome size
## -q: use minimum FDR 0.05 cutoff to call significant regions.
## -B: ask MACS2 to output bedGraph files for experiment.
## --nomodel --extsize 150: the subset data is not big enough (<1000 peak) for
## macs2 to generate a model. We manually feed one.
## because we used toy genome, the genome size we set as 10M
macs2 callpeak -t bam/${prefix}.$species.${tag[$i]}.srt.markDup.bam \
               -f BAM -g 10e6 -n ${prefix}.$species.${tag[$i]} \
               --outdir macs2/${tag[$i]} -q 0.05 \
               -B --nomodel --extsize 150

done

