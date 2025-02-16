#!/bin/bash

#read input

GEO=$(cat geo_acc.txt)

#This is a loop
for i in $GEO
	do
		SRR=$(grep $i  SraRunTable.txt | cut -d ',' -f 1|sort)
		SRR=$(echo $SRR | sed 's/ /.fastq.gz /g')
		SRR=$SRR.fastq.gz
		salmon quant -i gencode_v35_index -l A -r $SRR -o $i
	done

