#BSUB -M 128000
#BSUB -W 1440
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -q normal
#BSUB -J K3
#BSUB -e /users/reziw3/%J.err
#BSUB -o /users/reziw3/%J.err
module load sratoolkit

SRA=$(tail -n +2 SraRunTable2.txt | cut -d ',' -f 1)

for i in $SRA
	do
		prefetch $i -X 100G
		if [ -f $i.fastq.gz ]
                then
                        echo "$i already converted"
                else
                        echo "(o) Converting SRA entry: $i"
                        #Downlaod SRA entry
                        fasterq-dump $i
                        echo "(o) Done converting: $i"
	fi
done

