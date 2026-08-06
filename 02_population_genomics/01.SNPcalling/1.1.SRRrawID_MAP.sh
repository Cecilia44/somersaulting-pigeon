#!/bin/sh
cd /projects/sis/2.tumbler/02.SNP/1.bam/
REF=/projects/sis/2.domesticPigeon/00.ref/wangms/Pigeon_scaffolds.FINAL.fasta
for ID in xinjiang_pigeon
do
## map paired
bwa mem -t 10 -M -R  "@RG\tID:${ID}\tLB:${ID}\tPL:ILLUMINA\tSM:${ID}" ${REF} /projects/sis/2.tumbler/02.SNP/0.cleanData_link/${ID}_1.fq.gz /projects/sis/2.tumbler/02.SNP/0.cleanData_link/${ID}_2.fq.gz | samtools view -bS - -o ${ID}.bam

########## sort################
picard -Xmx25g -Djava.io.tmpdir=${ID}.pe SortSam I=${ID}.bam O=${ID}.sort.bam SORT_ORDER=coordinate
samtools index ${ID}.sort.bam

######## remove duplicate###########
picard -Xmx25g -Djava.io.tmpdir=${ID}.pe MarkDuplicates I=${ID}.sort.bam O=${ID}.sort.dedup.bam REMOVE_DUPLICATES=true METRICS_FILE=${ID}.dedup.metrics VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=50

samtools index ${ID}.sort.dedup.bam

done

