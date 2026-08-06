java_tmp=/projects/sis/java_tmp
gatk=/projects/software/anaconda3/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar
REF=/projects/sis/2.domesticPigeon/00.ref/wangms/Pigeon_scaffolds.FINAL.fasta
bam_dir=/projects/sis/2.tumbler/02.SNP/1.bam

java -Xmx30g -Djava.io.tmpdir=$java_tmp -jar $gatk HaplotypeCaller -R $REF -I $bam_dir/ID.sort.dedup.bam -ERC GVCF -O ID.raw.g.vcf
