java_tmp=/projects/sis/java_tmp
gatk=/projects/software/anaconda3/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar
REF=/projects/sis/2.domesticPigeon/00.ref/wangms/Pigeon_scaffolds.FINAL.fasta

java -Xmx60g -Djava.io.tmpdir=$java_tmp -jar $gatk GenotypeGVCFs -R $REF -V gendb://../3.DBImport/database_chrID -O chrID.vcf.gz --use-new-qual-calculator --verbosity ERROR
