java_tmp=/projects/sis/java_tmp
gatk=/projects/software/anaconda3/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar

java -Xmx75g -Djava.io.tmpdir=$java_tmp -jar $gatk \
    GenomicsDBImport \
    --genomicsdb-workspace-path database_CHR \
    --sample-name-map sample_name_map \
    --intervals CHR
