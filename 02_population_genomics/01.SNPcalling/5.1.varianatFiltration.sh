gatk=/projects/software/anaconda3/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar
REF=/projects/sis/2.domesticPigeon/00.ref/wangms/Pigeon_scaffolds.FINAL.fasta
VCF=/projects/sis/2.tumbler/02.SNP/4.genotypeGVCF

java -Xmx10g -Djava.io.tmpdir=$java_tmp -jar $gatk SelectVariants -R $REF -select-type SNP -V $VCF/chrID.vcf.gz -O raw/chrID.snp.vcf.gz
java -Xmx10g -Djava.io.tmpdir=$java_tmp -jar $gatk SelectVariants -R $REF -select-type INDEL -V $VCF/chrID.vcf.gz -O raw/chrID.indel.vcf.gz

java -Xmx10g -Djava.io.tmpdir=$java_tmp -jar $gatk VariantFiltration -V raw/chrID.snp.vcf.gz -R $REF --filter-expression "QD < 2.0" --filter-name "QD2" --filter-expression "QUAL < 40.0" --filter-name "QUAL40" --filter-expression "SOR > 3.0" --filter-name "SOR3" --filter-expression "FS > 60.0" --filter-name "FS60" --filter-expression "MQ < 40.0" --filter-name "MQ40" --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" --cluster-window-size 10  --cluster-size 3 -O filter/chrID.snp.filter.vcf.gz
java -Xmx10g -Djava.io.tmpdir=$java_tmp -jar $gatk VariantFiltration -V raw/chrID.indel.vcf.gz -R $REF --filter-expression "QD < 2.0" --filter-name "QD2" --filter-expression "QUAL < 40.0" --filter-name "QUAL40" --filter-expression "FS > 200.0" --filter-name "FS200" --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O filter/chrID.indel.filter.vcf.gz
bcftools view -f PASS filter/chrID.snp.filter.vcf.gz -o filtered/chrID.snp.filtered.vcf
bcftools view -f PASS filter/chrID.indel.filter.vcf.gz -o filtered/chrID.indel.filtered.vcf

java -Xmx10g -Djava.io.tmpdir=$java_tmp -jar $gatk MergeVcfs -I filtered/chrID.snp.filtered.vcf -I filtered/chrID.indel.filtered.vcf -O filtered/chrID.filtered.vcf

vcftools --gzvcf filtered/chrID.snp.filtered.vcf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --recode --out filtered_SNP_biAllele_missing0.5/chrID.filteredSNP_biAllele_missing0.5

vcftools --gzvcf filtered/chrID.snp.filtered.vcf --min-alleles 2 --max-alleles 2 --max-missing 1 --recode --out filtered_SNP_biAllele_missing0.5/chrID.filteredSNP_biAllele_missing1
