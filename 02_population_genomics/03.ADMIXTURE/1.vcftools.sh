vcftools --gzvcf /projects/sis/2.tumbler/02.SNP/5.varianatFiltration_combine_AliData/filtered_SNP_biAllele_missing1/pigeon388_filtered_biAllele_missing1_chrAuto.vcf.gz --keep list_239 --maf 0.01 --recode --out pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01
plink --vcf pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.vcf --allow-extra-chr --recode --out pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode
awk '{print $1 "\t" $1 ":" $4 "\t" $3 "\t" $4}' pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.map > pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.map_
rm pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.map
mv pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.map_ pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.map
plink --file pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode --allow-extra-chr --indep-pairwise 100 20 0.1
plink --file pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode --extract plink.prune.in --make-bed --allow-extra-chr --out pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.prune0.1
