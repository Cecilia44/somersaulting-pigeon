vcftools --vcf /projects/sis/2.tumbler/07.pca_ali/212/pigeon212_filtered_biAllele_missing1_maf0.01_scaffold60.recode.vcf --keep list_119 --maf 0.01 --recode --out pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60
plink --vcf pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.vcf --allow-extra-chr --recode --out pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode
awk '{print $1 "\t" $1 ":" $4 "\t" $3 "\t" $4}' pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.map > pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.map_
rm pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.map
mv pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.map_ pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.map
plink --file pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode --allow-extra-chr --indep-pairwise 100 20 0.1
plink --file pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode --extract plink.prune.in --make-bed --allow-extra-chr --out pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.prune0.1
