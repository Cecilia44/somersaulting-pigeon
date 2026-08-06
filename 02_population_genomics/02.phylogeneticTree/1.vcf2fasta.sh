python vcf2phylip_fasta_nexsus.py -i pigeon_filtered_biAllele_missing1_chrAuto.vcf.gz -r -o SRR516969-70-71 --output-folder ./ --output-prefix pigeon_filtered_biAllele_noMissing_chrAuto_avoidIUPAC.recode -f

python fa-rename.py --ids list_name pigeon_filtered_biAllele_noMissing_chrAuto_avoidIUPAC.recode.min4.fasta > pigeon_filtered_biAllele_noMissing_chrAuto_avoidIUPAC.recode.min4.rename.fasta
