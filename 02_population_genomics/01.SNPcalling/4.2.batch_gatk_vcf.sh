for CHR in `cat chr_list`

do
cp chrID_genotypeGVCF.sh ${CHR}_gatk_genotypeGVCF.sh

sed -i "s/chrID/${CHR}/g" ${CHR}_gatk_genotypeGVCF.sh
nohup ./${CHR}_gatk_genotypeGVCF.sh > ${CHR}_gatk_genotypeGVCF.out &
done
