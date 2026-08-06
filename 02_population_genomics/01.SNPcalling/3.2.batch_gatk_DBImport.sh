for CHR in `cat chr_list`

do
cp gatk_DBImport.sh ${CHR}_gatk_DBImport.sh

sed -i "s/CHR/${CHR}/g" ${CHR}_gatk_DBImport.sh

#nohup ./${CHR}_gatk_DBImport.sh > ${CHR}_gatk_DBImport.out &
done
