for CHR in `cat chr_list`

do
cp varianatFiltration.sh ${CHR}_varianatFiltration.sh

sed -i "s/chrID/${CHR}/g" ${CHR}_varianatFiltration.sh
nohup ./${CHR}_varianatFiltration.sh > ${CHR}_varianatFiltration.out &
done
