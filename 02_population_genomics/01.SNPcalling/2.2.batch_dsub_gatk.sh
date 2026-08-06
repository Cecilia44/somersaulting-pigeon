for sra in `cat xinjiang_pigeon`

do
cp SRRID_gatk.sh ${sra}_gatk.sh

sed -i "s/ID/${sra}/g" ${sra}_gatk.sh

#nohup ./${sra}_gatk.sh > ${sra}_gatk.out &
done
