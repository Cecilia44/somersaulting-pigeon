for ID in `cat xinjiang_pigeon`

do
cp SRRrawID_MAP.sh ${ID}_MAP.sh

sed -i "s/xinjiang_pigeon/${ID}/g" ${ID}_MAP.sh

nohup ./${ID}_MAP.sh > ${ID}_MAP.out &
done
