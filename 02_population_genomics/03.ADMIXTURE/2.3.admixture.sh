for K in {2..20}
do
admixture -s 2345 --cv ../pigeon239_filtered_biAllele_missing1_chrAuto_maf0.01.recode.prune0.1.bed $K -j20 |tee log${K}.out
done

grep -h CV log*.out | sort -nk4  -t ' ' > cross-validation_error.txt
