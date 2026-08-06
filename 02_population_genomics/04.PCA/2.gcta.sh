gcta64 --bfile pigeon119_filtered_biAllele_missing1_maf0.01_scaffold60.recode.prune0.1 --make-grm --autosome-num 60 --thread-num 50 --out pca_matrix
gcta64 --grm pca_matrix --pca 119 --thread-num 50 --out pca
