#python scCODA_addPlot.py --meta XN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_allcell_auto --reference automatic --plot > sccoda_cerebellum_plot_allcell_auto.log

python scCODA_addPlot.py --meta XN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_allcell_astrocytes --reference astrocytes --plot > sccoda_cerebellum_plot_allcell_astrocytes.log &

python scCODA_addPlot.py --meta XN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_allcell_bergmann --reference bergmann  --plot > sccoda_cerebellum_plot_allcell_bergmann.log &

python scCODA_addPlot.py --meta XN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_allcell_oligodendrocytes --reference oligodendrocytes --plot > sccoda_cerebellum_plot_allcell_oligodendrocytes.log &


python scCODA_addPlot.py --meta XN_meta_remove_choroid_endothelial_fibroblast_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_removecell_astrocytes --reference automatic --plot > sccoda_cerebellum_plot_removecell_auto.log &

python scCODA_addPlot.py --meta XN_meta_remove_choroid_endothelial_fibroblast_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_removecell_astrocytes --reference astrocytes --plot > sccoda_cerebellum_plot_removecell_astrocytes.log &

python scCODA_addPlot.py --meta XN_meta_remove_choroid_endothelial_fibroblast_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_removecell_bergmann --reference bergmann  --plot > sccoda_cerebellum_plot_removecell_bergmann.log &

python scCODA_addPlot.py --meta XN_meta_remove_choroid_endothelial_fibroblast_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell.cluster --control_level homing --case_level tumbler --outdir sccoda_cerebellum_plot_removecell_oligodendrocytes --reference oligodendrocytes --plot > sccoda_cerebellum_plot_removecell_oligodendrocytes.log &
