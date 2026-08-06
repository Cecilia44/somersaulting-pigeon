python scCODA_addPlot.py --meta DN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_allcell_auto --reference automatic --plot > sccoda_cerebrum_plot_allcell_auto.log

python scCODA_addPlot.py --meta DN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_allcell_astrocytes --reference astrocytes --plot > sccoda_cerebrum_plot_allcell_astrocytes.log &

python scCODA_addPlot.py --meta DN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_allcell_microglia --reference microglia  --plot > sccoda_cerebrum_plot_allcell_microglia.log &

python scCODA_addPlot.py --meta DN_meta.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_allcell_oligodendrocytes --reference oligodendrocytes --plot > sccoda_cerebrum_plot_allcell_oligodendrocytes.log &


python scCODA_addPlot.py --meta DN_meta_remove_endothelial_vas_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_removecell_astrocytes --reference automatic --plot > sccoda_cerebrum_plot_removecell_auto.log &

python scCODA_addPlot.py --meta DN_meta_remove_endothelial_vas_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_removecell_astrocytes --reference astrocytes --plot > sccoda_cerebrum_plot_removecell_astrocytes.log &

python scCODA_addPlot.py --meta DN_meta_remove_endothelial_vas_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_removecell_microglia --reference microglia  --plot > sccoda_cerebrum_plot_removecell_microglia.log &

python scCODA_addPlot.py --meta DN_meta_remove_endothelial_vas_rbc.csv --sample_col orig.ident --group_col group --celltype_col cell_type --control_level homing --case_level tumbler --outdir sccoda_cerebrum_plot_removecell_oligodendrocytes --reference oligodendrocytes --plot > sccoda_cerebrum_plot_removecell_oligodendrocytes.log &
