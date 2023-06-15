#!/usr/bin/env bash

#eval "$(conda shell.bash hook)"
#conda activate hclust2

hclust2.py \
	-i ../output_tables/metagenomics.metaphlan4.control.filtered.s.renamed2 \
	-o ../output_figures/metagenomics.metaphlan4.heatmap.filtered.s.2.png \
	--ftop 10 \
	--f_dist_f correlation \
	--no_fclustering \
	--s_dist_f euclidean \
	--flinkage average \
	--slinkage average \
	--flabel_size 8 \
	--slabel_size 9 \
	--no_slabels \
	--max_flabel_len 50 \
	--metadata_rows 1 \
	--legend_file ../output_figures/metagenomics.metaphlan4.heatmap.filtered.s.2.legend.png
