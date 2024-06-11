#!/bin/bash

TAB_IN=test.tsv
TAB_OUT=lof_svs_constrained.tsv

col1=$(find_col2 LOF $TAB_IN)
col2=$(find_col2 X_LOF_PLI_cds     $TAB_IN) 
col3=$(find_col2 X_LOF_LOEUF_cds   $TAB_IN)
col4=$(find_col2 X_LOF_FDR_ASD_cds $TAB_IN)
col5=$(find_col2 X_LOF_FDR_DD_cds  $TAB_IN)
col6=$(find_col2 X_LOF_FDR_NDD_cds $TAB_IN)

echo "col 1: $col1"
echo "col 2: $col2"
echo "col 3: $col3"
echo "col 4: $col4"
echo "col 5: $col5"
echo "col 6: $col6"

awk -v col1=$col1 \
	-v col2=$col2 \
	-v col3=$col3 \
	-v col4=$col3 \
	-v col5=$col5 \
	-v col6=$col6 \
'BEGIN{FS="\t";OFS="\t"}$1=="CHROM" || $col2==1 || $col3==1 || $col4==1 || $col5==1 || $col6==1' $TAB_IN > $TAB_OUT
