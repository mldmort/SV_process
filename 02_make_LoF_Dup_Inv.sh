#!/bin/bash

GENCODE=/expanse/projects/sebat1/miladm/UCSD/resources/annotations/GENCODE_v42/gencode.v42.chr_patch_hapl_scaff.basic.annotation.sorted.moreInfo.bed.gz
VEP_IN=vep.tsv
BED_OUT=dup_inv_lof.bed

echo -e "key\tLOF_DUP_INV\tLOF_DUP_INV_GENES" > $BED_OUT

# vep cols:
#CHROM	POS	END	ID	SVTYPE	PLATFORM	SVLEN	SRC	GENCODE

#gencode cols:
#chr1    11868   12227   exon    lncRNA  DDX11L2 ENSG00000290825.1

NCOL1=7

col_ft=4
col_bio=5
col_sym=6
col_id=7

bedtools intersect \
 -a <(awk 'BEGIN{FS="\t";OFS="\t"} \
$5 == "DUP" || $5 == "INV" \
{ \
	key = $1"_"$2"_"$3"_"$4; \
	print $1, $2-1, $2, "bp1", key, $4, $5; \
	print $1, $3-1, $3, "bp2", key, $4, $5; \
	print $1, $2-1, $3, "bp1_bp2", key, $4, $5; \
}' $VEP_IN) \
 -b $GENCODE -wa -wb | 
awk -v ncol1=$NCOL1 \
	-v col_ft=$col_ft \
	-v col_bio=$col_bio \
	-v col_sym=$col_sym \
	-v col_id=$col_id \
'BEGIN{FS="\t";OFS="\t"} \
($(ncol1 + col_ft) == "CDS" || $(ncol1 + col_ft) == "start_codon" || $(ncol1 + col_ft) == "stop_codon" || $(ncol1 + col_ft) == "gene") && ($(ncol1 + col_bio) == "protein_coding") \
{\
	bp_mode = $4;
	key = $5;
	svtype = $7;
	ft = $(ncol1 + col_ft);
	bio = $(ncol1 + col_bio);
	sym = $(ncol1 + col_sym);
	id = $(ncol1 + col_id);

	if (ft == "CDS" || ft == "start_codon" || ft == "stop_codon") {ft_mode = "cds"} else {ft_mode = "gene"};

	keys_dict[key] = 1;\

	if (ft_mode == "cds" && bp_mode == "bp1") {
		bp1_cds_genes[svtype][key][sym] = 1;\
	}

	if (ft_mode == "cds" && bp_mode == "bp2") {
		bp2_cds_genes[svtype][key][sym] = 1;\
	}

	if (ft_mode == "cds" && bp_mode == "bp1_bp2") {
		bp1_bp2_cds_genes[svtype][key][sym] = 1;\
	}

	if (ft_mode == "gene" && bp_mode == "bp1") {
		bp1_genes[svtype][key][sym] = 1;\
	}

	if (ft_mode == "gene" && bp_mode == "bp2") {
		bp2_genes[svtype][key][sym] = 1;\
	}
}\
END{\
	for (key in keys_dict) {\
		LOF = 0;\
		mode = "";\
		sep = "";\
		gene = "";\
		# check for DUPs
		svtype = "DUP";\
		if ((key in bp1_cds_genes[svtype]) && (key in bp2_cds_genes[svtype])) {\
			for (gene in  bp1_cds_genes[svtype][key]) {\
				if (gene in bp2_cds_genes[svtype][key]) {\
					LOF = 1;\
					mode = mode sep gene;\
					sep = ",";\
					break;\
				}\
			}\
		}\

		# check for INVs
		svtype = "INV";\
		# the case where bp1 is inside a gene and bp2 is outside that gene
		if (key in bp1_genes[svtype]) {\
			fail = 0;\
			for (gene in bp1_genes[svtype][key]) {\
				if ( (key in bp2_genes[svtype]) && (gene in bp2_genes[svtype][key]) ) {\
					fail = 1;\
				}\
			}\
			if (fail == 0) {\
				LOF = 1;\
				mode = mode sep gene;\
				sep = ",";\
			}\
		}\

		# the case where bp2 is inside a gene and bp1 is outside that gene
		if (key in bp2_genes[svtype]) {\
			fail = 0;\
			for (gene in bp2_genes[svtype][key]) {\
				if ( (key in bp1_genes[svtype]) && (gene in bp1_genes[svtype][key]) ) {\
					fail = 1;\
				}\
			}\
			if (fail == 0) {\
				LOF = 1;\
				mode = mode sep gene;\
				sep = ",";\
			}\
		}\

		# the case where bp1 and bp2 are inside the same gene and bp1--bp2 spans a cds exon
		if ((key in bp1_genes[svtype]) && (key in bp2_genes[svtype]) && (key in bp1_bp2_cds_genes[svtype])) {\
			for (gene in  bp1_genes[svtype][key]) {\
				if ((gene in bp2_genes[svtype][key]) && (gene in bp1_bp2_cds_genes[svtype][key])) {\
					LOF = 1;\
					mode = mode sep gene;\
					sep = ",";\
					break;\
				}\
			}\
			
		}\

		if (LOF == 1) {print key, 1, mode}
	}
}' >> $BED_OUT
