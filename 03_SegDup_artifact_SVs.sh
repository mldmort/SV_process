#!/bin/bash

SEGDUP=/expanse/projects/sebat1/miladm/UCSD/repeats/SegDups_GRCh38_pairs.bed
VEP_IN=vep.tsv
BED_OUT=sd_artifact_svs.bed

echo -e "key\tSD_ART\tSD_SEG1\tSD_SEG2" > $BED_OUT

# vep cols:
#CHROM	POS	END	ID	SVTYPE	PLATFORM	SVLEN	SRC	GENCODE

# segdup cols:
#chrom	chromStart	chromEnd	name	strand	otherChrom	otherStart	otherEnd	fracMatch	fracMatchIndel

NCOL1=7

col_chrom=1
col_start=2
col_end=3
col_otherChrom=6
col_otherStart=7
col_otherEnd=8


bedtools intersect \
 -a <(awk 'BEGIN{FS="\t";OFS="\t"} \
$5 == "DUP" || $5 == "INV" || $5 == "DEL" \
{ \
	key_sv = $1"_"$2"_"$3"_"$4; \
	print $1, $2-1, $2, "bp1", key_sv, $4, $5; \
	print $1, $3-1, $3, "bp2", key_sv, $4, $5; \
}' $VEP_IN) \
 -b $SEGDUP -wa -wb | 
awk -v ncol1=$NCOL1 \
	-v col_chrom=$col_chrom \
	-v col_start=$col_start \
	-v col_end=$col_end \
	-v col_otherChrom=$col_otherChrom \
	-v col_otherStart=$col_otherStart \
	-v col_otherEnd=$col_otherEnd \
'BEGIN{FS="\t";OFS="\t"} \
{\
	#######
	# key_sv is recorded in a dictionary. if bp1 is intersecting with an sd, we record the key_sd1 in bp1_dict. if bp2
	# is intersecting with an sd, we record the key_sd1 again in bp2_dict. in END we go over all key_sv that are recorded
	# in keys_sv_dict, and ask is there an sd1-sd2 pair which are intersecting with bp1 and bp2?
	#######

	bp_mode = $4;
	key_sv = $5;
	svtype = $7;
	key_sd1 = $(ncol1 + col_chrom)"_"$(ncol1 + col_start)"_"$(ncol1 + col_end)
	key_sd2 = $(ncol1 + col_otherChrom)"_"$(ncol1 + col_otherStart)"_"$(ncol1 + col_otherEnd)

	keys_sv_dict[key_sv] = 1;\
	
	keys_sd_dict[key_sd1] = key_sd2;\

	if (bp_mode == "bp1") {\
		bp1_dict[key_sv][key_sd1] = 1;\
	}\

	if (bp_mode == "bp2") {\
		bp2_dict[key_sv][key_sd1] = 1;\
	}\
}\
END{\
	for (key_sv in keys_sv_dict) {\
		sd_art = 0;\
		if (key_sv in bp1_dict) {\
			for (key_sd1 in bp1_dict[key_sv]) {\
				key_sd2 = keys_sd_dict[key_sd1];\
				if ( (key_sv in bp2_dict) && (key_sd2 in bp2_dict[key_sv]) ) {\
					sd_art = 1;\
					break;\
				}\
			}\
		}\
		if (sd_art == 1) {print key_sv, 1, key_sd1, key_sd2};\
	}
}' >> $BED_OUT
