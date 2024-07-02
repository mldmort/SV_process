#!/bin/bash

#bcftools +split-vep <vcf> -l | less
#0       Allele
#	1       Consequence
#	2       IMPACT
#	3       SYMBOL
#	4       Gene
#5       Feature_type
#6       Feature
#	7       BIOTYPE
#	8       EXON
#	9       INTRON
#10      HGVSc
#11      HGVSp
#12      cDNA_position
#13      CDS_position
#14      Protein_position
#15      Amino_acids
#16      Codons
#17      Existing_variation
#18      ALLELE_NUM
#	19      DISTANCE
#	20      STRAND
#21      FLAGS
#22      SYMBOL_SOURCE
#23      HGNC_ID
#24      SOURCE
#	25      OverlapBP
#	26      OverlapPC
#	27      gnomadV4_AF
#	28      gnomadV4_PC
#	29      gnomadV4_POPMAX_AF
#	30      gnomadV4_name
#	31      TSSDistance
#	32      GENCODE
#	33      GNOCCHI
#	34      FB_PR_ENH_M
#	35      FB_PR_ENH_F
#	36      FANTOM_ENH

######### NOTE: you should use this instead to keep the number of total calls
#bcftools +split-vep -c Gene,BIOTYPE,SYMBOL,IMPACT,EXON,INTRON,Consequence,DISTANCE,STRAND,OverlapBP,OverlapPC,gnomadV4_AF,gnomadV4_PC,gnomadV4_POPMAX_AF,gnomadV4_name,TSSDistance,GENCODE,GNOCCHI,FB_PR_ENH_M,FB_PR_ENH_F,FANTOM_ENH LR_IL_merged_markTR50_markSD50_HW_infoAnnot_denovo_inh_vep.vcf.gz | bcftools query -f '%CHROM\n' | wc -l

DIR_IN=/expanse/projects/sebat1/miladm/UCSD/LONG_READ_COHORT/process_IL_LR
VCF_IN=$DIR_IN/LR_IL_merged_markTR50_markSD50_HW_infoAnnot_addInfo_vep.vcf.gz

date

echo "Making VEP table..."

VEP_OUT=vep.tsv
#echo -n -e "CHROM\tPOS\tEND\tID\tSVTYPE\tConsequence\tIMPACT\tSYMBOL\tGENE\tBIOTYPE\tEXON\tINTRON\tDISTANCE\tSTRAND\tOverlapBP\tOverlapPC\tgnomadV4_AF\tgnomadV4_PC\tgnomadV4_POPMAX_AF\tgnomadV4_name\tTSSDistance\tGENCODE\tGNOCCHI\tFB_PR_ENH_M\tFB_PR_ENH_F\tFANTOM_ENH\n" > $VEP_OUT

#bcftools +split-vep -X -c Gene,BIOTYPE,SYMBOL,IMPACT,EXON,INTRON,Consequence,DISTANCE,STRAND,OverlapBP,OverlapPC,gnomadV4_AF,gnomadV4_PC,gnomadV4_POPMAX_AF,gnomadV4_name,TSSDistance,GENCODE,GNOCCHI,FB_PR_ENH_M,FB_PR_ENH_F,FANTOM_ENH -i '(SVLEN<2e6 && SVLEN>-2e6) || SVTYPE=="BND"' -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%BIOTYPE\t%EXON\t%INTRON\t%DISTANCE\t%STRAND\t%OverlapBP\t%OverlapPC\t%gnomadV4_AF\t%gnomadV4_PC\t%gnomadV4_POPMAX_AF\t%gnomadV4_name\t%TSSDistance\t%GENCODE\t%GNOCCHI\t%FB_PR_ENH_M\t%FB_PR_ENH_F\t%FANTOM_ENH\n' $VCF_IN >> $VEP_OUT
################################

echo -n -e "CHROM\tPOS\tEND\tID\tSVTYPE\tPLATFORM\tSVLEN\tSRC\tGENCODE\tdenovo_LR\tdenovo_LR_LC\tMAT_INH_LR\tMAT_INH_LR_LC\tPAT_INH_LR\tPAT_INH_LR_LC\tPMAT_INH_LR_LC\tdenovo_IL\tdenovo_IL_LC\tMAT_INH_IL\tPAT_INH_IL\tPMAT_INH_IL\tAC\tPASS_STRICT\tTR50\tSD50\tHWP\tSQ5_SAMPLES\tSQ10_SAMPLES\tSQ20_SAMPLES\tSQ30_SAMPLES\tSQ40_SAMPLES\tSQ50_SAMPLES\tSQ60_SAMPLES\tSQ70_SAMPLES\tAD2_SAMPLES\tAD3_SAMPLES\tAD4_SAMPLES\tAD5_SAMPLES\tAD1_P_ALT_SAMPLES\tAD1_P_REF_SAMPLES\tAD1_M_ALT_SAMPLES\tAD1_M_REF_SAMPLES\tAD1_N_ALT_SAMPLES\tAD1_N_REF_SAMPLES\tHET_SAMPLES\tHOMALT_SAMPLES\tZERO_COV_SAMPLES\tgnomadV4_AF\tgnomadV4_PC\tgnomadV4_POPMAX_AF\tgnomadV4_name\tGNOCCHI\tFB_PR_ENH_M\tFB_PR_ENH_F\tFANTOM_ENH\tConsequence\tIMPACT\tSYMBOL\tGENE\tBIOTYPE\tDISTANCE\tTSSDistance\tAF_FOUNDER_LR\tAF_OFFSPRING_LR\tAF_FOUNDER_IL\tAF_OFFSPRING_IL\n" > $VEP_OUT
bcftools +split-vep -X -c Gene,BIOTYPE,SYMBOL,IMPACT,EXON,INTRON,Consequence,DISTANCE,STRAND,OverlapBP,OverlapPC,gnomadV4_AF,gnomadV4_PC,gnomadV4_POPMAX_AF,gnomadV4_name,TSSDistance,GENCODE,GNOCCHI,FB_PR_ENH_M,FB_PR_ENH_F,FANTOM_ENH \
 -i '(SVLEN<2e6 && SVLEN>-2e6) || SVTYPE=="BND" || SVTYPE=="."' \
 -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%PLATFORM\t%SVLEN\t%SRC\t%GENCODE\t%denovo_LR\t%denovo_LR_LC\t%MAT_INH_LR\t%MAT_INH_LR_LC\t%PAT_INH_LR\t%PAT_INH_LR_LC\t%PMAT_INH_LR_LC\t%denovo_IL\t%denovo_IL_LC\t%MAT_INH_IL\t%PAT_INH_IL\t%PMAT_INH_IL\t%AC\t%PASS_STRICT\t%TR50\t%SD50\t%HWP\t%SQ5_SAMPLES\t%SQ10_SAMPLES\t%SQ20_SAMPLES\t%SQ30_SAMPLES\t%SQ40_SAMPLES\t%SQ50_SAMPLES\t%SQ60_SAMPLES\t%SQ70_SAMPLES\t%AD2_SAMPLES\t%AD3_SAMPLES\t%AD4_SAMPLES\t%AD5_SAMPLES\t%AD1_P_ALT_SAMPLES\t%AD1_P_REF_SAMPLES\t%AD1_M_ALT_SAMPLES\t%AD1_M_REF_SAMPLES\t%AD1_N_ALT_SAMPLES\t%AD1_N_REF_SAMPLES\t%HET_SAMPLES\t%HOMALT_SAMPLES\t%ZERO_COV_SAMPLES\t%gnomadV4_AF\t%gnomadV4_PC\t%gnomadV4_POPMAX_AF\t%gnomadV4_name\t%GNOCCHI\t%FB_PR_ENH_M\t%FB_PR_ENH_F\t%FANTOM_ENH\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%BIOTYPE\t%DISTANCE\t%TSSDistance\t%AF_FOUNDER_LR\t%AF_OFFSPRING_LR\t%AF_FOUNDER_IL\t%AF_OFFSPRING_IL\n' $VCF_IN >> $VEP_OUT
date
