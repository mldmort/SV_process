# you should be in an conda/micromamba environment with pandas
import pandas as pd
import sys

vep_file = 'vep.tsv'
#'variants_selected.tsv'
score_file = 'scores.tsv'
#'scores_selected.tsv'
dup_inv_lof_file = 'dup_inv_lof.bed'
sd_art_file = 'sd_artifact_svs.bed'
output_file = 'combined_table.tsv'

df_v = pd.read_table(vep_file, sep='\t', header=0)
df_v['key'] = df_v['CHROM'] + "_" + df_v['POS'].astype(str) + "_" + df_v['END'].astype(str) + "_" + df_v['ID']
print('vep table:')
print(df_v)

df_s = pd.read_table(score_file, sep='\t', header=0)
df_s['key'] = df_s['CHROM'] + "_" + df_s['POS'].astype(str) + "_" + df_s['END'].astype(str) + "_" + df_s['ID']
print('score table:')
print(df_s)

df_l = pd.read_table(dup_inv_lof_file, sep='\t', header=0)
print('dup/inv lof table:')
print(df_l)

df_v_l = pd.merge(df_v, df_l, on='key', how='outer')
df_v_l['LOF_DUP_INV'] = df_v_l['LOF_DUP_INV'].fillna(0).astype(int)
df_v_l['LOF_DUP_INV_GENES'] = df_v_l['LOF_DUP_INV_GENES'].fillna("").astype(str)
print('df_v_l:')
print(df_v_l)
print(f"sum df_v_l['LOF_DUP_INV']: {df_v_l['LOF_DUP_INV'].sum()}")
#print(df_v_l.loc[df_v_l.LOF_DUP_INV != 0])

df_sd = pd.read_table(sd_art_file, sep='\t', header=0)
print('SD artifact table:')
print(df_sd)

df_v_l_sd = pd.merge(df_v_l, df_sd, on='key', how='outer')
df_v_l_sd['SD_ART'] = df_v_l_sd['SD_ART'].fillna(0).astype(int)
df_v_l_sd['SD_SEG1'] = df_v_l_sd['SD_SEG1'].fillna("").astype(str)
df_v_l_sd['SD_SEG2'] = df_v_l_sd['SD_SEG2'].fillna("").astype(str)
print(df_v_l_sd)
print(f"sum df_v_l_sd['SD_ART']: {df_v_l_sd['SD_ART'].sum()}")

df = pd.merge(df_v_l_sd, df_s, on='key', how='inner')
print('merged table:')
print(df)

print('reduce redundant variant annotations...')
# for variant columns, with redundant ','-separated repeated values, reduce the column to one.
var_cols = ['gnomadV4_AF', 'gnomadV4_PC', 'gnomadV4_POPMAX_AF', 'gnomadV4_name', 'GENCODE', 'GNOCCHI', 'FB_PR_ENH_M', 'FB_PR_ENH_F', 'FANTOM_ENH']

for var_col in var_cols:
	#print("df['GENCODE']:")
	#print(df['GENCODE'])
	df[var_col] = df[var_col].str.split(',',expand=False).str.get(0)
	#print("df['GENCODE']:")
	#print(df['GENCODE'])

print('prioritize GENCODE/ENCODE annotations...')
# for GENCODE column, with functional annotations output the highest priority annotation.
# the rank of annotations is according to the following list, with decreasing priority.
#func_el_rank = ['start_codon', 'stop_codon', 'stop_codon_redefined_as_selenocysteine', 'CDS', 'TSS', 'exon', 'five_prime_UTR', 'three_prime_UTR', 'DNase-H3K4me3', 'PLS', 'pELS', 'dELS', 'CTCF-only']
func_el_rank = ['start_codon', 'stop_codon', 'stop_codon_redefined_as_selenocysteine', 'CDS', 'TSS', 'five_prime_UTR', 'three_prime_UTR', 'exon', 'DNase-H3K4me3', 'PLS', 'pELS', 'dELS', 'CTCF-only']
func2rank = {el: i for i, el in enumerate(func_el_rank)}

def get_priority(func_str):
	if func_str == '.':
		return '.'
	func_list = func_str.split('&')
	func_rank_list = [(f, func2rank[f]) for f in func_list]
	func_rank_list_sorted = sorted(func_rank_list, key=lambda x: x[1], reverse=False) # sorts from small to large
	return func_rank_list_sorted[0][0]

df['GENCODE'] = df['GENCODE'].apply(get_priority)

print('clean up gnomadV4 annotations...')
# for gnomAD SV annotations check if reciprocal overlap is greated than the threshold
# and choose the highest overlap data.
def get_gnomad_vals(row):
	if row['gnomadV4_PC'] == '.':
		return '.', '.', '.', '.'
	rec_ov_thr = 80
	pc_idx_list = [ (float(x), i) for i, x in enumerate(row['gnomadV4_PC'].split('&')) ]
	pc_idx_list_sorted = sorted(pc_idx_list, key=lambda x: x[0], reverse=True) # sorts from large to small
	if pc_idx_list_sorted[0][0] >= rec_ov_thr:
		gnomadV4_AF 		= row['gnomadV4_AF'].split('&')
		gnomadV4_POPMAX_AF 	= row['gnomadV4_POPMAX_AF'].split('&')
		gnomadV4_name 		= row['gnomadV4_name'].split('&')
		idx = pc_idx_list_sorted[0][1]
		ret = gnomadV4_AF[idx], pc_idx_list_sorted[0][0], gnomadV4_POPMAX_AF[idx], gnomadV4_name[idx]
	else:
		ret = '.', '.', '.', '.'
	return ret

df[['gnomadV4_AF', 'gnomadV4_PC' ,'gnomadV4_POPMAX_AF', 'gnomadV4_name']] = df.apply(get_gnomad_vals, axis=1, result_type='expand')

print('choose maximum Gnocchi value...')
# choose maximum Gnocchi value for the variant
def get_gnocchi(gnocchi_str):
	if gnocchi_str == '.':
		return '.'
	return max([float(x) for x in gnocchi_str.split('&')])

df['GNOCCHI'] = df['GNOCCHI'].apply(lambda gno: get_gnocchi(gno))

print('reduce redundant FB and FANTOM annotations...')
# collapse repeated annotation for the columns below
df['FB_PR_ENH_M'] = df['FB_PR_ENH_M'].apply(lambda x: '&'.join(set(x.split('&'))))
df['FB_PR_ENH_F'] = df['FB_PR_ENH_F'].apply(lambda x: '&'.join(set(x.split('&'))))
df['FANTOM_ENH'] = df['FANTOM_ENH'].apply(lambda x: '&'.join(set(x.split('&'))))

print('rename and drop extra columns...')
rename_mapper = {'CHROM_x': 'CHROM', 'POS_x': 'POS', 'END_x': 'END', 'ID_x': 'ID', 'SVTYPE_x': 'SVTYPE'}
df.rename(columns=rename_mapper, inplace=True)

drop_cols = ['key' , 'CHROM_y', 'POS_y', 'END_y', 'ID_y', 'SVTYPE_y']
df.drop(drop_cols, inplace=True, axis=1)

meta_file = '/expanse/projects/sebat1/miladm/UCSD/LONG_READ_COHORT/REACH_sample_info.tsv'
df_meta = pd.read_table(meta_file, sep='\t', header=0)
print(df_meta)

aff_dict = {}
for sample, aff in zip(df_meta['Sample_ID'].tolist(), df_meta['Affected'].tolist()):
	aff_dict[sample] = aff
#print(aff_dict)

def get_case(Ser):
	Ser_split = Ser.str.split(',')
	Ser_new = Ser_split.apply(lambda x: [aff_dict[xx.strip('_IL').strip('_PB').strip('_ONT')] if xx.strip('_IL').strip('_PB').strip('_ONT') in aff_dict != "." else "." for xx in x])
	return Ser_new.str.join(sep=',')

cols = ["SQ5_SAMPLES" ,"SQ10_SAMPLES" ,"SQ20_SAMPLES" ,"SQ30_SAMPLES" ,"SQ40_SAMPLES" ,"AD2_SAMPLES" ,"AD3_SAMPLES" ,"AD4_SAMPLES" ,"AD5_SAMPLES" ,"ZERO_COV_SAMPLES"]
df[["case_"+col for col in cols]] = df[cols].apply(get_case, axis=0)

print('write the final table...')
df.to_csv(output_file, sep='\t', header=True, index=False)

