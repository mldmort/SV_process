import numpy
import argparse


def get_args():
	# Create ArgumentParser object
	parser = argparse.ArgumentParser(description='Setting minimal arguments for setting the score table.')
	
	# Add arguments
	parser.add_argument('--pli_thr', help='Significant pLI threshold.', default=.99, type=float)
	parser.add_argument('--loeuf_thr', help='Significant LOEUF threshold.', default=.37, type=float)
	parser.add_argument('--fdr_thr', help='Significant FDR threshold.', default=.05, type=float)
	parser.add_argument('--file_out', help='Output file name.', default='scores.tsv', type=str)
	
	# Parse the command-line arguments
	args = parser.parse_args()

	return args

args = get_args()
sig_pli_thr = args.pli_thr
sig_loeuf_thr = args.loeuf_thr
sig_fdr_thr = args.fdr_thr
file_out = args.file_out
print(f'sig_pli_thr: {sig_pli_thr}')
print(f'sig_loeuf_thr: {sig_loeuf_thr}')
print(f'sig_fdr_thr: {sig_fdr_thr}')
print(f'file_out: {file_out}')

# Header
fh = open(file_out, 'w')
fh.write("CHROM\tPOS\tEND\tID\tSVTYPE\t")
fh.write("GENES_PLI\tMAX_PLI\tMAX_PLI_GENE\tSIG_PLIS\tSIG_PLI_GENES\tSIG_PLI_ENS\tSIG_PLI_BIOTYPE\tSIG_PLI_IMPACT\tSIG_PLI_CSQ\tSIG_PLI_DIST\tSIG_PLI_TSS_DIST\t")
fh.write("GENES_LOEUF\tMIN_LOEUF\tMIN_LOEUF_GENE\tSIG_LOEUFS\tSIG_LOEUF_GENES\tSIG_LOEUF_ENS\tSIG_LOEUF_BIOTYPE\tSIG_LOEUF_IMPACT\tSIG_LOEUF_CSQ\tSIG_LOEUF_DIST\tSIG_LOEUF_TSS_DIST\t")
fh.write("GENES_FDR_ASD\tMIN_FDR_ASD\tMIN_FDR_ASD_GENE\tSIG_FDR_ASDS\tSIG_FDR_ASD_GENES\tSIG_FDR_ASD_ENS\tSIG_FDR_ASD_BIOTYPE\tSIG_FDR_ASD_IMPACT\tSIG_FDR_ASD_CSQ\tSIG_FDR_ASD_DIST\tSIG_FDR_ASD_TSS_DIST\t")
fh.write("GENES_FDR_DD\tMIN_FDR_DD\tMIN_FDR_DD_GENE\tSIG_FDR_DDS\tSIG_FDR_DD_GENES\tSIG_FDR_DD_ENS\tSIG_FDR_DD_BIOTYPE\tSIG_FDR_DD_IMPACT\tSIG_FDR_DD_CSQ\tSIG_FDR_DD_DIST\tSIG_FDR_DD_TSS_DIST\t")
fh.write("GENES_FDR_NDD\tMIN_FDR_NDD\tMIN_FDR_NDD_GENE\tSIG_FDR_NDDS\tSIG_FDR_NDD_GENES\tSIG_FDR_NDD_ENS\tSIG_FDR_NDD_BIOTYPE\tSIG_FDR_NDD_IMPACT\tSIG_FDR_NDD_CSQ\tSIG_FDR_NDD_DIST\tSIG_FDR_NDD_TSS_DIST\t")
fh.write("GENES_S_HET\n")

file_tada = '/expanse/projects/sebat1/miladm/UCSD/resources/ASD_TADA_2022/fdr_pvals_tada.tsv'
#file_gnomadv4 = '/expanse/projects/sebat1/miladm/UCSD/resources/gnomAD/v4/Constraint/gnomad.v4.0.constraint_metrics.tsv'
file_gnomadv4 = '/expanse/projects/sebat1/miladm/UCSD/resources/gnomAD/v4/Constraint/gnomad.v4.0.constraint_metrics_chrX_chrY_v3.tsv'
file_s_het = '/expanse/projects/sebat1/miladm/UCSD/resources/S_HET/Supplementary_table_1.tsv'
### for TSCC
#file_tada =         '/tscc/projects/ps-sebat1/miladm/UCSD/resources/ASD_TADA_2022/fdr_pvals_tada.tsv'
##file_gnomadv4 =    '/tscc/projects/ps-sebat1/miladm/UCSD/resources/gnomAD/v4/Constraint/gnomad.v4.0.constraint_metrics.tsv'
#file_gnomadv4 =     '/tscc/projects/ps-sebat1/miladm/UCSD/resources/gnomAD/v4/Constraint/gnomad.v4.0.constraint_metrics_chrX_chrY_v3.tsv'
#file_s_het =        '/tscc/projects/ps-sebat1/miladm/UCSD/resources/S_HET/Supplementary_table_1.tsv'

file_vep = 'vep.tsv'

ensg_s_het = {}
with open(file_s_het, 'r') as f:
	#feild_num	name
	#1	ensg
	#7 post_mean
	f.readline()
	for line in f:
		line_parts = line.rstrip().split("\t")
		ensg = line_parts[0]
		post_mean = line_parts[6]
		ensg_s_het[ensg] = post_mean

gene_pli = {}
gene_loeuf = {}
with open(file_gnomadv4, 'r') as f:
	#feild_num	name
	#1	gene
	#17 lof.pLI
	#21 lof.oe_ci.upper
	f.readline()
	for line in f:
		line_parts = line.rstrip().split("\t")
		gene = line_parts[0]
		pli = line_parts[16]
		loeuf = line_parts[20]
		#gene_pli[gene] = pli
		if pli != 'NA':
			if gene not in gene_pli:
				gene_pli[gene] = pli
			elif float(pli) > float(gene_pli[gene]):
				gene_pli[gene] = pli
		#gene_loeuf[gene] = loeuf
		if loeuf != 'NA':
			if gene not in gene_loeuf:
				gene_loeuf[gene] = loeuf
			elif float(loeuf) < float(gene_loeuf[gene]):
				gene_loeuf[gene] = loeuf
		
gene_fdr_asd = {}
gene_fdr_dd  = {}
gene_fdr_ndd = {}
with open(file_tada, 'r') as f:
	f.readline()
	for line in f:
		gene, fdr_tada_asd, fdr_tada_dd, fdr_tada_ndd, p_tada_asd, p_tada_dd, p_tada_ndd = line.rstrip().split("\t")
		gene_fdr_asd[gene] = fdr_tada_asd
		gene_fdr_dd[gene] = fdr_tada_dd
		gene_fdr_ndd[gene] = fdr_tada_ndd

# Get max pLI, min LOEUF, min FDRs 
with open(file_vep, 'r') as f:
	header = f.readline().rstrip().split("\t")
	i_CHROM = header.index('CHROM')
	i_POS = header.index('POS')
	i_END = header.index('END')
	i_ID = header.index('ID')
	i_SVTYPE = header.index('SVTYPE')
	i_Consequence = header.index('Consequence')
	i_IMPACT = header.index('IMPACT')
	i_SYMBOL = header.index('SYMBOL')
	i_GENE = header.index('GENE')
	i_BIOTYPE = header.index('BIOTYPE')
	i_DISTANCE = header.index('DISTANCE')
	i_TSSDistance = header.index('TSSDistance')
	
	for line in f:
		linesplit = line.split("\t")
		chrom, pos, END, ID, svtype = linesplit[i_CHROM], linesplit[i_POS], linesplit[i_END], linesplit[i_ID], linesplit[i_SVTYPE]
		Consequence = linesplit[i_Consequence]
		IMPACT = linesplit[i_IMPACT]
		SYMBOL = linesplit[i_SYMBOL]
		GENE = linesplit[i_GENE]
		BIOTYPE = linesplit[i_BIOTYPE]
		DISTANCE = linesplit[i_DISTANCE]
		TSSDistance = linesplit[i_TSSDistance]
		
		genes = SYMBOL.split(",")
		genes_ens = GENE.split(",")
		biotypes = BIOTYPE.split(",")
		dists = DISTANCE.split(",")
		tss_dists = TSSDistance.split(",")
		csqs = Consequence.split(",")
		impacts = IMPACT.split(",")

		gene_plis = []
		gene_loeufs = []
		gene_fdr_asds = []
		gene_fdr_dds = []
		gene_fdr_ndds = []
		gene_s_het = []

		sig_plis = []
		sig_pli_genes = []
		sig_pli_genes_ens = []
		sig_pli_biotypes = []
		sig_pli_impacts = []
		sig_pli_csqs = []
		sig_pli_dists = []
		sig_pli_tss_dists = []

		sig_loeufs = []
		sig_loeuf_genes = []
		sig_loeuf_genes_ens = []
		sig_loeuf_biotypes = []
		sig_loeuf_impacts = []
		sig_loeuf_csqs = []
		sig_loeuf_dists = []
		sig_loeuf_tss_dists = []

		sig_fdr_asds = []
		sig_fdr_asd_genes = []
		sig_fdr_asd_genes_ens = []
		sig_fdr_asd_biotypes = []
		sig_fdr_asd_impacts = []
		sig_fdr_asd_csqs = []
		sig_fdr_asd_dists = []
		sig_fdr_asd_tss_dists = []

		sig_fdr_dds = []
		sig_fdr_dd_genes = []
		sig_fdr_dd_genes_ens = []
		sig_fdr_dd_biotypes = []
		sig_fdr_dd_impacts = []
		sig_fdr_dd_csqs = []
		sig_fdr_dd_dists = []
		sig_fdr_dd_tss_dists = []

		sig_fdr_ndds = []
		sig_fdr_ndd_genes = []
		sig_fdr_ndd_genes_ens = []
		sig_fdr_ndd_biotypes = []
		sig_fdr_ndd_impacts = []
		sig_fdr_ndd_csqs = []
		sig_fdr_ndd_dists = []
		sig_fdr_ndd_tss_dists = []

		max_pli     = -1
		max_pli_gene = "."
		min_loeuf   = 1000
		min_loeuf_gene = "."
		min_fdr_asd = 1000
		min_fdr_asd_gene = "."
		min_fdr_dd  = 1000
		min_fdr_dd_gene = "."
		min_fdr_ndd = 1000
		min_fdr_ndd_gene = "."

		for ensg in genes_ens:
			if ensg in ensg_s_het:
				gene_s_het.append(ensg_s_het[ensg])
			else:
				gene_s_het.append('.')
				
		for i_g, gene in enumerate(genes):
			if gene in gene_pli: 
				gene_plis.append(gene_pli[gene])
				#if gene_pli[gene] != "NA":
				if float(gene_pli[gene]) > max_pli: 
					max_pli = float(gene_pli[gene])
					max_pli_gene = gene
				if float(gene_pli[gene]) >= sig_pli_thr: 
					sig_pli_genes.append(gene)
					sig_plis.append(gene_pli[gene])
					sig_pli_genes_ens.append(genes_ens[i_g])
					sig_pli_biotypes.append(biotypes[i_g])
					sig_pli_impacts.append(impacts[i_g])
					sig_pli_csqs.append(csqs[i_g])
					sig_pli_dists.append(dists[i_g])
					sig_pli_tss_dists.append(tss_dists[i_g])
			else: gene_plis.append(".")
			
			if gene in gene_loeuf: 
				gene_loeufs.append(gene_loeuf[gene])
				#if gene_loeuf[gene] != "NA":
				if float(gene_loeuf[gene]) < min_loeuf: 
					min_loeuf = float(gene_loeuf[gene])
					min_loeuf_gene = gene
				if float(gene_loeuf[gene]) <= sig_loeuf_thr: 
					sig_loeuf_genes.append(gene)
					sig_loeufs.append(gene_loeuf[gene])
					sig_loeuf_genes_ens.append(genes_ens[i_g])
					sig_loeuf_biotypes.append(biotypes[i_g])
					sig_loeuf_impacts.append(impacts[i_g])
					sig_loeuf_csqs.append(csqs[i_g])
					sig_loeuf_dists.append(dists[i_g])
					sig_loeuf_tss_dists.append(tss_dists[i_g])
			else: gene_loeufs.append(".")
			
			if gene in gene_fdr_asd: 
				gene_fdr_asds.append(gene_fdr_asd[gene])
				if float(gene_fdr_asd[gene]) < min_fdr_asd: 
					min_fdr_asd = float(gene_fdr_asd[gene])
					min_fdr_asd_gene = gene
					if float(gene_fdr_asd[gene]) <= sig_fdr_thr: 
						sig_fdr_asd_genes.append(gene)
						sig_fdr_asds.append(gene_fdr_asd[gene])
						sig_fdr_asd_genes_ens.append(genes_ens[i_g])
						sig_fdr_asd_biotypes.append(biotypes[i_g])
						sig_fdr_asd_impacts.append(impacts[i_g])
						sig_fdr_asd_csqs.append(csqs[i_g])
						sig_fdr_asd_dists.append(dists[i_g])
						sig_fdr_asd_tss_dists.append(tss_dists[i_g])
			else: gene_fdr_asds.append(".")
			
			if gene in gene_fdr_dd: 
				gene_fdr_dds.append(gene_fdr_dd[gene])
				if float(gene_fdr_dd[gene]) < min_fdr_dd: 
					min_fdr_dd = float(gene_fdr_dd[gene])
					min_fdr_dd_gene = gene
					if float(gene_fdr_dd[gene]) <= sig_fdr_thr: 
						sig_fdr_dd_genes.append(gene)
						sig_fdr_dds.append(gene_fdr_dd[gene])
						sig_fdr_dd_genes_ens.append(genes_ens[i_g])
						sig_fdr_dd_biotypes.append(biotypes[i_g])
						sig_fdr_dd_impacts.append(impacts[i_g])
						sig_fdr_dd_csqs.append(csqs[i_g])
						sig_fdr_dd_dists.append(dists[i_g])
						sig_fdr_dd_tss_dists.append(tss_dists[i_g])
			else: gene_fdr_dds.append(".")
			
			if gene in gene_fdr_ndd: 
				gene_fdr_ndds.append(gene_fdr_ndd[gene])
				if float(gene_fdr_ndd[gene]) < min_fdr_ndd: 
					min_fdr_ndd = float(gene_fdr_ndd[gene])
					min_fdr_ndd_gene = gene
					if float(gene_fdr_ndd[gene]) <= sig_fdr_thr: 
						sig_fdr_ndd_genes.append(gene)
						sig_fdr_ndds.append(gene_fdr_ndd[gene])
						sig_fdr_ndd_genes_ens.append(genes_ens[i_g])
						sig_fdr_ndd_biotypes.append(biotypes[i_g])
						sig_fdr_ndd_impacts.append(impacts[i_g])
						sig_fdr_ndd_csqs.append(csqs[i_g])
						sig_fdr_ndd_dists.append(dists[i_g])
						sig_fdr_ndd_tss_dists.append(tss_dists[i_g])
			else: gene_fdr_ndds.append(".")

		fh.write(f"{chrom}\t{pos}\t{END}\t{ID}\t{svtype}\t")
		fh.write(f"{','.join(gene_plis)}\t{max_pli}\t{max_pli_gene}\t{','.join(sig_plis)}\t{','.join(sig_pli_genes)}\t{','.join(sig_pli_genes_ens)}\t{','.join(sig_pli_biotypes)}\t{','.join(sig_pli_impacts)}\t{','.join(sig_pli_csqs)}\t{','.join(sig_pli_dists)}\t{','.join(sig_pli_tss_dists)}\t")
		fh.write(f"{','.join(gene_loeufs)}\t{min_loeuf}\t{min_loeuf_gene}\t{','.join(sig_loeufs)}\t{','.join(sig_loeuf_genes)}\t{','.join(sig_loeuf_genes_ens)}\t{','.join(sig_loeuf_biotypes)}\t{','.join(sig_loeuf_impacts)}\t{','.join(sig_loeuf_csqs)}\t{','.join(sig_loeuf_dists)}\t{','.join(sig_loeuf_tss_dists)}\t")
		fh.write(f"{','.join(gene_fdr_asds)}\t{min_fdr_asd}\t{min_fdr_asd_gene}\t{','.join(sig_fdr_asds)}\t{','.join(sig_fdr_asd_genes)}\t{','.join(sig_fdr_asd_genes_ens)}\t{','.join(sig_fdr_asd_biotypes)}\t{','.join(sig_fdr_asd_impacts)}\t{','.join(sig_fdr_asd_csqs)}\t{','.join(sig_fdr_asd_dists)}\t{','.join(sig_fdr_asd_tss_dists)}\t")
		fh.write(f"{','.join(gene_fdr_dds)}\t{min_fdr_dd}\t{min_fdr_dd_gene}\t{','.join(sig_fdr_dds)}\t{','.join(sig_fdr_dd_genes)}\t{','.join(sig_fdr_dd_genes_ens)}\t{','.join(sig_fdr_dd_biotypes)}\t{','.join(sig_fdr_dd_impacts)}\t{','.join(sig_fdr_dd_csqs)}\t{','.join(sig_fdr_dd_dists)}\t{','.join(sig_fdr_dd_tss_dists)}\t")
		fh.write(f"{','.join(gene_fdr_ndds)}\t{min_fdr_ndd}\t{min_fdr_ndd_gene}\t{','.join(sig_fdr_ndds)}\t{','.join(sig_fdr_ndd_genes)}\t{','.join(sig_fdr_ndd_genes_ens)}\t{','.join(sig_fdr_ndd_biotypes)}\t{','.join(sig_fdr_ndd_impacts)}\t{','.join(sig_fdr_ndd_csqs)}\t{','.join(sig_fdr_ndd_dists)}\t{','.join(sig_fdr_ndd_tss_dists)}\t")
		fh.write(f"{','.join(gene_s_het)}\n")

fh.close()
