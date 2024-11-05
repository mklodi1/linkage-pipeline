import re
import csv
import gzip
from collections import defaultdict
from datetime import date

today = date.today()


gtex_data = {}
brainspan_data = {}
hdbr_data = {}
impc_knockout = {}
gnomad = {}
known_NDD = []
LI_SNV_genes = []
LI_CNV_genes = []
LI_SV_genes = []
LI_all_genes = []
final_data = []
HGNC_genes = {}
reg_genes = {}
SNP_dom_fam = {}
SNP_rec_fam = {}
SV_fam = {}
CNV_fam = {}
SNP_unique = {}
SNP_unique = defaultdict(list)
CNV_unique = {}
SV_unique = {}
all_unique = {}
SNV_count = defaultdict(list)
SV_var_count = {}
CNV_var_count = {}
pvalue = {}

import sys
results_path = sys.argv[1]
dom_var_res = results_path + "/SNP_Analysis/Results/LI_ASD_dom_var.txt"
rec_var_res = results_path + "/SNP_Analysis/Results/LI_ASD_dom_var.txt"
dom_genes = results_path + "/SNP_Analysis/Results/LI_dom_genes_fam.txt"
rec_genes = results_path + "/SNP_Analysis/Results/LI_rec_genes_fam.txt"
SNP_genes = results_path + "/SNP_Analysis/Results/LI_SNP_genes.txt"

CNV_genes = results_path + "/CNV_Analysis/Results/LI_CNV_genes.txt"
CNV_gene_fam = results_path + "/CNV_Analysis/Results/LI_CNV_gene_fam.txt"

SV_non_benign = results_path + "/SV_Analysis/Results/LI_SV_inh_non_benign.txt"
SV_nb_genes = results_path + "/SV_Analysis/Results/LI_SV_genes_filt.txt"

gene_table = results_path + "/Gene_Analysis/LI_gene_table" 

# create dictionary of transcript ID and gene
with open('/lab01/Projects/Mudassir_Projects/Data/HGNC/HGNC_Data.txt') as f9:
	for line in f9:
		cols = line.split('\t')
		gene = cols[0].strip()
		transcript = cols[2].strip()
		HGNC_genes[transcript] = gene	

# list of known NDD genes
with open('/lab01/Projects/Rohan_Projects/CNV_Project/2023/CNV_Pipeline/Data/NDD_genes.tsv') as f4:
	for line in f4:
		cols = line.split('\t')
		if 'gene' in line:
			continue
		else:
			known_NDD.append(cols[0])
			
# GTEx Max TPM value
with gzip.open('/lab02/Data_Raw/Xiaolong/HumanCommon/GTEx_v8/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',mode='rt') as f2:
	for line in f2:
		cols = line.split('\t')
		if 'Name' in cols:
			continue
		else:
			if len(cols) > 2:
				gene_name = cols[1]
				data = cols[9:21] # brain related expression data
				tpm_vals = [eval(i) for i in data]
				max_tpm_val = max(tpm_vals)
				gtex_data[gene_name] = max_tpm_val

# BrainSpan Max TPM value
with open('/lab02/Data_Raw/Xiaolong/HumanCommon/Expression/BrainSpan/20190508BrainSpanGeneExpression.csv') as f3:
	for line in f3:
		if 'gene' in line:
			continue
		cols = line.split('\t')
		if len(cols) > 2:
			gene_name = cols[0]
			data = cols[2:]
			bs_tpm_vals = [eval(i) for i in data]
			bs_max_tpm_val = max(bs_tpm_vals)
			brainspan_data[gene_name] = bs_max_tpm_val	

# HDBR Max TPM value
with open('/lab02/Data_Raw/Xiaolong/HumanCommon/Expression/brainDevelopment/tpms.tsv') as f7:
	for line in f7:
		if '#' in line:
			continue
		elif 'basal' in line:
			continue
		else:
			cols = line.split('\t')
			gene_name = cols[1]
			data = cols[2:]
			for elem in data:
				elem = elem.strip('\n')
				final_data.append(elem)
			while('' in final_data):
				final_data.remove('')
			hdbr_tpm_vals = [eval(i) for i in final_data]
			hdbr_max_tpm_val = max(hdbr_tpm_vals)
			hdbr_data[gene_name] = hdbr_max_tpm_val
		final_data = []


# IMPC knockout phenotype
with gzip.open('/nfs/xing_lab/DataSets/IMPC/release-20.1/results/phenotypeHitsPerGene.csv.gz',mode='rt') as f5:
	impc = csv.reader(f5)
	for line in impc:
		if 'Gene Symbol' in line:
			continue
		else:
			gene_name = line[0]
			gene_name = gene_name.upper()
			phenotype = line[3]
			impc_knockout[gene_name] = phenotype


# LOEUF, Missense Z-score, pLI from gnomAD
with open('/lab02/Data_Raw/Xiaolong/HumanCommon/gnomAD/gnomad.v2.1.1.lof_metrics.by_gene.txt') as f8:
	for line in f8:
		if 'Gene' in line:
			continue
		else:
			cols = line.split('\t')
			gene_name = cols[0]
			pli = cols[20]
			oe_lof = cols[23]
			loeuf = cols[29]
			misz = cols[32]
			
			gene_data = [pli, oe_lof, loeuf, misz]
			gnomad[gene_name] = gene_data
		
with open (SNP_genes) as f:
	for line in f:
		cols = line.split('\t')
		gene = cols[0]
		transcript = cols[1].strip()
		reg_genes[gene] = transcript
		LI_all_genes.append(gene)
		LI_SNV_genes.append(gene)

with open (CNV_genes) as f:
	for line in f:
		cols = line.split('\t')
		gene = cols[0].strip()
		count = str(cols[1])
		#transcript = line[1]
		#HGNC_gene = transcript_gene[transcript]
		LI_all_genes.append(gene)
		LI_CNV_genes.append(gene)
		#HGNC_genes.append(HGNC_gene)
		CNV_var_count[gene] = count

with open (SV_nb_genes) as f:
	for line in f:
		cols = line.split('\t')
		gene = cols[0]
		transcript = cols[1].strip()
		count = str(cols[2]).strip()
		reg_genes[gene] = transcript
		LI_all_genes.append(gene)
		LI_SV_genes.append(gene)
		SV_var_count[gene] = count

with open (dom_genes) as f:
	for line in f:
		cols = line.split('\t')
		gene = cols[0].strip()
		fam = cols[1].strip()
		samp = cols[2].strip()
		SNP_dom_fam[gene] = [fam,samp]
		fams = re.findall(r'\d+',fam)
		fams = list(set(fams))
		SNP_unique[gene] = fams
		if gene in all_unique.keys():
			current_fams = all_unique[gene]
			fams = current_fams + fams
			fams = list(set(fams))
			all_unique[gene] = fams
		else:
			all_unique[gene] = fams


with open (rec_genes) as f:
	for line in f:
		cols = line.split('\t')
		gene = cols[0].strip()
		fam = cols[1].strip()
		samp = cols[2].strip()
		SNP_rec_fam[gene] = [fam,samp]
		fams = re.findall(r'\d+',fam)
		fams = list(set(fams))
		if gene in SNP_unique.keys():
			dom_fams = SNP_unique[gene]
			fams = dom_fams + fams
			fams = list(set(fams))
			SNP_unique[gene] = fams
		else:
			SNP_unique[gene] = fams
		if gene in all_unique.keys():
			current_fams = all_unique[gene]
			fams = current_fams + fams
			fams = list(set(fams))
			all_unique[gene] = fams
		else:
			all_unique[gene] = fams

with open (CNV_gene_fam) as f:
	for line in f:
		cols = line.split('\t')
		gene = cols[0].strip()
		fam = cols[1].strip()
		samp = cols[2].strip()
		CNV_fam[gene] = [fam,samp]
		fams = re.findall(r'\d+',fam)
		fams = list(set(fams))
		CNV_unique[gene] = fams
		if gene in all_unique.keys():
			current_fams = all_unique[gene]
			fams = current_fams + fams
			fams = list(set(fams))
			all_unique[gene] = fams
		else:
			all_unique[gene] = fams


with open (SV_non_benign) as f:
	for line in f:
		cols = line.split('\t')
		gene = cols[0].strip()
		dom_fam = cols[1].strip()
		dom_samp = cols[2].strip()
		rec_fam = cols[3].strip()
		rec_samp = cols[4].strip()
		denovo_fam = cols[5].strip()
		denovo_samp = cols[6].strip()
		SV_fam[gene] = [dom_fam,dom_samp,rec_fam,rec_samp,denovo_fam,denovo_samp]
		fams = re.findall(r'\d+',dom_fam) + re.findall(r'\d+',rec_fam) + re.findall(r'\d+',denovo_fam)
		fams = list(set(fams))
		SV_unique[gene] = fams
		if gene in all_unique.keys():
			current_fams = all_unique[gene]
			fams = current_fams + fams
			fams = list(set(fams))
			all_unique[gene] = fams
		else:
			all_unique[gene] = fams

with open(dom_var_res) as f:
	for line in f:
		cols = line.split('\t')
		if cols[0] != '':
			gene = cols[1]
			pos = cols[15]
			pval = float(cols[7])
			pvalue[gene] = pval
			SNV_count[gene].append(pos)
		else:
			pos = cols[15]
			SNV_count[gene].append(pos)

with open(rec_var_res) as f:
	for line in f:
		cols = line.split('\t')
		if cols[0] != '':
			gene = cols[1]
			pos = cols[15]
			SNV_count[gene].append(pos)
		else:
			pos = cols[15]
			SNV_count[gene].append(pos) 

for gene in SNV_count.keys():
	SNV_count[gene] = list(set(SNV_count[gene]))

HGNC_manual = {'AGPHD1':'HYKK','C15orf38':'C15orf38-AP3S2','INSYN1':'C15orf59','TTLL13P':'TTLL13'}
print(SNV_count)

#print(all_unique)
#print(SV_unique['MYO9A'])
#print(all_unique['MYO9A'])
LI_all_genes = list(set(LI_all_genes))

# create final candidate gene table by gene
with open (gene_table, 'w') as f6:
	f6.write('Gene Name' + '\t' + 'Known_NDD' + '\t' + 'SNV' + '\t' + 'CNV' + '\t' + 'SV' + '\t' + 'GTEx' + '\t' + 'Brainspan' + '\t' + 'HDBR' + '\t' + 'pLI' + '\t' + 'oe_lof' + '\t' + 'LOEUF' + '\t' + 'mis_z_score' + '\t' + 'SNP_var_count' + '\t' + 'SNP_fam_count' + '\t' + 'SNP_dom_families' + '\t' + 'SNP_dom_samples' + '\t' + 'SNP_rec_families' + '\t' + 'SNP_rec_samples' + '\t' + 'CNV_var_count' + '\t' + 'CNV_fam_count' + '\t' + 'CNV_families' + '\t' + 'CNV_samples' + '\t' + 'SV_var_count' + '\t' + 'SV_fam_count' + '\t' + 'SV_dom_families' + '\t' + 'SV_dom_samples' + '\t' + 'SV_rec_families' + '\t' + 'SV_rec_samples' + '\t' + 'SV_denovo_families' + '\t' + 'SV_denovo_samples' + '\t' + 'Unique_families' + '\t' + 'Unique_fam_count' + '\t' + 'IMPC_Knockout' + '\n')
	for gene_name in LI_all_genes:
		f6.write(gene_name + '\t')
		if gene_name in known_NDD:
			f6.write('1' + '\t')
		else:
			f6.write('0' + '\t')
		'''
		if gene_name in LI_SNV_genes and gene_name in pvalue.keys():
				if pvalue[gene_name] < 0.000162:
					f6.write('1' + '\t')
				else:
					f6.write('0' + '\t')
		else:
                        f6.write('0' + '\t')
		'''
		if gene_name in LI_SNV_genes:
			f6.write('1' + '\t')
		else:
			f6.write('0' + '\t')

		if gene_name in LI_CNV_genes:
			f6.write('1' + '\t')
		else:
			f6.write('0' + '\t')

		if gene_name in LI_SV_genes:
			f6.write('1' + '\t')
		else:
			f6.write('0' + '\t')
		
		# check for GTEx
		if gene_name in gtex_data.keys():
			tpm_val = gtex_data[gene_name]
			if tpm_val > 5:
				tpm_val = str(tpm_val)
				f6.write('1' + '\t')
			else:
				tpm_val = str(tpm_val)
				f6.write('0' + '\t')
		else: 
			if gene_name in reg_genes.keys():
				print('GTEx', gene_name, reg_genes[gene_name])
				transcript = reg_genes[gene_name]
				if transcript in HGNC_genes.keys():
					HGNC_gene = HGNC_genes[transcript]
					if HGNC_gene in gtex_data.keys():
						tpm_val = gtex_data[HGNC_gene]
						print('GTEx', gene_name, HGNC_gene, transcript, tpm_val)
						if tpm_val > 5:
							tpm_val = str(tpm_val)
							f6.write('1' + '\t')
						else:
							f6.write('0' + '\t')
					else:
						f6.write('NA' + '\t')
				else:
					if gene_name in HGNC_manual.keys():
						gene_name = HGNC_manual[gene_name]
						print('gtex', gene_name)
						if gene_name in gtex_data.keys():
							tpm_val = gtex_data[gene_name]
							if tpm_val > 5:
								tpm_val = str(tpm_val)
								f6.write('1' + '\t')
							else:
								f6.write('0' + '\t')
						else:
							f6.write('NA' + '\t')
					else:
						f6.write('NA' + '\t')
			else:
				f6.write('NA' + '\t')

		# check for BrainSpan
		if gene_name in brainspan_data.keys():
			tpm_val = brainspan_data[gene_name]
			if tpm_val > 5:
				tpm_val = str(tpm_val)
				f6.write('1' + '\t')
			else:
				tpm_val = str(tpm_val)
				f6.write('0' + '\t')
		else:
			if gene_name in reg_genes.keys():
				print('BS', gene_name, reg_genes[gene_name])
				transcript = reg_genes[gene_name]
				if transcript in HGNC_genes.keys():
					HGNC_gene = HGNC_genes[transcript]
					if HGNC_gene in brainspan_data.keys():
						tpm_val = brainspan_data[HGNC_gene]
						print('BS', gene_name, HGNC_gene, transcript, tpm_val)
						if tpm_val > 5:
							tpm_val = str(tpm_val)
							f6.write('1' + '\t')
						else:
							f6.write('0' + '\t')
					else:
						f6.write('NA' + '\t')

				else:
					if gene_name in HGNC_manual.keys():
						gene_name = HGNC_manual[gene_name]
						print('bs', gene_name)
						if gene_name in brainspan_data.keys():
							tpm_val = brainspan_data[gene_name]
							if tpm_val > 5:
								tpm_val = str(tpm_val)
								f6.write('1' + '\t')
							else:
								f6.write('0' + '\t')
					else:
						f6.write('NA' + '\t')
			else:
				f6.write('NA' + '\t')
		
		# HDBR
		if gene_name in hdbr_data.keys():
			tpm_val = hdbr_data[gene_name]
			if tpm_val > 5:
				tpm_val = str(tpm_val)
				f6.write('1' + '\t')
			else:
				tpm_val = str(tpm_val)
				f6.write('0' + '\t')
		else:
			if gene_name in reg_genes.keys():
				#print('HDBR', gene_name, reg_genes[gene_name])
				transcript = reg_genes[gene_name]
				if transcript in HGNC_genes.keys():
					HGNC_gene = HGNC_genes[transcript]
					if HGNC_gene in hdbr_data.keys():
						tpm_val = hdbr_data[HGNC_gene]
						#print('HDBR', gene_name, HGNC_gene, transcript, tpm_val)
						if tpm_val > 5:
							tpm_val = str(tpm_val)
							f6.write('1' + '\t')
						else:
							f6.write('0' + '\t')
					else:
						f6.write('NA' + '\t')

				else:
					if gene_name in HGNC_manual.keys():
						gene_name = HGNC_manual[gene_name]
						print('hdbr', gene_name)
						if gene_name in hdbr_data.keys():
							print('hdbr', gene_name)
							tpm_val = hdbr_data[gene_name]
							if tpm_val > 5:
								tpm_val = str(tpm_val)
								f6.write('1' + '\t')
							else:
								f6.write('0' + '\t')
						else:
							f6.write('NA' + '\t')
					else:
						f6.write('NA' + '\t')
			else:
				f6.write('NA' + '\t')
		
		# gnomad
		if gene_name in gnomad.keys():
			gene_data = gnomad[gene_name]
			for elem in gene_data:
				if elem == '':
					f6.write('NA' + '\t')
				else:
					f6.write(elem + '\t')
		else:
			if gene_name in reg_genes.keys():
				print('gnomad', gene_name, reg_genes[gene_name])
				transcript = reg_genes[gene_name]
				if transcript in HGNC_genes.keys():
					HGNC_gene = HGNC_genes[transcript]
					print('gnomad', gene_name, HGNC_gene, transcript, tpm_val)
					if HGNC_gene in gnomad.keys():
						gene_data = gnomad[HGNC_gene]
						for elem in gene_data:
							if elem == '':
								f6.write('NA' + '\t')
							else:
								f6.write(elem + '\t')
					else:
						f6.write('NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t')

				else:
					if gene_name in HGNC_manual.keys():
						gene_name = HGNC_manual[gene_name]
						if gene_name in gnomad.keys():
							gene_data = gnomad[gene_name]
							print('gnomad', gene_name)
							for elem in gene_data:
								if elem == '':
									f6.write('NA' + '\t')
								else:
									f6.write(elem + '\t')
					else:
						f6.write('NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t')
			else:
				f6.write('NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t')
		
		# SNV Variant Count
		if gene_name in SNV_count.keys():
			var_count = len(SNV_count[gene_name])
			f6.write(str(var_count) + '\t')
		else:
			f6.write('0' + '\t')

		# SNV Count
		if gene_name in SNP_unique.keys():
			f6.write(str(len(SNP_unique[gene_name])) + '\t')
		
		else:
			f6.write('0' + '\t')

		# SNV dominant
		if gene_name in SNP_dom_fam.keys():
			fam_info = SNP_dom_fam[gene_name]
			for i in fam_info:
				f6.write(i + '\t')
		else:
			f6.write('NA' + '\t' + 'NA' + '\t')
		
		# SNV recessive
		if gene_name in SNP_rec_fam.keys():
			fam_info = SNP_rec_fam[gene_name]
			for i in fam_info:
				f6.write(i + '\t')
		else:
			f6.write('NA' + '\t' + 'NA' + '\t')
		
		# CNV Variant Count
		if gene_name in CNV_var_count.keys():
			f6.write(CNV_var_count[gene_name] + '\t')
		else:
			f6.write('0' + '\t')

		# CNV count
		if gene_name in CNV_unique.keys():
			f6.write(str(len(CNV_unique[gene_name])) + '\t')

		else:
			f6.write('0' + '\t')

		# CNV
		if gene_name in CNV_fam.keys():
			fam_info = CNV_fam[gene_name]
			if fam_info[0] == '':
				fam_info[0] = 'NA'
			if fam_info[1] == '':
				fam_info[1] = 'NA'
			for i in fam_info:
				f6.write(i + '\t')
		else:
			f6.write('NA' + '\t' + 'NA' + '\t')

		# SV Variant Count
		if gene_name in SV_var_count.keys():
			f6.write(str(SV_var_count[gene_name]) + '\t')
		else:
			f6.write('0' + '\t')

		# SV count
		if gene_name in SV_unique.keys():
			f6.write(str(len(SV_unique[gene_name])) + '\t')
		else:
			f6.write('0' + '\t')

		# SV
		if gene_name in SV_fam.keys():
			fam_info = SV_fam[gene_name]
			for i in fam_info:
				if i == '':
					i = 'NA'
				f6.write(i + '\t')
		else:
			f6.write('NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t')

		# Unique
		if gene_name in all_unique.keys():
			fams = all_unique[gene_name]
			f6.write('[')
			f6.write(','.join(fams))
			f6.write(']' + '\t')
			f6.write(str(len(all_unique[gene_name])) + '\t')

		else:
			f6.write('NA' + '\t' + '0' + '\t')

                # IMPC Knockout
		if gene_name in impc_knockout.keys():
			pheno = impc_knockout[gene_name]
			if pheno == '':
				f6.write('NA' + '\n')
			else:
				f6.write(pheno + '\n')
		else:
			f6.write('NA' + '\n')

