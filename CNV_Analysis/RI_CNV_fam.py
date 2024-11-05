import re

wgs_id = {}
gene_fam = {}
wgsids = []

import sys
results_path = sys.argv[1]
linkage_tsv = results_path + "/CNV_Analysis/Results/RI_CNV_linkage.tsv"
CNV_genes = results_path + "/CNV_Analysis/Results/RI_CNV_genes.txt"
CNV_gene_fam = results_path + "/CNV_Analysis/Results/RI_CNV_gene_fam.txt"

with open('/lab01/Projects/Mudassir_Projects/linkage/Data/cohort_summary_edited_2023-06-27.txt') as f:
	for line in f:
		if line == '\n':
			continue
		else:
			cols = line.split('\t')
			partid = cols[1]
			wgsid = cols[22]
			wgs_id[partid] = wgsid

with open(linkage_tsv) as f1:
	for line in f1:
		if 'chrom' in line:
			continue
		else:
			cols = line.split('\t')
			genes = cols[6].split('|')
			var_info = cols[7].split(',')
			var_info = var_info[3:]
			famid = re.findall(r'\d+',var_info[0])
			for c in var_info: 
				aff_status_id = c.split('|')
				print(aff_status_id)
				aff_status = aff_status_id[1]
				aff_status = aff_status.replace(')','')
				aff_status = aff_status.replace(']','')
				if 'affected' == aff_status:
					aff_sample = re.findall(r'\d+',c)
					if len(aff_sample) > 1:
						aff_sample = aff_sample[1]
					else:
						aff_sample = aff_sample[0]
					wgid = wgs_id[aff_sample]
					if wgid in wgsids:
						continue
					else:
						wgsids.append(wgid)
			family = famid[0]
			for gene in genes:
				gene_fam[gene] = [family,wgsids]
			wgsids = []

with open(CNV_gene_fam,'w') as f2:
	for gene in gene_fam:
		samps = gene_fam[gene]
		f2.write(gene + '\t' + samps[0] + '\t')
		for ind in samps[1]:
			f2.write(ind)
		f2.write('\n')
