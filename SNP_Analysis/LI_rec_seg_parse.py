"""
LI_rec_seg_parse.py

> under recessive inheritance pattern, identifies samples and families containing each language impairment SNP candidate risk gene

"""
import sys
import re
import csv
from collections import defaultdict

reflist = {}
ind_fam = {}
gene_families = {}
true_p = 0
false_p = 0
true_n = 0
false_n = 0
unknown_cases = 0
fams = []
samps = []
gene_fams_samp = {}


results_path = sys.argv[1]
reflist_f = results_path + "/SNP_Analysis/Results/pVAAST/LI_ASD/rec/LI_ASD_rec.reflist"
csv_f = results_path + "/SNP_Analysis/Results/pVAAST/LI_ASD/rec/LI_ASD_rec.csv"
var_res = results_path + "/SNP_Analysis/Results/LI_ASD_rec_var.txt"
genes = results_path + "/SNP_Analysis/Results/LI_rec_genes_fam.txt"

with open(reflist_f) as f1:
	for line1 in f1:
		if 'PED_FILE' in line1:
			continue
		else:
			cols = line1.split('\t')
			sampid = cols[0]
			aff_status = cols[5]
			ped_file = cols[1]
			ped_ID = cols[2]
			ids = re.findall(r'\d+', ped_file)
			family = ids[1]
			reflist.setdefault(family, [])
			reflist[family].append([sampid, aff_status, ped_ID])
			ind_fam[sampid] = family	

with open(csv_f) as f:
	with open(var_res,'w') as f2:
		gene_families = defaultdict(list)
		data = csv.reader(f)
		next(f)
		next(f)
		for line in data:
			genotype = line[18]
			if line[1] != '':
				gene = line[1]
				gene_families[gene] = [[],[]]
			aff_samples = [int(s) for s in re.findall(r'\b\d+\d+\b', genotype)]
			#print(gene, aff_samples)
			ranges = re.findall(r'\b\d+[-]\d+\b', genotype)
			for r in ranges: # for each range found in genotype
				bounds = re.findall(r'\b\d+\b', r)	
				difference = int(bounds[1]) - int(bounds[0])
				if difference == 1:
					continue
				else:
					for i in range(difference):
						num = int(bounds[1]) - i
						if num in aff_samples:
							continue
						else:
							aff_samples.append(num)
			ind_genotypes = genotype.split(';')
			#print(ind_genotypes)
			for g in ind_genotypes:
				if '^' in g:
					unknown_cases = unknown_cases + 1
				else:
					sample_ids = re.findall(r'\d+', g)
					proband = sample_ids[0]
					family = ind_fam[proband]
					#print(family, gene, fams)
					gene_families[gene][0].append(family)
					#print(gene, fams)	
					for samples in reflist[family]:
						if proband in samples:
							ped_ID = samples[2]
							ped_ID = ped_ID.strip('.vat.gvf')
					if ped_ID in samps:
						continue
					else:
						gene_families[gene][1].append(ped_ID)
					for ind in reflist[family]:
						ind_id = int(ind[0])
						status = ind[1].strip()
						if ind_id in aff_samples and status == 'Aff':
							true_p = true_p + 1
						#if ind_id in aff_samples and status == 'Unaff':
						#	false_p = false_p + 1
						#if ind_id not in aff_samples and status == 'Aff':
						#	false_n = false_n + 1
						#if ind_id not in aff_samples and status == 'Unaff':
						#	true_n = true_n + 1
			fams = list(set(fams))
			#print(gene, gene_families[gene])			
			#gene_families[gene] = ([fams,samps])
			#print('clearing families')
			#fams = []
			#samps = []
			# write to file
			for elem in line:
				f2.write(elem + '\t')
			f2.write(str(true_p) + '\n')
			true_p = 0
			#true_n = 0
			#false_p = 0
			#false_n = 0				


for gene in gene_families.keys():
	gene_families[gene][0] = list(set(gene_families[gene][0]))
	gene_families[gene][1] = list(set(gene_families[gene][1]))
					
with open(genes,'w') as f3:
	for gene in gene_families:
		f3.write(gene + '\t')
		fams = gene_families[gene][0]
		samples = gene_families[gene][1]
		for fam in fams:
			if fams.index(fam) == len(fams) - 1:
				f3.write(fam)
			else:
				f3.write(fam + ',')
		f3.write('\t')
		for samp in samples:
			if samples.index(samp) == len(samples) - 1:
				f3.write(samp)
			else:
				f3.write(samp + ',')
		f3.write('\n')
