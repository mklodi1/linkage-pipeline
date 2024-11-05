inheritance = {}
LI_genes = []
LI_var_count = {}
import re
import csv
from collections import defaultdict

import sys
results_path = sys.argv[1]
linkage_tsv = results_path + "/SV_Analysis/Results/RI_SV_linkage.tsv"
SV_genes = results_path + "/SV_Analysis/Results/RI_SV_genes.txt"
SV_inheritance = results_path + "/SV_Analysis/Results/RI_SV_inheritance.txt"
linkage_filt = results_path + "/SV_Analysis/Results/RI_SV_linkage_filt.tsv"
non_benign = results_path + "/SV_Analysis/Results/RI_SV_inh_non_benign.txt"
SV_nb_genes = results_path + "/SV_Analysis/Results/RI_SV_genes_filt.txt"

with open(linkage_filt) as f:
	inheritance = defaultdict(list)
	for line in f:
		if 'AnnotSV_ID' in line:
			continue
		else:
			line = line.split('\t')
			length = len(line)
			chrom = line[1]
			start = line[2].strip()
			end = line[3].strip()
			start = int(start)
			end = int(end)
			gene = line[287]
			transcript = line[289].strip()

			inh = [line[367],line[370],line[372],line[375],line[377],line[380]]
			
			families = [line[367],line[372],line[377]]
			samples = [line[370],line[375],line[380]]

			#SV_fam[gene] = [dom_fam,dom_samp,rec_fam,rec_samp,denovo_fam,denovo_samp]
			#fams = re.findall(r'\d+',dom_fam) + re.findall(r'\d+',rec_fam) + re.findall(r'\d+',denovo_fam)

			
			if gene in inheritance.keys():
				inheritance[gene][0].append(line[367])
				inheritance[gene][1].append(line[370])
				inheritance[gene][2].append(line[372])
				inheritance[gene][3].append(line[375])
				inheritance[gene][4].append(line[377])
				inheritance[gene][5].append(line[380])
				
				cur_count = LI_var_count[gene]
				cur_count = cur_count + 1
				LI_var_count[gene] = cur_count

			else:
				inheritance[gene] = [[],[],[],[],[],[]]
				inheritance[gene][0].append(line[367])
				inheritance[gene][1].append(line[370])
				inheritance[gene][2].append(line[372])
				inheritance[gene][3].append(line[375])
				inheritance[gene][4].append(line[377])
				inheritance[gene][5].append(line[380])
				
				LI_var_count[gene] = 1

			#print(gene, families)
			#print(gene, inheritance[gene])
			add = [gene,transcript]
			if add in LI_genes:
				continue
			else:
				LI_genes.append(add)
for elem in LI_genes:
        if elem[1] == '':
                print(elem[0])
                LI_genes.remove(elem)
#print('interest:', inheritance['THSD4'])
#print('count:', LI_var_count['THSD4'])	
#print(LI_genes)


final_inh = {}
final_inh = defaultdict(list)

for gene in inheritance.keys():
	final_inh[gene] = [[],[],[],[],[],[]]
	dom_families = inheritance[gene][0]
	rec_families = inheritance[gene][2]
	dn_families = inheritance[gene][4]
	dom_samples = inheritance[gene][1]
	rec_samples = inheritance[gene][3]
	dn_samples = inheritance[gene][5]
	for fam in dom_families:
		if fam == '':
			continue
		else:
			final_inh[gene][0].append(fam)
	for fam in dom_samples:
		if fam == '':
			continue
		else:
			final_inh[gene][1].append(fam)

	for fam in rec_families:
		if fam == '':
			continue
		else:
			final_inh[gene][2].append(fam)

	for fam in rec_samples:
		if fam == '':
			continue
		else:
			final_inh[gene][3].append(fam)

	for fam in dn_families:
		if fam == '':
			continue
		else:
			final_inh[gene][4].append(fam)

	for fam in dn_samples:
		if fam == '':
			continue
		else:
			final_inh[gene][5].append(fam)
	
	#for samp in samps:
	#	for s in samp:
	#		if s == '':
	#			continue
	#		else:
	#			final_inh[gene][1].append(s)
	
#re.findall(r'\d+',final_inh['THSD4'][0])

#print('THSD4',final_inh['THSD4'])			
with open(SV_nb_genes,'w') as f2:
	for gene in LI_genes:
		count = LI_var_count[gene[0]]
		count = str(count)
		f2.write(gene[0] + '\t' + gene[1]  + '\t' + count + '\n')

with open(non_benign,'w') as f3:
	print(final_inh)
	for gene in final_inh:
		dom_families = final_inh[gene][0]
		dom_samples = final_inh[gene][1]
		rec_families = final_inh[gene][2]
		rec_samples = final_inh[gene][3]
		dn_families = final_inh[gene][4]
		dn_samples = final_inh[gene][5]
		f3.write(gene + '\t')
		for fam in dom_families:
			f3.write(fam + ',')
		f3.write('\t')
		for samp in dom_samples:
			f3.write(samp + ',')
		f3.write('\t')
		for fam in rec_families:
			f3.write(fam + ',')
		f3.write('\t')
		for samp in rec_samples:
			f3.write(samp + ',')
		f3.write('\t')
		for fam in dn_families:
			f3.write(fam + ',')
		f3.write('\t')
		for samp in dn_samples:
			f3.write(samp + ',')
		f3.write('\n')
