"""
RI_SV_analysis.py

> extracts SVs and corresponding genes within reading impairment loci (from linkage analysis)

"""

import csv
from collections import defaultdict
RI_genes = []
inheritance = {}
RI_var_count = {}

import sys
results_path = sys.argv[1]
linkage_tsv = results_path + "/SV_Analysis/Results/RI_SV_linkage.tsv"
SV_genes = results_path + "/SV_Analysis/Results/RI_SV_genes.txt"
SV_inheritance = results_path + "/SV_Analysis/Results/RI_SV_inheritance.txt"

with open('/lab01/Projects/Sammy_Projects/NJLAGS_SVs/data/inh_candidates.csv') as f:
	with open(linkage_tsv,'w') as f1:
		data = csv.reader(f)
		inheritance = defaultdict(list)
		for line in data:
			if 'AnnotSV_ID' in line:
				header = line
				header_len = len(header)
				for elem in header:
					if (header.index(elem)) == header_len-1:
						f1.write(elem + '\n')
					else:
						f1.write(elem + '\t')
			else:
				length = len(line)
				chrom = line[1]
				start = line[2].strip()
				end = line[3].strip()
				start = int(start)
				end = int(end)
				genes = line[287]
				transcript = line[289].strip()
				if line[286] == 'split':
					if chrom == '16':
						if start > 17000000 and end < 27000000:
							for elem in line:
								f1.write(elem + '\t')
							f1.write('\n')
							inh = (line[397],line[400],line[402],line[405],line[407],line[410])
							candidates = genes.split(';')
							for gene in candidates:
								if gene in RI_var_count.keys():
									cur_count = RI_var_count[gene]
									cur_count = cur_count + 1
									RI_var_count[gene] = cur_count
								else:
									RI_var_count[gene] = 1

								if gene == ' ' or transcript == '':
									continue
								else:
									inheritance[gene] = inh
									add = [gene,transcript]
									if add in RI_genes:
										continue
									else:
										RI_genes.append(add)		
for elem in RI_genes:
	if elem[1] == '':
		print(elem[0])
		RI_genes.remove(elem)

#print(LI_genes)
#print(inheritance)

with open(SV_genes,'w') as f2:
	for gene in RI_genes:
		count = RI_var_count[gene[0]]
		count = str(count)
		f2.write(gene[0] + '\t' + gene[1]  + '\t' + count + '\n')

with open(SV_inheritance,'w') as f3:
	for gene in inheritance:
		patterns = inheritance[gene]
		f3.write(gene + '\t')
		for pattern in patterns:
			#pattern = pattern.replace('[','')
			#pattern = pattern.replace(']','')
			#print(pattern)
			i = patterns.index(pattern)
			f3.write(pattern + '\t')
		f3.write('\n')

