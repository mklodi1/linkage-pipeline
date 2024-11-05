import sys
results_path = sys.argv[1]
linkage = results_path + "/CNV_Analysis/Results/RI_CNV_linkage.tsv"
CNV_genes = results_path + "/CNV_Analysis/Results/RI_CNV_genes.txt"

RI_genes = []

with open('/lab01/Projects/Rohan_Projects/CNV_Project/2023/CNV_Pipeline/Results/ASD_RI_results.tsv') as f:
	with open(linkage,'w') as f2:
		for line in f:
			if 'start' in line:
				f2.write(line)
			else:
				cols = line.split('\t')
				chrom = cols[0]
				start = cols[1].strip()
				end = cols[2].strip()
				start = int(start)
				end = int(end)
				genes = cols[6]
				if chrom == '16':
					if start > 17000000 and end < 27000000:
						f2.write(line)
						candidates = genes.split('|')
						print(candidates)
						for gene in candidates:
							if gene not in RI_genes:
								RI_genes.append([gene,'2'])

				

with open(CNV_genes,'w') as f2:
	for gene in RI_genes:
		f2.write(gene[0] + '\t' + gene[1] + '\n')
