
import sys
results_path = sys.argv[1]
linkage_tsv = results_path + "/CNV_Analysis/Results/LI_CNV_linkage.tsv"
CNV_genes = results_path + "/CNV_Analysis/Results/LI_CNV_genes.txt"

LI_genes = []

with open('/lab01/Projects/Rohan_Projects/CNV_Project/2023/CNV_Pipeline/Results/ASD_LI_results.tsv') as f:
	with open(linkage_tsv,'w') as f1:
		for line in f:
			if 'start' in line:
				f1.write(line)
			else:
				cols = line.split('\t')
				chrom = cols[0]
				start = cols[1].strip()
				end = cols[2].strip()
				start = int(start)
				end = int(end)
				genes = cols[6]
				if chrom == '15':
					if start > 70000000 and end < 93000000:
						f1.write(line)
						candidates = genes.split('|')
						print(candidates)
						for gene in candidates:
							LI_genes.append([gene, '1'])

				

with open(CNV_genes,'w') as f2:
	for gene in LI_genes:
		f2.write(gene[0] + '\t' + gene[1] + '\t'  + '\n')
