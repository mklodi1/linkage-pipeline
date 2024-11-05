"""
RI_SNP_analysis.py

> extracts SNP candidate genes identified by pVAAST and stores into RI_SNP_genes.txt for further analysis

"""
import sys
import csv
RI_genes = []

results_path = sys.argv[1]
dom_csv = results_path + "/SNP_Analysis/Results/pVAAST/RI_ASD/dom/RI_ASD_dom.csv"
rec_csv = results_path + "/SNP_Analysis/Results/pVAAST/RI_ASD/rec/RI_ASD_rec.csv"
var_res = results_path + "/SNP_Analysis/Results/RI_ASD_dom_var.txt"
genes = results_path + "/SNP_Analysis/Results/RI_dom_genes_fam.txt"
final_SNP_genes = results_path + "/SNP_Analysis/Results/RI_SNP_genes.txt"


with open(dom_csv) as f:
	data = csv.reader(f)
	next(f)
	for line in data:
		if 'Gene' in line:
			continue
		else:
			gene_name = line[1]
			if gene_name == '':
				continue
			else:
				transcript_id = line[4]
				add = [gene_name,transcript_id]
				if add in RI_genes:
					continue
				else:
					RI_genes.append(add)

with open(rec_csv) as f2:
	data2 = csv.reader(f2)
	next(f2)
	for line in data2:
		if 'Gene' in line:
			continue
		else:
			gene_name = line[1]
			if gene_name == '':
				continue
			else:
				transcript_id = line[4]
				add = [gene_name,transcript_id]
				if add in RI_genes:
					continue
				else:
					RI_genes.append(add)

print(len(RI_genes))

with open(final_SNP_genes,'w') as f1:
	for gene in RI_genes:
		f1.write(gene[0] + '\t' + gene[1] + '\n')
