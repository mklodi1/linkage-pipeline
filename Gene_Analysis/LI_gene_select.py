
all_genes = []
t1 = []
t2 = []
t3 = []
t4 = []
v = []
b = []

with open('/lab01/Projects/Mudassir_Projects/linkage/Variant_Analysis/Gene_Analysis/Results/LI_gene_table_non_benign_new.txt') as f:
	for line in f:
		if 'Gene' in line:
			continue
		else:
			cols = line.split('\t')
			gene = cols[0]
			known_ndd = cols[1]
			snv = cols[2]
			cnv = cols[3]
			sv = cols[4]
			variant_types = [snv,cnv,sv]
			variants = variant_types.count('1')
			gtex = cols[5]
			bs = cols[6]
			hdbr = cols[7]
			brain = [gtex,bs,hdbr]
			brain_count = brain.count('1')
			unique_fam = int(cols[31])
			
			if variant_types[0] == '1':
				v.append('SNV')
			if variant_types[1] == '1':
				v.append('CNV')
			if variant_types[2] == '1':
				v.append('SV')				
			
			if brain[0] == '1':
				b.append('GTEx')
			if brain[1] == '1':
				b.append('BrainSpan')
			if brain[2] == '1':
				b.append('HDBR')			

	
			var = (',').join(v)
			br = (',').join(b)

			


			if variants >= 2 and brain_count >= 2 and unique_fam >= 2:
				t1.append([gene, 'Tier 1', known_ndd, str(variants), var, str(brain_count), br, str(unique_fam)])
			
			elif variants >= 1 and brain_count >= 1 and unique_fam >= 2:
				t2.append([gene, 'Tier 2', known_ndd, str(variants), var, str(brain_count), br, str(unique_fam)])
			else:
				t3.append([gene, 'Tier 3', known_ndd, str(variants), var, str(brain_count), br, str(unique_fam)])
			all_genes.append(gene)

		v = []
		var = ''
		b = []
		br = ''

with open('/lab01/Projects/Mudassir_Projects/linkage/Variant_Analysis/Gene_Analysis/Results/LI_gene_summary.txt','w') as f2:	
	f2.write('Total Genes: ' + str(len(all_genes)) + '\n')
	f2.write('Tier 1: ' + str(len(t1)) + ' genes' + '\n')
	f2.write('Tier 2: ' + str(len(t2)) + ' genes' + '\n')
	f2.write('Tier 3: ' + str(len(t3)) + ' genes' + '\n')

with open('/lab01/Projects/Mudassir_Projects/linkage/Variant_Analysis/Gene_Analysis/Results/LI_gene_tiers.txt','w') as f1:
	f1.write('Gene' + '\t' + 'Tier' + '\t' + 'Known_NDD' + '\t' + 'var_type_count' + '\t' + 'var_types' + '\t' + 'brain_exp_count' + '\t' + 'brain_exp_db' + '\t' + 'unique_fam_count' + '\n')
	for gene in t1:
		f1.write('\t'.join(gene))
		f1.write('\n')
	for gene in t2:
		f1.write('\t'.join(gene))
		f1.write('\n')
	for gene in t3:
		f1.write('\t'.join(gene))
		f1.write('\n')
