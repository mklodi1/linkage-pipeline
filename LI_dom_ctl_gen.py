"""
LI_dom_ctl_gen.py

> generates CTL file for LI dominant run

"""


import os

dir_path = r'/lab01/Projects/Mudassir_Projects/linkage/Pedigree/dom/LI_ASD'
count = 0
with open('/lab01/Projects/Mudassir_Projects/linkage/CTL/LI_ASD_dom.ctl', 'w') as f:
	f.write('#-------------------------   Basic options----------------------------\n')
	f.write('input_ped_cdr_files:    ')
	for file_path in os.listdir(dir_path):
			fam = file_path[4:8]
			fam_number = int(fam)
			cdr_file = '/nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/CDR/fam_%s.filtered.cdr' % (fam)
			if os.path.exists(cdr_file) == True:
				f.write('/lab01/Projects/Mudassir_Projects/linkage/Pedigree/dom/LI_ASD/' + file_path + '\t')
				cdr_file = '/nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/CDR/fam_%s.filtered.cdr' % (fam)
				f.write(cdr_file + '\t')
	f.write('\n')
	f.write('pedigree_representatives:\n')
	f.write('unknown_representatives:        yes\n')
	f.write('additional_cases:\n')
	f.write('inheritance_model:      dominant               # dominant / recessive\n')
	f.write('\n')
	f.write('#--------------------   Performance Tuning ----------------------------\n')
	f.write('informative_site_selection:     3\n')
	f.write('simulate_genotyping_error:      yes\n')
	f.write('genotyping_error_rate:  1e-4\n')
	f.write('# penetrance_lower_bound:       0.99\n')
	f.write('penetrance_lower_bound: 0.5\n')
	f.write('penetrance_upper_bound:         1\n')
	f.write('genotype_caching_limit: 1000000\n')
	f.write('\n')
	f.write('#----------------   Gene and Variant Filtering ------------------------\n')
	f.write('\n')
	f.write('# Disease prevalence: maximum rate at which the observed phenotype occurs in the population\n')
	f.write('# Setting this option contrains the search space explored during LOD score calculation max_prevalence_filter:  0.01\n')
	f.write('max_prevalence_filter:  0.01\n')
	f.write('lod_score_filter:       yes\n')
	f.write('clrt_score_filter:      yes\n')
	f.write('nocall_filter:  yes\n')
	f.write('nocall_filter_cutoff:   2\n')
	f.write('inheritance_error_filter:       no # set to no for scoring de novo mutations\n')
	f.write('\n')
	f.write('#-----------------------  Developer Options ---------------------------\n')
	f.write('clrt_randomization_round:       100\n')
	f.write('locus_heterogeneity_penalty:    -1\n')
	f.write('incomplete_penetrance_penalty:  -1\n')
	f.write('mcmc_use_functional_score: yes # use a functional score in the MCMC algorithm')


