# Linkage Analysis
This repository contains the code and description for the paper "Narrow Versus Broad Phenotype Definitions Affect Genetic Analysis of Language More Than Other Broad Autism Phenotype Traits".

## Workflow

![alt text][logo]

[logo]: https://github.com/mklodi1/linkage-pipeline/blob/main/Linkage_Analysis_Workflow.png "Analysis workflow"

## Pipeline Structure
The workflow.sh script will run the entire pipeline from start to finish, with the five main components listed below run sequentially. Each of the five main components within this pipeline perform various tasks and analyses contributing to the final results, and will also generate individual results folders pertaining to each specific analysis. The purpose of each file is described below in accordance to its component. 

### Code (Pedigree and CTL File Generation)
#### Objective
Generates pedigree files following dominant and recessive inheritance patterns for each given family. After this, CTL files are generated for each of the four analyses (LI_dom, LI_rec, RI_dom, RI_rec) using the generated corresponding pedigree files. 
#### Files
1. LI_ped_dom.py - generates family pedigree files focusing on LI* patients following dominant inheritance model

2. LI_ped_rec.py - generates family pedigree files focusing on LI* patients following recessive inheritance model

3. RI_ped_dom.py - generates family pedigree files focusing on RI* patients following dominant inheritance model

4. RI_ped_rec.py - generates family pedigree files focusing on RI* patients following recessive inheritance model

5. LI_dom_ctl_gen.py - generates CTL file for LI* dominant run

6. LI_rec_ctl_gen.py - generates CTL file for LI* recessive run

7. RI_dom_ctl_gen.py - generates CTL file for RI* dominant run

8. RI_rec_ctl_gen.py - generates CTL file for RI* recessive run

### SNP_Analysis
#### Objective
Performs analysis on SNP variants identified by pVAAST and generates SNP-related candidate gene list. 
#### Files
1. LI_SNP_analysis.py - extracts LI* SNP candidate genes identified by pVAAST

2. LI_dom_seg_parse.py - under dominant inheritance pattern, identifies samples and families containing each language impairment SNP candidate risk gene

3. LI_rec_seg_parse.py - under recessive inheritance pattern, identifies samples and families containing each language impairment SNP candidate risk gene

4. RI_SNP_analysis.py - extracts RI* SNP candidate genes identified by pVAAST

5. RI_dom_seg_parse.py - under dominant inheritance pattern, identifies samples and families containing each reading impairment SNP candidate risk gene

6. RI_rec_seg_parse.py - under recessive inheritance pattern, identifies samples and families containing each reading impairment SNP candidate risk gene

### CNV_Analysis
#### Objective
Performs analysis on CNV variants obtained from a previous study (https://github.com/JXing-Lab/NJLAGS_SV) and compiles CNV-related candidate gene list. 
#### Files
1. LI_CNV_analysis.py - extracts CNVs and corresponding genes within language impairment loci (from linkage analysis)

2. LI_CNV_fam.py - compiles list of affected LI* individuals and corresponding families containing mutation within each CNV candidate gene

3. RI_CNV_analysis.py - extracts CNVs and corresponding genes within reading impairment loci (from linkage analysis)

4. RI_CNV_fam.py - compiles list of affected RI* individuals and corresponding families containing mutation within each CNV candidate gene

### SV_Analysis
#### Objective
Performs analysis on SV variants obtained from a previous study (https://github.com/JXing-Lab/NJLAGS_SV) and compiles SV-related candidate gene list.
#### Files
1. LI_SV_analysis.py - extracts SVs and corresponding genes within language impairment loci (from linkage analysis)
   
2. LI_SV_benign_genes.py - compiles final list of genes and corresponding families/samples for LI*
   
3. LI_benign_filt.py - filters out benign SVs
   
4. RI_SV_analysis.py - extracts SVs and corresponding genes within reading impairment loci (from linkage analysis)
   
5. RI_SV_benign_genes.py - compiles final list of genes and corresponding families/samples for LI*

6. RI_benign_filt.py - filters out benign SVs

### Gene_Analysis
#### Objective
Compiles results from SNV analysis, CNV analysis, and SV analysis, and performs filtration based on criteria presented in Workflow to compile final candidate gene lists for LI* and RI*. 
#### Files
1. LI_gene_table.py - generates final table containing LI* genes from SNV, CNV, and SV analyses

2. LI_gene_select.py - automated selection of LI* genes based on above described criteria 

3. RI_gene_table.py - generates final table containing RI* genes from SNV, CNV, and SV analyses

4. RI_gene_select.py - automated selection of RI* genes based on above described criteria
