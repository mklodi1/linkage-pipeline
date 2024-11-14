# Linkage Analysis
This repository contains the code and description for the paper "Narrow Versus Broad Phenotype Definitions Affect Genetic Analysis of Language More Than Other Broad Autism Phenotype Traits".

## Workflow

![alt text][logo]

[logo]: https://github.com/mklodi1/linkage-pipeline/blob/main/Linkage_Analysis_Workflow.png "Analysis workflow"

## Pipeline Structure
The workflow.sh script will run the entire pipeline from start to finish, with five main components. Each of the five main components within this pipeline perform various tasks and analyses contributing to the final results. The purpose of each file is described below in accordance to its component. 

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
1. LI_SNP_analysis.py - generates family pedigree files focusing on LI* patients following dominant inheritance model

2. LI_dom_seg_parse.py - generates family pedigree files focusing on LI* patients following recessive inheritance model

3. LI_rec_seg_parse.py - generates family pedigree files focusing on RI* patients following dominant inheritance model

4. RI_SNP_analysis.py - generates family pedigree files focusing on RI* patients following recessive inheritance model

5. RI_dom_seg_parse.py - 

6. RI_rec_seg_parse.py -

### CNV_Analysis
#### Objective
Performs analysis on CNV variants obtained from a previous study (https://github.com/JXing-Lab/NJLAGS_SV) and compiles CNV-related candidate gene list. 
#### Files
1. LI_CNV_analysis.py -

2. LI_CNV_fam.py -

3. RI_CNV_analysis.py -

4. RI_CNV_fam.py - 
5. 
