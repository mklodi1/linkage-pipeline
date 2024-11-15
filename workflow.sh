#!/bin/sh

# Statistics::Distributions Module
export PATH=/usr/bin:$PATH

dir_path=/lab01/Projects/Mudassir_Projects/linkage
var_results=Variant_Analysis_$(date +"%m%d%Y")
full_path="$dir_path"/"$var_results"

echo "Results are located in "$var_results""
 
mkdir Variant_Analysis_$(date +"%m%d%Y")

# Pedigree File Generation (LI_dom, LI_rec, RI_dom, RI_rec)
python3 "$dir_path"/Code/LI_ped_dom.py
python3 "$dir_path"/Code/LI_ped_rec.py
python3 "$dir_path"/Code/RI_ped_dom.py
python3 "$dir_path"/Code/RI_ped_rec.py

echo "Exception Handling"
python3 "$dir_path"/Code/exception.py

echo "Finished pedigree file generation"

rm "$dir_path"/Pedigree/dom/RI_ASD/fam_1009_RI_ASD_dom.ped
rm "$dir_path"/Pedigree/dom/LI_ASD/fam_1009_LI_ASD_dom.ped 
rm "$dir_path"/Pedigree/dom/LI_ASD/fam_2008_LI_ASD_dom.ped  

# CTL File Generation
python "$dir_path"/Code/LI_dom_ctl_gen.py
python "$dir_path"/Code/LI_rec_ctl_gen.py
python "$dir_path"/Code/RI_dom_ctl_gen.py
python "$dir_path"/Code/RI_rec_ctl_gen.py

echo "Finished CTL file generation"

# Generate Results folder for SNP_Analysis
mkdir "$var_results"/SNP_Analysis/
mkdir "$var_results"/SNP_Analysis/Results/
mkdir "$var_results"/SNP_Analysis/Results/pVAAST/
mkdir "$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/
mkdir "$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/dom/
mkdir "$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/rec/
mkdir "$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/
mkdir "$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/dom/
mkdir "$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/rec/

# Execute pVAAST Command
echo "Running LI_ASD_dom pVAAST analysis"

/lab01/Tools/VAAST_3.0/bin/vaast_tools/../VAAST -m pvaast -o "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/dom/LI_ASD_dom -pv_control /lab01/Projects/Mudassir_Projects/linkage/CTL/LI_ASD_dom.ctl -region /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/LI.bed   /lab01/Projects/VAAST_Projects/Data/Features/refGene_hg19.gff3 /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/bg_gtex.cdr -gw 1e6 -min_suc 50 -p 22 -o /lab01/Projects/Mudassir_Projects/linkage/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/dom/LI_ASD_dom -max_time 300 &> "$var_results"_LI_ASD_dom.txt &

echo "Running LI_ASD_rec pVAAST analysis"

/lab01/Tools/VAAST_3.0/bin/vaast_tools/../VAAST -m pvaast -o "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/rec/LI_ASD_rec -pv_control /lab01/Projects/Mudassir_Projects/linkage/CTL/LI_ASD_rec.ctl -region /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/LI.bed   /lab01/Projects/VAAST_Projects/Data/Features/refGene_hg19.gff3 /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/bg_gtex.cdr -gw 1e6 -min_suc 50 -p 22 -o /lab01/Projects/Mudassir_Projects/linkage/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/rec/LI_ASD_rec -max_time 300 &> "$var_results"_LI_ASD_rec.txt &

echo "Running RI_ASD_dom pVAAST analysis"

/lab01/Tools/VAAST_3.0/bin/vaast_tools/../VAAST -m pvaast -o "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/dom/RI_ASD_dom -pv_control /lab01/Projects/Mudassir_Projects/linkage/CTL/RI_ASD_dom.ctl -region /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/RI.bed   /lab01/Projects/VAAST_Projects/Data/Features/refGene_hg19.gff3 /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/bg_gtex.cdr -gw 1e6 -min_suc 50 -p 22 -o /lab01/Projects/Mudassir_Projects/linkage/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/dom/RI_ASD_dom -max_time 300 &> "$var_results"_RI_ASD_dom.txt &

echo "Running RI_ASD_rec pVAAST analysis"

/lab01/Tools/VAAST_3.0/bin/vaast_tools/../VAAST -m pvaast -o "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/rec/RI_ASD_rec -pv_control /lab01/Projects/Mudassir_Projects/linkage/CTL/RI_ASD_rec.ctl -region /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/RI.bed   /lab01/Projects/VAAST_Projects/Data/Features/refGene_hg19.gff3 /nfs/xing_lab/brz_vcf/VAAST_JUN_3_2019/bg_gtex.cdr -gw 1e6 -min_suc 50 -p 22 -o /lab01/Projects/Mudassir_Projects/linkage/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/rec/RI_ASD_rec -max_time 300 &> "$var_results"_RI_ASD_rec.txt &

wait

echo "Finished pVAAST results"

echo "Creating pVAAST csv files"
# Creates csv files based on .vaast results 
python3 /lab01/Projects/VAAST_Projects/CommonCodes/pVaastParser_LOD.py "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/dom/LI_ASD_dom.vaast "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/dom/LI_ASD_dom.csv
python3 /lab01/Projects/VAAST_Projects/CommonCodes/pVaastParser_LOD.py "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/rec/LI_ASD_rec.vaast "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/LI_ASD/rec/LI_ASD_rec.csv
python3 /lab01/Projects/VAAST_Projects/CommonCodes/pVaastParser_LOD.py "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/dom/RI_ASD_dom.vaast "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/dom/RI_ASD_dom.csv
python3 /lab01/Projects/VAAST_Projects/CommonCodes/pVaastParser_LOD.py "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/rec/RI_ASD_rec.vaast "$dir_path"/"$var_results"/SNP_Analysis/Results/pVAAST/RI_ASD/rec/RI_ASD_rec.csv

echo "Finished creating csv files"

# SNP Analysis
python "$dir_path"/SNP_Analysis/LI_SNP_analysis.py "$full_path"
python "$dir_path"/SNP_Analysis/RI_SNP_analysis.py "$full_path"
python "$dir_path"/SNP_Analysis/LI_dom_seg_parse.py "$full_path"
python "$dir_path"/SNP_Analysis/LI_rec_seg_parse.py "$full_path"
python "$dir_path"/SNP_Analysis/RI_dom_seg_parse.py "$full_path"
python "$dir_path"/SNP_Analysis/RI_rec_seg_parse.py "$full_path"

echo "Finished SNP Analysis"

# Generate Results folder for CNV_Analysis
mkdir "$var_results"/CNV_Analysis/
mkdir "$var_results"/CNV_Analysis/Results/

# CNV Analysis
python "$dir_path"/CNV_Analysis/LI_CNV_analysis.py "$full_path"
python "$dir_path"/CNV_Analysis/LI_CNV_fam.py "$full_path"
python "$dir_path"/CNV_Analysis/RI_CNV_analysis.py "$full_path"
python "$dir_path"/CNV_Analysis/RI_CNV_fam.py "$full_path"

echo "Finished CNV Analysis"

# Generate Results folder for SV_Analysis
mkdir "$var_results"/SV_Analysis/
mkdir "$var_results"/SV_Analysis/Results/

# SV Analysis
python "$dir_path"/SV_Analysis/LI_SV_analysis.py "$full_path"
python "$dir_path"/SV_Analysis/RI_SV_analysis.py "$full_path"

# SV Benign Analysis
python "$dir_path"/SV_Analysis/LI_benign_filt.py "$full_path"
python "$dir_path"/SV_Analysis/LI_SV_benign_genes.py "$full_path"
python "$dir_path"/SV_Analysis/RI_benign_filt.py "$full_path"
python "$dir_path"/SV_Analysis/RI_SV_benign_genes.py "$full_path"

echo "Finished SV Analysis"

#mkdir "$var_results"/Gene_Analysis/

# Gene Table
python "$dir_path"/Gene_Analysis/LI_gene_table.py "$full_path"
python "$dir_path"/Gene_Analysis/RI_gene_table.py "$full_path"

echo "Finished Gene Analysis"

# Gene Selection
python "$dir_path"/Gene_Analysis/LI_gene_select.py "$full_path"
python "$dir_path"/Gene_Analysis/RI_gene_select.py "$full_path"


