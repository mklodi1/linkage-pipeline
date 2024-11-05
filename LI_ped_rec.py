"""
LI_ped_rec.py

> generates pedigree files for Language Impairment patients under recessive inheritance model

"""

import re
from collections import Counter
from collections import defaultdict

sample_gvcf = dict()
snv_sub = {}
family_data = {}
individuals = {}
scores = {}
parentIDs = []
subpeds = {}
fathers = []
mothers = []
maxKeys = []
top = 0
final_ped = []
count = 0


#with open("pedigree_summary.csv", "r") as f:
with open("/lab01/Projects/Mudassir_Projects/linkage/Data/cohort_summary_edited_2023-06-27.txt", "r") as f1:
    #with open("LI_rec_ped.txt",'w') as f2:
    #    f2.write("Family" + "\t" + "Subject" + "\t" + "Father" + "\t" + "Mother" + "\t" + "Sex" + "\t" + "Phenotype" + "\n")
	family_data = defaultdict(dict)
	subpeds = defaultdict(dict)
	for line1 in f1:
            PED_columns = line1.split("\t")
            if len(PED_columns) > 1 and 'Father' not in line1:
                family = PED_columns[0]
                subject = PED_columns[1]
                SNV_ID = PED_columns[22]
                snv_sub[subject] = SNV_ID

                PED_family = PED_columns[0]
                PED_subject = PED_columns[1]
                PED_father = PED_columns[2]
                fathers.append(PED_father)
                PED_mother = PED_columns[3]
                mothers.append(PED_mother)
                sex = PED_columns[8]
                LI_phenotype = PED_columns[10] # LI
                RI_phenotype = PED_columns[11] # RI
                ASD_phenotype = PED_columns[9] # ASD
                sequenced = PED_columns[20] # sequencing, 1 or 0
                                

                if LI_phenotype == '2' or ASD_phenotype == '2':
                   phenotype = '2'
                elif LI_phenotype == 'x' or ASD_phenotype == 'x':
                   phenotype = '0'
                else:
                   phenotype = '1'

                ind = [PED_family, PED_subject, PED_father, PED_mother, sex, phenotype]

                # create nested dictionary of family and subjects within family
                family_data[family][PED_subject] = [PED_father, PED_mother, sex, phenotype, sequenced]
                fam = [PED_father, PED_mother, PED_subject]
                parents = (PED_father, PED_mother)
                filePath = '/lab01/Projects/Mudassir_Projects/linkage/Pedigree/rec/LI_ASD/fam_%s_LI_ASD_rec.ped' % (family)
                f2 = open(filePath, 'w')
                
                # create nested dictionary of all possible sub-pedigrees within family
                if parents in subpeds[family].keys() and PED_father != '0' and PED_mother != '0': 
                   subpeds[family][(PED_father,PED_mother)].append(PED_subject)
                else:
                   subpeds[family][(PED_father,PED_mother)] = [PED_father, PED_mother, PED_subject]
        
        # choose sub-pedigree
	for family in subpeds:
            filePath = '/lab01/Projects/Mudassir_Projects/linkage/Pedigree/rec/LI_ASD/fam_%s_LI_ASD_rec.ped' % (family)
            f2 = open(filePath, 'w')
#            print(family)
            subped_perms = subpeds[family].values()
            for perm in subped_perms:
                for ind in perm:
                    if ind in family_data[family].keys() and ind != '0' and family_data[family][ind][4] == '1':
                       count = count + 1
                    if ind in family_data[family].keys() and ind != '0' and family_data[family][ind][3] == '2':
                       count = count + 1
                perm_key = tuple(perm) 
                scores[perm_key] = count
                count = 0
            final_ped = list(max(scores, key=scores.get))

            # determines if there are multiple sub-pedigrees with the max score
            top = max(scores.items(), key=lambda x: x[1])
            top_score = top[1]
            for key, value in scores.items():
                if value == top_score:
                   maxKeys.append(key)
#            if len(maxKeys) > 1:
#               print('Max Keys', maxKeys)
            scores.clear()
            maxKeys.clear()
            top = 0
            top_score = 0

#            print('Final PED', final_ped)
            for p in subpeds[family].keys():
                if subpeds[family][p] == final_ped:
                   parents = p
#            print('Parents', parents)
            
            # format and write to file
            for subject in final_ped:
                if subject == '0':
                   continue
                else:
                   subject_gvf = subject
                   if subject not in family_data[family].keys():
                      father_gvf = '0'
                      mother_gvf = '0'
                      if subject in fathers:
                         sex = '1'
                         phenotype = '0'
                      if subject in mothers:
                         sex = '2'
                         phenotype = '0'
                   else:
                      father_gvf = family_data[family][subject][0]
                      mother_gvf = family_data[family][subject][1]
                      sex = family_data[family][subject][2]
                      phenotype = family_data[family][subject][3]
                      if subject in snv_sub.keys() and subject != '' and snv_sub[subject] != '':
                         subject_gvf = str(snv_sub[subject] + '.vat.gvf')
                      if subject in parents:
                         father_gvf = '0'
                         mother_gvf = '0'
                      else:
                        if family_data[family][subject][0] in snv_sub.keys() and father_gvf != '' and snv_sub[father_gvf] != '':
                                father_gvf = str(snv_sub[family_data[family][subject][0]] + '.vat.gvf')
                        else:
                                father_gvf = family_data[family][subject][0]
                        if family_data[family][subject][1] in snv_sub.keys() and mother_gvf != '' and snv_sub[mother_gvf] != '':
                                mother_gvf = str(snv_sub[family_data[family][subject][1]] + '.vat.gvf') 
                        else:
                                mother_gvf = family_data[family][subject][1]
                f2.write(family + '\t' + subject_gvf + '\t' + father_gvf + '\t' + mother_gvf + '\t' + sex  + '\t' + phenotype + '\n')
            final_ped.clear()            
                

	''' 
                # establish parent IDs
                for ind in individuals:
                    parentIDs.append(ind[2])
                    parentIDs.append(ind[3])

                # establish family founders (if their parental IDs are '0')
                for ind in individuals:
                    if ind[2] == '0' and ind[3] == '0':
                        parentID = ind[1]
                        ancestors.append(parentID)
                #print('Ancestors', ancestors)          

                # determines number of possible sub-pedigrees by counting the number of fathers in the file
                for ind in individuals:
                    fathers.append(ind[2])
                fathers = list(set(fathers))

                for f in fathers:
                    if f != '0':
                       permutations = permutations + 1

                # creates new list pertaining to each sub-pedigree by identifying all subjects with a particular father ID
                for i in range(permutations):
                    subpeds.append([])

                print(subpeds)

                #for ped in subpeds:
                #    print(ped)
                           
                # determine which family founders yield the most number of subjects in the pedigree - change to most number of AFFECTED individuals
                for anc in ancestors:
                    count = parentIDs.count(anc)
                    if count >= top:
                        top = count
                        head.append(anc)
                print(head)
                count = 0
                top = 0

                # establish first level of pedigree (includes grandparents and parents)
                for ind in individuals:
                    if ind[2] in head or ind[3] in head:
                                level.append(ind[1])
                level = list(set(level))
                for ind in individuals:
                    if ind[1] in level:
                        level.append(ind[2])
                        level.append(ind[3])
                level = list(set(level))
                print("level1:", level)

                # add individuals to Head if either parent IDs are in Head
                for ind in individuals:
                    if ind[1] in parentIDs:
                        if ind[2] in head or ind[3] in head:
                            head.append(ind[1])
                print('head', head)

                # determine which couple yields the most subjects in the pedigree (father and mother)
                for ind in head:
                    if ind in ancestors:
                        continue
                    else:
                        count = parentIDs.count(ind)
                        if count > top:
                            top = count
                            head.append(ind)
                head = list(set(head))

                # if the above results in more than 2 IDs (as expected), tiebreaker will be how many subjects in the sub-pedigree are affected individuals
                if len(head) > 2:
                    for ind in individuals:
                        if ind[2] in head and ind[5] == '2':
                                agg.append(ind[2])
                        if ind[3] in head and ind[5] == '2':
                                agg.append(ind[3])
                    head = agg
                count = 0
                top = 0
                print("Final", head)

                # generate list of final subjects in sub-pedigree		
                for ind in individuals:
                    if ind[1] in head: # if it is father
                        final_fam.append(ind[1])
                    if ind[2] in head: # if child
                        final_fam.append(ind[1]) # adds child ID
                        final_fam.append(ind[3]) # adds mother ID
                        head.append(ind[3])
                final_fam = list(set(final_fam))
                head = list(set(head))

                if len(final_fam) == 1:
                        for ind in individuals:
                                if ind[1] == final_fam[0]:
                                        final_fam.append(ind[2])
                                        final_fam.append(ind[3])

                print("Final Fam:", final_fam)

                # format and print pedigree
                for ind in individuals:
                    subject_gvf = ind[1]
                    if ind[1] in head:
                        father_gvf = '0'
                        mother_gvf = '0'
                    else:
                        father_gvf = ind[2]
                        mother_gvf = ind[3]
                    if ind[1] in final_fam:
                        if ind[1] != '0' and ind[1] in snv_sub.keys() and snv_sub[ind[1]] != '':
                            subject_gvf = str(snv_sub[ind[1]] + '.vat.gvf')
                        if ind[2] != '0' and father_gvf != '0':
                                if ind[2] in snv_sub.keys() and snv_sub[ind[2]] != '':
                                        father_gvf = str(snv_sub[ind[2]] + '.vat.gvf')
                                else:
                                        father_gvf = '0'
                        if ind[3] != '0' and mother_gvf != '0':
                                if ind[3] in snv_sub.keys() and snv_sub[ind[3]] != '':
                                        mother_gvf = str(snv_sub[ind[3]] + '.vat.gvf')
                                else:
                                        mother_gvf = '0'
                        f2.write(ind[0] + '\t' + subject_gvf + '\t' + father_gvf + '\t' + mother_gvf + '\t' + ind[4] + '\t' + ind[5] + '\n')
		
                head.clear()
                level.clear()
                final_head.clear()
                final_fam.clear()
                ancestors.clear()
                agg.clear()
                count = 0
                top = 0'''
	#subpeds.clear()
	fathers.clear()
	permutations = 0
