"""
RI_benign_filt.py

> filters out benign SVs

"""
import sys
results_path = sys.argv[1]
linkage_tsv = results_path + "/SV_Analysis/Results/RI_SV_linkage.tsv"
SV_genes = results_path + "/SV_Analysis/Results/RI_SV_genes.txt"
SV_inheritance = results_path + "/SV_Analysis/Results/RI_SV_inheritance.txt"
linkage_filt = results_path + "/SV_Analysis/Results/RI_SV_linkage_filt.tsv"

identifier = 0

with open(linkage_tsv) as f:
	with open(linkage_filt,'w') as f1:
		for line in f:
			if 'AnnotSV_ID' in line:
				f1.write(line)
			else:
				cols = line.split('\t')
				sv_type = cols[5]
				b_gain = cols[319]
				b_loss = cols[321]
				b_ins = cols[323]
				start = cols[2]
				end = cols[3]
			
				if sv_type == 'DEL':
					if b_loss == '':
						identifier = 0
					else:
						coords = b_loss.split(';')
						for coord in coords:
							print(coord)
							pos_coord = coord.split(':')
							pos_coord = pos_coord[1].split('-')
							start_coord = pos_coord[0]
							end_coord = pos_coord[1]
							print(cols[0], start_coord, end_coord, start, end)
							if start_coord == start and end_coord == end:
								print(sv_type, 'loss_yes')
								identifier = 1	
				if sv_type == 'INS':
					if b_gain == '':
						identifier = 0
					else:
						gain_coords = b_gain.split(';')
						for coord in gain_coords:
							pos_coord = coord.split(':')
							pos_coord = pos_coord[1].split('-')
							start_coord = pos_coord[0]
							end_coord = pos_coord[1]
							if start_coord == start and end_coord == end:
								print(sv_type, 'gain_yes')
								identifier = 1
					if b_ins == '':
						identifier = 0
					else:
						ins_coords = b_ins.split(';')
						for coord in ins_coords:
							pos_coord = coord.split(':')
							pos_coord = pos_coord[1].split('-')
							start_coord = pos_coord[0]
							end_coord = pos_coord[1]
							if start_coord == start and end_coord == end:
								print(sv_type, 'ins_yes')
								identifier = 1
				print(identifier)
				if identifier == 0:
					f1.write(line)
			identifier = 0
