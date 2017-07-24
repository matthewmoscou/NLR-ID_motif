#! /usr/bin/python

# parse_fimo_data.py

import re
import sets
import string

# Import fimo output from individual motifs
patterns = range(1,41)
pattern_sequence = {}
gene_position_pattern = {}

for i in patterns:
	pattern_sequence[i] = []

fimo_output = open('fimo_out/fimo.txt', 'r')
fimo_output_data = fimo_output.readlines()
truth = False

for line in fimo_output_data:
	sline = string.split(line)

	if truth:
		gene_position_pattern[sline[1]] = {}
	else:
		truth = True

truth = False

TSV = []

for pattern in patterns:
	if int(sline[0]) < 10:
		TSV.append(open('NLR_ID_FIMO_I0' + str(pattern) + '.tsv', 'w')) 
	else:
		TSV.append(open('NLR_ID_FIMO_I' + str(pattern) + '.tsv', 'w')) 

for line in fimo_output_data:
	sline = string.split(line)

	if truth:
		i = int(sline[0])

		if float(sline[7]) < 0.05:
			if int(sline[2]) in gene_position_pattern[sline[1]].keys():
				print "Here be dragons - 2"

			if int(sline[0]) < 10:
				gene_position_pattern[sline[1]][int(sline[2])] = 'I0' + sline[0]
			else:
				gene_position_pattern[sline[1]][int(sline[2])] = 'I' + sline[0]

			pattern_sequence[i].append(sline[1])

			# Export for QKdomain analysis
			if int(sline[0]) < 10:
				TSV[int(sline[0]) - 1].write(sline[1] + '\t' + 'xxx' + '\t' + 'xxx' + '\t' + 'FIMO' + '\t' + 'I0' + sline[0] + '\t' + 'I0' + sline[0] + '\t' + sline[2] + '\t' + sline[3] + '\t' + sline[7] + '\t' + 'T' + '\t' + '12-07-2017' + '\n')
			else:
				TSV[int(sline[0]) - 1].write(sline[1] + '\t' + 'xxx' + '\t' + 'xxx' + '\t' + 'FIMO' + '\t' + 'I' + sline[0] + '\t' + 'I' + sline[0] + '\t' + sline[2] + '\t' + sline[3] + '\t' + sline[7] + '\t' + 'T' + '\t' + '12-07-2017' + '\n')

	else:
		truth = True

patterns.reverse()

for pattern in patterns:
	TSV[pattern - 1].close()

patterns.reverse()

fimo_output.close()

# Export for iTOL
for pattern in patterns:
	if pattern < 10:
		ID = 'I0' + str(pattern)
		output = open('NLR_ID_FIMO_I0' + str(pattern) + '_ITOL.txt', 'w')
	else:
		ID = 'I' + str(pattern)
		output = open('NLR_ID_FIMO_I' +str(pattern) + '_ITOL.txt', 'w')
	
	output.write('DATASET_BINARY' + '\n')
	output.write('SEPARATOR TAB' + '\n')
	output.write('DATASET_LABEL\t' + ID + '\n')
	output.write('COLOR\t#801b22' + '\n')
	output.write('FIELD_SHAPES\t2' + '\n')
	output.write('FIELD_LABELS\tf1' + '\n')
	output.write('DATA' + '\n')

	for gene in pattern_sequence[pattern]:
		output.write(gene + '\t' + '1' + '\n')
	
	output.close()
