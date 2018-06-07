import io

vcfFile = open('../../data/merge.vcf', 'r')
newFile = open('../../data/new_vcf.vcf', 'w')
vcf_buffer = io.BufferedReader(vcfFile)

for line in vcfFile:

	# Goals:
	# * Delete all FORMATs except GT
	# * Delete INFO except SF, AC and AN (0,2)

	#line = line.decode('utf-8')
	line = line.replace('\n', '')
	if line[0] == '#' and line[1] != '#':
		columns_labels = line.split('	')			
		variants_dict = dict.fromkeys(columns_labels)
	elif line[1] != '#':
		variants_info = line.split('	')

		if str(variants_info[7].split(';')[2][0:1]) == 'AN':
			variants_info[7] = str(variants_info[7].split(';')[0]) + ';' + str(variants_info[7].split(';')[2]) 
		else:
			variants_info[7] = str(variants_info[7].split(';')[0])

		variants_info[8] = str(variants_info[8].split(':')[0])

		for i, gt in enumerate(variants_info[9:]):
			index = i + 9
			variants_info[index] = gt.split(':')[0]

		line = ''
		for i, element in enumerate(variants_info):
			
			if i != 0:
				line = str(line) + '	' + str(element)
			else:
				line = str(element)

	newFile.write(line + '\n')

vcfFile.close()
newFile.close()