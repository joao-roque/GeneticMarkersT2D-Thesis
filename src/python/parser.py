import vcf
import numpy as np 
import pandas as pd 
import argparse
import csv
import sys
import gzip as gz 
import os
import collections
import io
import psutil
import multiprocessing as mp
from numba import cuda

# VCF Parser - parser.py  <input>  <output> 
# Storing raw control data locally 

def main():

	# Argument parser
	arg_parser = argparse.ArgumentParser(description = 'VCFv4.1 parser that is focused on gzip usage to reduce storage space. \
		It can convert a vcf to csv, clean it to be organized by samples x features, and from a group of control vcfs gather the same variants\
		present in the case file (aswell as cleaning it up too). It also appends a cases to a controls file')
	arg_parser.add_argument('-i', help = 'Point to input VCF file or CSV (gzip)', type = str)
	arg_parser.add_argument('-t', help = 'Clean for ML use (y) - This means organized by samples x features', type = str)	
	arg_parser.add_argument('-o', help = 'Name for new output csv file', type = str)
	arg_parser.add_argument('-c', help = 'Point to folder of control samples', type = str)
	arg_parser.add_argument('-f', help = 'Point to folder of control samples with correct variants to add clean and generate single file')
	arg_parser.add_argument('-s', help = 'Name for new controls file output', type = str)
	arg_parser.add_argument('-cases', help = 'Cases file to append to controls', type = str)
	arg_parser.add_argument('-controls', help = 'Controls file to append to cases', type = str)
	arg_parser.add_argument('-save', help = 'Save appended dataset')
	arg_parser.add_argument('-missing', help = 'Delete columns with more than 10%% missing data (y) to -i file')
	# input in controls and cases uncleaned files (raw csvs from vcf)
	arg_parser.add_argument('-validate', help = 'Validate if controls and cases have the same genotypes for all variants (y)')
	args = arg_parser.parse_args()

	# loads up file
	if args.i:
		print('Loading file...')
		variants_df  = pd.DataFrame()
		if args.i.split('.')[-2] == 'vcf':
			variants_df, vcf_header = vcf_to_csv(args.i)
			print('\nFile successfully loaded')

			if args.o:
				variants_df.to_csv(args.o, compression = 'gzip')

		elif args.i.split('.')[-2] == 'csv':
			variants_df = pd.read_csv(args.i, compression = 'gzip')
			print(variants_df.head())
			print('\nFile successfully loaded')

		# Saves clean variants or raw vcf data both in csv with gzip compression
		if args.t == 'y' or args.t == 'Y':
			print('\nCleaning...')		
			clean_variants = clean_it(args, variants_df)
			

		# returns controls with same variants as case
		if args.c:
			print('Started control parsing...')
			control_variants = parse_controls(args, variants_df)
			'''
			if args.s:
				control_variants.to_csv(args.s + '.gz', compression = 'gzip')
				print('csv saved with gzip compression at: ' + args.s)
			else:
				control_variants.to_csv('controls.gz', compression = 'gzip')
				print('csv saved with gzip compression at: ' + 'controls.gz')
		'''

		if args.missing:
			cleaning_columns(variants_df, args)


	else:
		print('Input file not loaded. Following operations will not be performed: -t, -c, -missing')

	if args.f:
		parse_final_controls(args)

	if args.cases and args.controls:

		if args.validate:
			validate(args)
		else:
			new_dataset = append_dataset(args)
			print('Saving...')
			if args.save:			
				new_dataset.to_csv(args.save + '.csv.gz', compression='gzip')
				print('csv saved with gzip compression at: ' + args.save)
			else:
				new_dataset.to_csv('dataset' + '.gz', compression='gzip')
				print('csv saved with gzip compression at current folder')
	else:
		print('This operation requires both cases and controls')


# Creates a csv and returns dataframe of vcf file
# Also returns the header
def vcf_to_csv(variants_file):

	#vcfFile = open(args.i, 'r')
	#print('Loading...')
	controls_df = pd.DataFrame()
	try:
		vcfFile = None
		vcfFile = gz.open(variants_file, 'rb')
	except Exception as e:
		print(variants_file + ' didnt open properly\n')
		print(e)

	# buffer can be up to 3 times faster
	if vcfFile:
		counter = 0
		variants_buffer = io.BufferedReader(vcfFile)

		for line in variants_buffer:
			counter += 1
			line = line.decode('utf-8')
			if line[0] == '#' and line[1] != '#':
				columns_labels = line.split('	')			
				variants_dict = dict.fromkeys(columns_labels)
			elif line[1] != '#':
				variants_info = line.split('	')

				# Add new line to dict
				for i in range(0, len(columns_labels)):
					if variants_dict[columns_labels[i]] == None:
						variants_dict[columns_labels[i]] = [variants_info[i]]
					else:	  
						variants_dict[columns_labels[i]].append(variants_info[i]) 

			# Check memory usage
			process = psutil.Process(os.getpid())
			sys.stdout.write('\r')
			sys.stdout.write('Memory usage:' + str(process.memory_info().rss/1000000) + 'MB')
			sys.stdout.flush()
		
		variants_df = pd.DataFrame.from_dict(variants_dict)
		print(variants_df.head())
		return((variants_df, columns_labels))

		'''
		if args[1]:
			variants_df.to_csv(args.o, compression='gzip')
		else:
			variants_df.to_csv('csv_from_vcf.csv', compression='gzip')
		'''
	else:
		return((None, None))

# Cleans csv so that it can be used for ML
# Example:
#
#    -   chr1:1010  ,  chr1: 1232, ...
# exome1    1				3
# exome2	1				2
#  ...
def clean_it(args, variants_df):
	
	#print(variants_df.head())

	#variants_df = variants_df.replace('|', '/')

	# Get exome indexes
	exome_index_start = False
	for index, name in enumerate(list(variants_df)):
		if 'exome' in name:
			exome_index_start = index
			break

	# Control file samples always start at 4
	if not exome_index_start:
		exome_index_start = 5

	clean_df = pd.DataFrame()
	# Possible Genotypes
	possible_genotypes = ['0/0', '0/1', '1/1', '0/2', '1/2', '2/2',\
	 '0/3', '1/3', '2/3', '3/3', '0/4', '1/4', '2/4', '3/4', '4/4',\
	 '0/5', '1/5', '2/5', '3/5', '4/5', '5/5', '0/6', '1/6', '2/6',\
	 '3/6', '4/6', '5/6', '6/6']

	# Make new column for each variant
	maxLen = len(variants_df)
	for index, row in variants_df.iterrows():
		
		new_variant = list()
		alt_number = len(row['ALT'].split(','))
		for item_index in range(exome_index_start, variants_df.shape[1]):			
			
			max_variants = sum(list(range(1, alt_number + 1)))			
		
			# make column to add
			# it shouldn't be bigger than 5

			# Adjust so smaller number is always first
			# 2/1 --> 1/2
			# Since genotypes are unphased there is no difference
			# Empty genotypes are set to None

			# Parse genotypes in this format:
			# GT:...:...
			if ':' in str(row[item_index]):
				row[item_index] = row[item_index].split(':')[0]

			if str(row[item_index]) == 'nan' or row[item_index] == None\
			or str(row[item_index]) == '.':
				genotype = None
			elif len(str(row[item_index])) <= 2:
				genotype = str(row[item_index])
				genotype = genotype.strip('\n')
				if genotype == '.':
					genotype = None
			elif len(str(row[item_index])) <= 6:
				genotype = str(row[item_index])
				genotype = genotype.strip("'")
				genotype = genotype.strip('\n')
				genotype = genotype.replace('|', '/')

				try:
					if int(genotype.split('/')[0]) > int(genotype.split('/')[1]):
						genotype = genotype.split('/')[1] + '/' + genotype.split('/')[0] 

					genotype = int(possible_genotypes.index(genotype))
				except:
					genotype = None
					print(' 1. Incorrect genotype at: ' + str(row['#CHROM']) + ':' + str(row['POS']) + ' ' + str(row[item_index]))
			else:
				print(' 2. Incorrect genotype at: ' + str(row['#CHROM']) + ':' + str(row['POS']) + ' ' + str(row[item_index]))
				genotype = None

			new_variant.append(genotype)
		
		if 'chr' not in str(row['#CHROM']):
			clean_df['chr' + str(row['#CHROM']) + ':' + str(row['POS'])] = new_variant
		else:
			clean_df[str(row['#CHROM']) + ':' + str(row['POS'])] = new_variant

		# Progress "bar"
		#process = psutil.Process(os.getpid())
		sys.stdout.write('\r')
		to_finish = 100 * index/maxLen
		sys.stdout.write("%.2f" % to_finish + '%')
		#sys.stdout.write('Memory usage:' + str(process.memory_info().rss/1000000) + 'MB')
		sys.stdout.flush()

	print(clean_df.head())

	if args.o:
		if args.o == 'n' or args.o == 'N':
			print('File wasnt saved')
		else:
			clean_df.to_csv(args.o + '.gz', compression='gzip')
			print('csv saved with gzip compression at: ' + args.o)
	else:
		clean_df.to_csv('CleanVariants.csv', compression='gzip')
		print('csv saved with gzip compression at: ' + 'CleanVariants.gz')

	return(clean_df)

# Compares control files to cases vcf (or csv) and gathers similar variants
# Also outputs details.txt with info of ordered samples and files not loaded
def parse_controls(args, variants_df):

	# Total file with relevant control variants
	control_variants = pd.DataFrame()

	# load each 6, and then save and close files
	already_processed = []
	filesList = os.listdir(args.c)

	try:
		# multiprocessing (with carefull read to avoid memory shortage)
		# mp.cpu_count() . Only using 6 cores because of memory
		pool = mp.Pool(processes = 6)
		results = [pool.apply_async(read_and_compare_controls, (filename, variants_df, args)) for filename in filesList]
		output = [p.get() for p in results]
		pool.close()

	except Exception as e:
		print(e)
		print('Error at: ' + filename)
		pool.close() 

	print('Finished parsing controls.\n')
	#control_variants = pd.concat(output, axis = 0)

	print('Cleaning...')
	parse_final_controls(args)


def parse_final_controls(args):

	filesList = os.listdir(args.f)
	files = []

	print('Loading files...')
	for filename in filesList:
		try:	
			vcfFile = None
			vcfFile_df = pd.read_csv(args.f + '/' + filename, compression = 'gzip')
		except Exception as e:
			print(filename + ' didnt open properly\n')
			with open('details.txt') as details_file:
				details_file(filename + " failed to open because: " + e + "\n")
			print(e)

		files.append(vcfFile_df)

	print('Concatenating files...')
	control_variants = pd.concat(files, axis = 0)
	control_variants.to_csv('check_controls.csv.gz', encoding = 'utf8', compression = 'gzip')

	print('Cleaning controls file...')
	control_variants = clean_it(args, control_variants)

	control_variants.to_csv(args.s + '/parsed_controls.csv', compression = 'gzip')
	print('Parsed controls file saved at: ' + args.f + '/parsed_controls.csv')

# This function will only store relevant lines to avoid memory issues

def read_and_compare_controls(filename, variants_df, args):

	print('\nLoading: ' + filename)
	controls_df = None
	# set is faster than list
	with gz.open(args.c + '/' + filename, 'rb') as f:	
		for lines_no, l in enumerate(f):
			pass
		lines_no = lines_no + 1

	variants_set = set(list(variants_df))

	try:	
		vcfFile = None
		vcfFile = gz.open(args.c + '/' + filename, 'rb')
		variants_buffer = io.BufferedReader(vcfFile)
		print('\nLoaded')
		print('Processing...')
	except Exception as e:
		print(filename + ' didnt open properly\n')

		with open('details.txt') as details_file:
			details_file(filename + " failed to open because: " + e + "\n")
		print(e)

	try:
		if vcfFile:

			# first column will be of samples
			controls_df = pd.DataFrame()
			current_line = 0

			for line in variants_buffer:
				
				line = line.decode('utf-8')
				if line[0] == '#' and line[1] != '#':
					columns_labels = line.split('	')			
					
					#variants_dict = dict.fromkeys(columns_labels)				
					controls_df = pd.DataFrame(columns = columns_labels)

					for index, label in enumerate(columns_labels):
						if label == 'POS':
							pos = index
						elif label == '#CHROM':
							chrom = index						
				elif line[1] != '#':
					line = line.split('	')

					# see if this variant is present in variants_df
					if 'chr' + line[chrom] + ':' + line[pos] in variants_set:
						# Control samples start in position 9
						# controls_df['chr' + line[chrom] + ':' + line[pos]] = line[9:len(line)-1]
						controls_df.loc[len(controls_df)] = line

				if current_line%50000 == 0:
					print('\n' + filename + ' current progress: ' + str(current_line/lines_no))

				current_line += 1

			print(filename)
			print(controls_df.head())

		filename = filename + '.csv'

		if args.f:
			controls_df.to_csv(args.f + '/' + filename , compression='gzip')
		else:
			controls_df.to_csv('../../../LocalData/Parsed_Controls/' + filename, compression = 'gzip')
		print('Saved ' + filename)

	except Exception as e:
		return(None)
		with open('details.txt') as details_file:
			details_file(filename + " failed to open because: " + e + "\n")
		print(e)

	return(controls_df)		

def append_dataset(args):

	print('Loading cases...')
	cases = pd.read_csv(args.cases, compression = 'gzip')
	
	print('Loading controls...')
	controls = pd.read_csv(args.controls, compression = 'gzip')

	print('Cases size: ' + str(len(cases.columns)))
	print('Controls size: ' + str(len(controls.columns)))

	# Adding labels
	cases['labels'] = np.ones(shape = (len(cases.index), 1), dtype = int)
	controls['labels'] = np.zeros(shape = (len(controls.index), 1), dtype = int)

	full_dataset = cases.append(controls)
	print('Full Dataset Metrics:\n')
	print(full_dataset.head(n = 5))
	print(len(full_dataset.columns))

	return(full_dataset)

def cleaning_columns(dataset, args):

	#dataset = pd.read_csv(args.missing, compression = 'gzip')
	dataset_len = len(list(dataset))

	number_of_values = dataset.shape[0] * dataset.shape[1]
	number_of_nan = dataset.isnull().sum().sum()
	print('Missing data percentage: ' + str(number_of_nan / number_of_values))
	print('Dataset size: ' + str(dataset_len))

	# Removing variants with more than 10% missing data
	print('Cleaning...')
	count = 0
	for column in dataset:

		if dataset[column].isnull().sum() / len(dataset[column]) >= 0.1:
			dataset = dataset.drop(column, axis = 1)

		count += 1
		sys.stdout.write('\r')
		to_finish = 100 * count/dataset_len
		sys.stdout.write("Progress: %.2f" % to_finish + '%')
		#sys.stdout.write('Memory usage:' + str(process.memory_info().rss/1000000) + 'MB')
		sys.stdout.flush()

	number_of_values = dataset.shape[0] * dataset.shape[1]
	number_of_nan = dataset.isnull().sum().sum()
	print('Missing data percentage after cleaning: ' + str(number_of_nan / number_of_values))
	print('New dataset size: ' + str(len(list(dataset))))

	print('Saving...')
	dataset.to_csv('../../data/full_dataset/full_clean_dataset.csv.gz', compression = 'gzip')

def validate(args):

	print('Loading cases...')
	cases = pd.read_csv(args.cases, compression = 'gzip', low_memory = False)
	cases = cases.drop(cases.columns[0], axis = 1)

	print('Loading controls...')
	controls = pd.read_csv(args.controls, compression = 'gzip', low_memory = False)
	controls = controls.drop(controls.columns[0], axis = 1)

	print('Cases size: ' + str(len(cases.columns)))
	print('Controls size: ' + str(len(controls.columns)))

	new_cases = pd.DataFrame(columns = list(cases))
	# Adding to new DF
	# Only adding cases because when appending it will resolve itself 

	# dictionary with controls chrm and positions
	print('Building dictionary...')
	controls_dict = {}

	for index_controls, row_controls in controls.iterrows():		
		
		if row_controls['#CHROM'] in controls_dict:
			controls_dict[row_controls['#CHROM']].add(row_controls['POS'])
		else:
			controls_dict[row_controls['#CHROM']] = set()
			controls_dict[row_controls['#CHROM']].add(row_controls['POS'])

		sys.stdout.write('\r')
		to_finish = 100 * (index_controls/len(controls.index))
		sys.stdout.write("Progress: %.2f" % to_finish + '%')
		#sys.stdout.write('Memory usage:' + str(process.memory_info().rss/1000000) + 'MB')
		sys.stdout.flush()


	print('\nParsing...')
	row_counter = 0
	for index_cases, row_cases in cases.iterrows():

		# sometimes it's in the format chr*
		row_cases['#CHROM'] = row_cases['#CHROM'].replace('chr', '')

		# verify if row is in controls
		if row_cases['POS'] in controls_dict[row_cases['#CHROM']]:

			# find occurrence of case in controls and verify if ALT is the same
			try:				

				controls_alt = controls.loc[(controls['#CHROM'] == row_cases['#CHROM']) & (controls['POS'] == row_cases['POS'])]
				if str(controls_alt['ALT'].item()) == str(row_cases['ALT']):
					new_cases.loc[row_counter] = row_cases
					row_counter += 1

			except Exception as e:
				print(str(e) + '\n')
				print('Error comparing ALTs in chrom: ' + str(row_cases['#CHROM']) + ' and pos: ' + str(row_cases['POS']) )

		sys.stdout.write('\r')
		to_finish = 100 * (index_cases/len(cases.index))
		sys.stdout.write("Progress: %.2f" % to_finish + '%')
		#sys.stdout.write('Memory usage:' + str(process.memory_info().rss/1000000) + 'MB')
		sys.stdout.flush()		

	new_cases.to_csv('../../data/cases/new_cases.csv.gz', compression = 'gzip')


if __name__ == "__main__":
	main()

