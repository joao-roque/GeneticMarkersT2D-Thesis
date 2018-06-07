import numpy as np
import pandas as pd 
import argparse
import scipy.stats as stats
from sklearn.preprocessing import Imputer
import statsmodels.api as sm

def main():

	# Argument parser
	arg_parser = argparse.ArgumentParser(description = 'Feature Selection')
	arg_parser.add_argument('-i', help = 'Input file to select features', type = str)
	arg_parser.add_argument('-important', help = 'File with important variants')
	args = arg_parser.parse_args()


	print('Loading...')
	dataset = pd.read_csv(args.i, compression = 'gzip')
	dataset = dataset.drop(dataset.columns[0], axis = 1)

	print('Performing imputation...')
	dataset = imputation_most_freq(dataset)
	important_variants_list = pd.read_csv(args.important, compression = 'gzip')
	important_variants_list = important_variants_list.drop(important_variants_list.columns[0], axis = 1)

	ld_between_features(dataset, important_variants_list)


# measures of wilcoxon between variant and label and LD between variants
# (LD is not working properly) ** yet
def ld_between_features(dataset, important_variants_list):
	
	chisquared_df = pd.DataFrame(columns = ['variant', 'chisquared_p-value', 'chi2'])
	linkage_df = pd.DataFrame(columns = ['correlation', 'LD'])

	# measure LD between highest ranking variants
	important_variants = dataset[important_variants_list['variant']]
	dataset_len = len(dataset.index)

	count = 0
	ld_count = 0
	print('Computing metrics...')

	for variant1 in list(important_variants):
		print('Current variant: ' + str(variant1))

		# first, test correlation between variant 1 and labels
		# Applying X^2 

		# Making contigency table (GTs_no x 2)  
		table = pd.crosstab(important_variants[variant1], dataset['labels'])
		table = table.loc[:, [0, 1]]

		chi2, p, dof, ex = stats.chi2_contingency(table, correction = False)
		chisquared_df.loc[count] = [variant1, p, chi2]


		for variant2 in list(important_variants):

			if variant1 != variant2:

				# Secondly, measure LD between variants				
				# Can only use SNPs for LD
				if max(important_variants[variant1]) <= 2 and max(important_variants[variant2]) <= 2:

					# Frequencies of variant 1
					
					fA1 = (count_no(important_variants[variant1],0) / dataset_len \
						+ 1/2 * count_no(important_variants[variant1],1) / dataset_len)					
					fA2 = (count_no(important_variants[variant1],2) / dataset_len \
					    + 1/2 * count_no(important_variants[variant1], 1) / dataset_len)

					# frequencies of variant 2
					fB1 = (count_no(important_variants[variant2], 0) / dataset_len \
						+ 1/2 * count_no(important_variants[variant2], 1) / dataset_len) 
					fB2 = (count_no(important_variants[variant2], 2) / dataset_len \
						+ 1/2 * count_no(important_variants[variant2], 1) / dataset_len)

					# D = f(A1_B1) * f(A2_B2) - f(A1_B2) * f(A2_B1)
					D = ((fA1 * fB1) * (fA2 * fB2)) - ((fA1 * fB2) * (fA2 * fB1))

					# r² = D² / (fA1 * fA2 * fB1 * fB2)
					r_squared = (D ** 2) / (fA1 * fA2 * fB1 * fB2)

					linkage_df.loc[ld_count] = [str(variant1) + '-' + str(variant2), r_squared]
					ld_count += 1
					

		count += 1

	chisquared_df.to_csv('../../data/tests/LD_and_chisquared/chisquared.csv.gz', compression = 'gzip')
	linkage_df.to_csv('../../data/tests/LD_and_chisquared/LD.csv.gz', compression = 'gzip')

# normal count() was throwing errors...
def count_no(data, number_to_count):

	total = 0
	for value in data:
		if value == number_to_count:
			total += 1

	return total


def imputation_most_freq(dataset):

	dataset_len = len(list(dataset))
	dataset_list = np.array(list(dataset))

	# Imputing missing data with most frequent in column
	mf_imputation = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
	dataset = mf_imputation.fit_transform(dataset)
	dataset = dataset.astype(int)

	dataset = pd.DataFrame(dataset, columns = dataset_list )
	#print(dataset.head())

	return dataset

if __name__ == "__main__":
	main()

