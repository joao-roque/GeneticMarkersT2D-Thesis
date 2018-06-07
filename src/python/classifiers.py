import pandas as pd
import numpy as np 
import sklearn as sk
import argparse
import sys
import matplotlib.pyplot as plt
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.preprocessing import Imputer
from sklearn import metrics
from sklearn.metrics import f1_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns 
import operator
import random
from sklearn import svm

def main():

	# Argument parser
	arg_parser = argparse.ArgumentParser(description = 'Feature Selection')
	arg_parser.add_argument('-dataset', help = 'Input file to classify or plot', type = str)
	arg_parser.add_argument('-label', help = 'Column Name of the labels variable')
	arg_parser.add_argument('-classifier', help = 'Input guide: "rf" for random forest,\
	 "glob" for a group of Extra Trees Classifier, "redu" for feature reduction/accuracy loss plot,\
	 "linear" for a SVM classifier.')
	arg_parser.add_argument('-imputation', help = 'If (y), it will perform imputation using the most frequent in column')
	arg_parser.add_argument('-plots', help = 'Input Guide: "frequencies" for frequency per run graph, \
		"redu" for accuracy loss for feature reduction')
	arg_parser.add_argument('-ranking', help = 'Make variables ranking (y)')
	arg_parser.add_argument('-trees', help = 'Test optimal nummber of trees (pick number)', type = int)
	arg_parser.add_argument('-top', help = 'Runs Extra trees with top variants (rank file)')
	arg_parser.add_argument('-redu_features', help = 'File of top important features')
	args = arg_parser.parse_args()

	if args.dataset:
		# read_csv adds an Unnamed columns of indexes
		print('Loading dataset...')
		dataset = pd.read_csv(args.dataset)#, compression = 'gzip')
		dataset = dataset.drop(dataset.columns[0], axis = 1)

		if args.imputation:

			print('Performing imputation...')
			dataset = imputation_most_freq(dataset)

			dataset.to_csv('../../data/full_dataset/imputed_clean_full_dataset.csv.gz', compression = 'gzip')
		
	if args.classifier:

		if args.classifier == 'glob':
			glob_of_forests(dataset, args)
		elif args.classifier == 'rf':
			(frequencies, f1s) = extra_trees(dataset, args)
		elif args.classifier == 'redu':
			if args.redu_features:
				redu(dataset, args)
			else:
				print('Redu function requires important features file.')
		elif args.classifier == 'linear':
			#if args.redu_features:
			linear_classifier(dataset, args)
			#else:
			#	print('A linear classifier requires top features file (-redu_features)')

	if args.plots:
		if args.plots  == 'frequencies':
			plot_glob(dataset)
		if args.plots == 'redu':
			random_data = pd.read_csv('../../data/plots/accuracy_by_feature_reduction/random_f1scores.csv.gz', compression = 'gzip')
			targeted_data = pd.read_csv('../../data/plots/accuracy_by_feature_reduction/targeted_f1scores.csv.gz', compression = 'gzip')
			plot_redu(random_data, targeted_data)

	if args.ranking:
		important_variants(dataset, args)

	if args.trees:
		dataset = imputation_most_freq(dataset)
		no_trees_test(dataset, args)

	if args.top:

		print('Performing imputation...')
		dataset = imputation_most_freq(dataset)

		important_variants_list = pd.read_csv(args.top, compression = 'gzip')		
		new_dataset = dataset[important_variants_list]

		new_dataset['labels'] = dataset['labels']

		(frequencies, f1s) = extra_trees(new_dataset, args)
		print(frequencies)

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

# runs extra trees classifier in dataset
def extra_trees(dataset, args):

	dataset_len = len(list(dataset))
	dataset_list = np.array(list(dataset))

	# Run several classifiers and find wich variants appear more
	f1_scores = []
	top_variants = np.array([])
	for run in range(0, 100):
		# Dividing train and test
		pre_X = dataset.drop(args.label, axis = 1)
		X_train, X_test, y_train, y_test = train_test_split(pre_X, dataset[args.label], test_size = 0.3)

		print('Fitting model number: '  + str(run))

		# Random Forest Model
		rf_model = ExtraTreesClassifier(100, n_jobs = -1)
		rf_model.fit(X_train, y_train)

		# To save model
		#joblib.dump(rf_model, '../../data/rf_model.pkl') 

		print('Predicting...')
		
		# Predicting
		# y_predicted = rf_model.predict(X_test)
		y_predicted = cross_val_predict(rf_model, X_test, y_test, cv=5)		
		# print(confusion_matrix(y_test, y_predicted))
		
		print('F1-score: ' + str(f1_score(y_test, y_predicted, average = 'micro')))

		# Saving f1 score and all indexes
		importances = rf_model.feature_importances_
		importances_indexes = np.argsort(importances)[::-1]
		importances_indexes = importances_indexes[0:100]

		f1_scores.append(f1_score(y_test, y_predicted, average = 'micro'))
		top_variants = np.concatenate((top_variants, dataset_list[importances_indexes]))


	# Results and important features
	f1s = mean(f1_scores)
	print('Final f1 score: ' + str(f1s))	 
	frequencies = Counter(top_variants)

	return (frequencies, f1s)

def glob_of_forests(dataset, args):

	frequencies_plot_info = pd.DataFrame(columns = ['variant', 'frequency', 'run_no'])

	# calling extra_trees for each 
	count = 0
	for runs_no in range(0,5):

		print('Glob of forests number ' + str(runs_no))
		(frequencies, f1s) = extra_trees(dataset,args)

		for key in frequencies.keys():
			frequencies_plot_info.loc[count] = [key, frequencies[key], runs_no]
			count += 1

	frequencies_plot_info.to_csv('../../data/plots/frequencies_of_variants/frequencies_ofvariant_per_run.csv.gz', compression = 'gzip')
	important_variants(frequencies_plot_info, args)

	plot_glob(frequencies_plot_info)	

# Extract important variants
def important_variants(dataset, args):

	# make a rank of variables from the average of frequencies
	variants_ranking = pd.DataFrame(columns = ['rank', 'variant', 'average_frequency']) 
	variants_dictionary = {}

	# it will repeat variants but since they're added to the dictionary the last value will be correct
	print('Looking for variants...')
	for current_variant in dataset['variant']:

		# select all
		variant_dataset = dataset.loc[dataset['variant'] == current_variant]
		frequencies_average = variant_dataset['frequency'].mean()
		variants_dictionary[current_variant] = frequencies_average

	count = 0
	dict_empty = False
	print('Saving variants in ranking...')
	while dict_empty == False:

		max_variant = max(variants_dictionary, key=variants_dictionary.get)
		variants_ranking.loc[count] = [count, max_variant, variants_dictionary[max_variant]]
		count +=1

		# Now the variant is deleted from the dictionary
		del variants_dictionary[max_variant]

		# exit condition
		if not variants_dictionary:
			dict_empty = True

	print('Saving ranking to file')
	variants_ranking.to_csv('../../data/plots/frequencies_of_variants/variants_ranking.csv.gz', compression = 'gzip')
	return variants_ranking


def plot_glob(frequencies_plot_info):

	print('Plotting...')

	frequencies_plot_info.head()

	plt.figure()
	plot = sns.swarmplot(x = "variant", y = "frequency", hue = "run_no", data = frequencies_plot_info)
	plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
	#plot.get_xaxis().set_visible(False)
	plt.xlabel('Variants')
	plt.tight_layout()
	plt.show()

# Testing optimal number of trees
def	no_trees_test(dataset, args):


	dataset_len = len(list(dataset))
	dataset_list = np.array(list(dataset))

	# Run several classifiers and find wich variants appear more
	scores_df = pd.DataFrame(columns = ['no_trees', 'f1', 'auc']) 
	for trees in range(1, args.trees + 1):
		
		# Dividing train and test
		pre_X = dataset.drop(args.label, axis = 1)
		X_train, X_test, y_train, y_test = train_test_split(pre_X, dataset[args.label], test_size = 0.3)

		print('Fitting model ' + str(trees) + ' ...')

		# Random Forest Model
		rf_model = ExtraTreesClassifier(trees, n_jobs = -1)
		rf_model.fit(X_train, y_train)

		# To save model
		#joblib.dump(rf_model, '../../data/rf_model.pkl') 

		print('Predicting...')
		
		# Predicting
		y_predicted = cross_val_predict(rf_model, X_test, y_test, cv=5)		
		#y_predicted = rf_model.predict(X_test)

		print(confusion_matrix(y_test, y_predicted))
		
		# Saving f1 and auc scores 
		fpr, tpr, thresholds = metrics.roc_curve(y_test, y_predicted, pos_label=2)
		scores_df.loc[trees] = [trees, f1_score(y_test, y_predicted, average = 'micro'),\
								   metrics.auc(fpr, tpr)]

	scores_df.to_csv('../../data/tests/accuracy_by_treesno/accuracy_by_treesno.csv.gz', compression = 'gzip')
	
	return scores_df

def plot_no_trees(scores):

	plt.figure()
	plt.scatter(scores['no_trees'], scores['f1'], color = 'r')
	#plt.scatter(scores['no_trees'], scores['auc'], color = 'g')
	plt.show()

	#g = sns.FacetGrid(, hue="time", col="sex", size=4,\
    #	              hue_kws={"marker": ["s", "D"]})
	#g.map(qqplot, "total_bill", "tip", s=40, edgecolor="w")
	#g.add_legend();

# run classifiers while reducing amount of features
# make it by random feature reducing or targeted reducing (leaving important features)
def redu(dataset, args):

	important_features_list = pd.read_csv(args.redu_features, compression = 'gzip')
	important_features_list = important_features_list['variant']
	features_list_len = len(important_features_list)

	# randomly
	# reduce 10 by 10
	
	count = 0
	print('Randomly selecting features...')
	random_data = pd.DataFrame(columns = ['features_num', 'f1-score'])
	for features_num in range(features_list_len, 10, -10):
		random_features = random.sample(list(important_features_list), k = features_num)
		random_features.append('labels')

		temporary_dataset = dataset[random_features]
		#print(temporary_dataset.head())
		(frequencies, f1s) = extra_trees(temporary_dataset, args)

		random_data.loc[count] = [features_num, f1s]
		count += 1

	random_data.to_csv('../../data/plots/accuracy_by_feature_reduction/random_f1scores.csv.gz', compression = 'gzip')
	
	
	# targeted
	print('Targeting top features...')
	count = 0
	targeted_data = pd.DataFrame(columns = ['features_num', 'f1-score'])
	for features_num in range(features_list_len, 10, -10):
		targeted_features = list(important_features_list.iloc[0:features_num])
		targeted_features.append('labels')

		temporary_dataset = dataset[targeted_features]
		#print(temporary_dataset.head())
		(frequencies, f1s) = extra_trees(temporary_dataset, args)

		targeted_data.loc[count] = [features_num, f1s]
		count += 1

	targeted_data.to_csv('../../data/plots/accuracy_by_feature_reduction/targeted_f1scores.csv.gz', compression = 'gzip')


def plot_redu(random_data, targeted_data):

	plt.figure()
	plt.title('Accuracy by feature reduction (number of features)')
	plt.xlabel('Number of features')
	plt.ylabel('f1-score')
	plt.scatter(random_data['features_num'], random_data['f1-score'], color = 'r')
	plt.scatter(targeted_data['features_num'], targeted_data['f1_score'], color = 'g')
	plt.legend(['Randomly selected features', 'Best features left'])
	plt.show()

def	linear_classifier(dataset, args):

	#important_features_list = pd.read_csv(args.redu_features, compression = 'gzip')
	#important_features_list = important_features_list['variant']
	#features_list_len = len(important_features_list)

	#targeted_features = list(important_features_list.iloc[1000:1005])
	#targeted_features.append('labels')

	#temporary_dataset = dataset[targeted_features]
	#print(temporary_dataset)

	#pre_X = temporary_dataset.drop(args.label, axis = 1)
	f1_scores = []
	for run in range(0,100):
		X_train, X_test, y_train, y_test = train_test_split(dataset, dataset[args.label], test_size = 0.3)

		svm_classifier = svm.SVC()
	
		print('Fitting model number: '  + str(run))
		svm_classifier.fit(X_train, y_train)

		print('Predicting...')
		#y_predicted = svm_classifier.predict(X_test)
		y_predicted = cross_val_predict(svm_classifier, X_test, y_test, cv=5)		
		print('F1 Score:' + str(f1_score(y_test, y_predicted, average = 'micro')))

		f1_scores.append(f1_score(y_test, y_predicted, average = 'micro'))

	f1s = mean(f1_scores)
	print('Final f1 score: ' + str(f1s))	
	frequencies = Counter(top_variants)

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

if __name__ == '__main__':
	main()
