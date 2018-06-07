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


class Classify(object):
	
	"""docstring for ClassName"""
	def __init__(self, dataset):
		self.dataset = pd.read_csv(dataset)
		self.dataset = self.dataset.drop(self.dataset.columns[0], axis = 1)
		
	def classify_it(self):

		dataset_len = len(list(self.dataset))
		dataset_list = np.array(list(self.dataset))

		# Run several classifiers and find wich variants appear more
		f1_scores = []
		top_variants = np.array([])
		for run in range(0, 100):
			# Dividing train and test
			pre_X = self.dataset.drop('labels', axis = 1)
			X_train, X_test, y_train, y_test = train_test_split(pre_X, self.dataset['labels'], test_size = 0.3)

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
			
			# Saving f1 score and all indexes
			importances = rf_model.feature_importances_
			importances_indexes = np.argsort(importances)[::-1]
			importances_indexes = importances_indexes[0:100]

			print('F1-score: ' + str(f1_score(y_test, y_predicted, average = 'micro')))
			f1_scores.append(f1_score(y_test, y_predicted, average = 'micro'))
			top_variants = np.concatenate((top_variants, dataset_list[importances_indexes]))


		# Results and important features
		f1s = self.mean(f1_scores)
		print('Final f1 score: ' + str(f1s))	 
		frequencies = Counter(top_variants)

		return(frequencies, f1s)
	

	def mean(self, numbers):
		return float(sum(numbers)) / max(len(numbers), 1)



	def plot_glob(self, frequencies_plot_info):

		print('Plotting...')

		frequencies_plot_info.head()

		plt.figure()
		plot = sns.swarmplot(x = "variant", y = "frequency", hue = "run_no", data = frequencies_plot_info)
		plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
		#plot.get_xaxis().set_visible(False)
		plt.xlabel('Variants')
		plt.tight_layout()
		plt.show()