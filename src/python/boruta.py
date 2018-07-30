import importlib.util
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score


# import Boruta module
spec = importlib.util.spec_from_file_location("BorutaPy", "../../external/boruta_py/boruta/boruta_py.py")
boruta = importlib.util.module_from_spec(spec)
spec.loader.exec_module(boruta)


#data
print('Loading data...')
#dataset = pd.read_csv('../../data/full_dataset/imputed_clean_full_dataset.csv.gz', compression = 'gzip')
dataset = pd.read_csv('../../data/full_dataset/full_data_new_features_top_50.csv')
dataset = dataset.drop(dataset.columns[0], axis = 1)

labels = dataset['labels']
dataset = dataset.drop('labels', axis = 1)

# initializing tree
print('Building Extra Tree...')
rf = ExtraTreesClassifier(bootstrap=False, class_weight=None, criterion='gini',
           max_depth=None, max_features=None, max_leaf_nodes=20,
           min_samples_split=4, min_weight_fraction_leaf=0.0,
           n_estimators=50, n_jobs=-1, oob_score=False, random_state=None,
           verbose=0, warm_start=False)


# define Boruta feature selection method
print('Defining Boruta...')
feat_selector = boruta.BorutaPy(rf, n_estimators='auto', verbose=5, random_state=1)

# find all relevant features - 5 features should be selected
print('Fitting...')
feat_selector.fit(dataset.as_matrix(), labels.as_matrix())

# check selected features - first 5 features are selected
feat_selector.support_

# check ranking of features
feat_selector.ranking_

print(dataset.columns[feat_selector.ranking_.astype(int)])
# call transform() on X to filter it down to selected features
#X_filtered = feat_selector.transform(X)

X_filtered = feat_selector.transform(dataset.as_matrix())

rf_model = ExtraTreesClassifier(bootstrap=False, class_weight=None, criterion='gini',
           max_depth=None, max_features=None, max_leaf_nodes=20,
           min_samples_split=4, min_weight_fraction_leaf=0.0,
           n_estimators=50, n_jobs=-1, oob_score=False, random_state=None,
           verbose=0, warm_start=False)

y_predicted = cross_val_predict(rf_model, X_filtered, labels, cv=5)

print('F1-score: ' + str(f1_score(labels, y_predicted)) + '\nAccuracy: ' + \
 str(accuracy_score(y_pred = y_predicted, y_true = labels)))
