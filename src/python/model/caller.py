import model
import classify

# constants
path_to_topgenes = '../../../data/knownvariants/risk_genes.csv'
path_to_dataset = '../../../data/full_dataset/imputed_clean_full_dataset.csv.gz'
path_to_test = '../../../data/full_dataset/data_for_test.csv.gz'
path_to_risk_genes_list = '../../../data/knownvariants/genes_list.csv'
path_to_risk_genes_list_test = '../../../data/knownvariants/genes_list_test.csv'
path_to_new_features = '../../../data/full_dataset/full_data_new_features.csv'
#path_to_dataset = '../../../data/full_dataset/full_treated_dataset_no_top100.csv.gz'

# call model
#run = model.model(path_to_dataset, path_to_topgenes)
#run.topGenes()

# call classifier
classifiers = classify.Classify(path_to_new_features)
(frequencies, f1s) = classifiers.classify_it()

print(frequencies)