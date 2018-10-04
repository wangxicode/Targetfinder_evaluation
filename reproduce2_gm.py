import sys,os
import numpy as np
import random
import pandas as pd

from random import randint
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_curve, auc, f1_score, precision_recall_curve
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import GradientBoostingClassifier

nonpredictors = ['enhancer_chrom', 'enhancer_start', 'enhancer_end', 'promoter_chrom', 'promoter_start', 'promoter_end', 'window_chrom', 'window_start', 'window_end', 'window_name', 'active_promoters_in_window', 'interactions_in_window', 'enhancer_distance_to_promoter', 'bin', 'label']

training_df = pd.read_hdf('training.gm.h5', 'training').set_index(['enhancer_name', 'promoter_name'])
predictors_df = training_df.drop(nonpredictors, axis = 1)
labels = training_df['label']

predictors_df_sorted = pd.concat([predictors_df[0:2113].reset_index().sort_values('promoter_name').set_index(['enhancer_name', 'promoter_name']),predictors_df[2113:44313].reset_index().set_index(['enhancer_name', 'promoter_name'])])

labels_sorted = pd.concat([pd.DataFrame(labels)[0:2113].reset_index().sort_values('promoter_name').set_index(['enhancer_name', 'promoter_name']),pd.DataFrame(labels)[2113:44313].reset_index().set_index(['enhancer_name', 'promoter_name'])])['label']

true_lab_group = pd.DataFrame(labels)[0:2113].groupby(['promoter_name']).size().tolist()

test_fold = []

for i in true_lab_group:
	rnum = randint(0,4)
	while(test_fold.count(rnum) >= 211*2):
		rnum = randint(0,4)
	test_fold = test_fold + [rnum for j in range(i)]

neg_lab = []

for i in range(5):
	neg_lab = neg_lab + [i for j in range(4220*2)]

random.shuffle(neg_lab)
test_fold = np.asarray(test_fold + neg_lab)

#for fset in ['epw','w','ep']:
#	for kern in ['linear','rbf','tree100','tree4k']:

print(0)

def main(argv = sys.argv):

	kern = argv[1]
	fset = argv[2]

	avg_auroc=0
	avg_auprc=0
	avg_f1=0

	for i in range(5):
		print(1)
		if kern=='rbf' or kern=='linear':
			clf=SVC(kernel=kern,probability=True)
		elif kern=='tree100':
			clf=GradientBoostingClassifier(n_estimators = 100, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 0)
		elif kern=='tree4k':
			clf=GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = 0)
		print(clf)
		if fset=='epw':
			f_fit = predictors_df_sorted.iloc[test_fold != i]
			f_test = predictors_df_sorted.iloc[test_fold == i]
		elif fset=='w':
			f_fit = predictors_df_sorted.iloc[test_fold != i, [i*3+2 for i in range(97)]+[291,294,297]]
			f_test = predictors_df_sorted.iloc[test_fold == i, [i*3+2 for i in range(97)]+[291,294,297]]
		elif fset=='ep':
			f_fit = predictors_df_sorted.iloc[test_fold != i, [i for i in range(298) if i not in [i*3+2 for i in range(97)]+[291,294,297]]]
			f_test = predictors_df_sorted.iloc[test_fold == i, [i for i in range(298) if i not in [i*3+2 for i in range(97)]+[291,294,297]]]

		l_fit = labels_sorted[test_fold != i]
		l_test = labels_sorted[test_fold == i]
			
		clf.fit(f_fit, l_fit)
		prob=clf.predict_proba(f_test)
		fpr, tpr, thresholds = roc_curve(l_test, prob[:, 1])
		roc_auc = auc(fpr, tpr)
		precision,recall,thresholds = precision_recall_curve(l_test, prob[:, 1])
		prc_auc = auc(recall,precision)
		f1=0
		for i in range(len(precision)):
			f1c=2*precision[i]*recall[i]/(precision[i]+recall[i])
			if f1c>f1:
				f1 = f1c
		avg_auroc=avg_auroc+roc_auc
		avg_f1=avg_f1+f1
		avg_auprc=avg_auprc+prc_auc
	print('reproduce2',fset,kern,avg_auroc/5,avg_auprc/5,avg_f1/5)
#	print('reproduce2 cv: {1} kern: {2} auroc: {3} auprc: {4} f1: {5}'.format(fset,kern,avg_auroc,avg_auprc,avg_f1))

main()
