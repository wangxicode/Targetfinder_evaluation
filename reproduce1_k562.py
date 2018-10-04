import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime

from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import roc_curve, auc, f1_score, precision_recall_curve

def count_chr(df,train,test):
    print(train)
    print(test)
    y=df['label'].values
    chr1=df['promoter_chrom'].values
    z=y[train]
    c1=chr1[train]
    ct={}
    np,nn=0,0
    for i in range(len(z)):
        if z[i]==1:
            np=np+1
        if z[i]==0:
            nn=nn+1
        if c1[i] in ct:
            ct[c1[i]]=ct[c1[i]]+1
        else:
            ct[c1[i]]=1
    print('train np:',np,' nn:',nn)
    for x in ct.keys():
        print(x,ct[x])
    z=y[test]
    c1=chr1[test]
    ct={}
    np,nn=0,0
    for i in range(len(z)):
        if z[i]==1:
            np=np+1
        if z[i]==0:
            nn=nn+1
        if c1[i] in ct:
            ct[c1[i]]=ct[c1[i]]+1
        else:
            ct[c1[i]]=1
    print('test np:',np,' nn:',nn)
    for x in ct.keys():
        print(x,ct[x])
                

def main(argv=sys.argv):
    cvset=eval(argv[1])   # 0..4
    fset=argv[2]          # epw, ep, w
    shuf=eval(argv[3])          # 0 1   1=shuffle
    kern=argv[4]          # linear, rbf
    rseed=eval(argv[5])         # random seed
    nonpredictors = ['enhancer_chrom', 'enhancer_start', 'enhancer_end', 'promoter_chrom', 'promoter_start', 'promoter_end', 'window_chrom', 'window_start', 'window_end', 'window_name', 'active_promoters_in_window', 'interactions_in_window', 'enhancer_distance_to_promoter', 'bin', 'label']

    training_df = pd.read_hdf('training.k562.h5', 'training').set_index(['enhancer_name', 'promoter_name'])

    training_df_sorted = training_df.reset_index().sort_values(by=['promoter_chrom','promoter_start']).set_index(['enhancer_name', 'promoter_name'])

    predict_df_epw = training_df_sorted.drop(nonpredictors, axis = 1)
#    predict_df_w = training_df_sorted.drop(nonpredictors, axis = 1).iloc[:,[i*3+2 for i in range(136)]]
#    predict_df_ep = training_df_sorted.drop(nonpredictors, axis = 1).iloc[:,[i for i in range(408) if (i+1)%3!=0]]
    predict_df_w = training_df_sorted.drop(nonpredictors, axis = 1).iloc[:,[i*3+2 for i in range(97)]+[291,294,297]]
    predict_df_ep = training_df_sorted.drop(nonpredictors, axis = 1).iloc[:,[i for i in range(298) if i not in [i*3+2 for i in range(97)]+[291,294,297]]]

    if fset=='epw':
        x=predict_df_epw.values
    if fset=='w':
        x=predict_df_w.values
    if fset=='ep':
        x=predict_df_ep.values
#    y=training_df['label'].values
    y=training_df_sorted['label'].values
#    print(x[0,:])
#    print(x[1,:])
#    print(y[0])
#    print(y[1])
#    cut=list(range(0,1000))+list(range(len(y)-1000,len(y)))
#    print(cut)
#    xc=x[cut]
#    yc=y[cut]
#    print(xc)
#    print(yc)
#    x=xc
#    y=yc

    print('training:')
    if kern=='rbf' or kern=='linear':
        clf=SVC(kernel=kern,probability=True)
    elif kern=='tree100':
        clf=GradientBoostingClassifier(n_estimators = 100, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = rseed)
    elif kern=='tree4k':
        clf=GradientBoostingClassifier(n_estimators = 4000, learning_rate = 0.1, max_depth = 5, max_features = 'log2', random_state = rseed)

    if shuf==0:
        skf = StratifiedKFold(y, 5,random_state=rseed,shuffle=False)
    if shuf==1:
        skf = StratifiedKFold(y, 5,random_state=rseed,shuffle=True)
    avg_auroc=0
    avg_auprc=0
    avg_f1=0
    icv=0

#    for train,test in skf:
#        count_chr(training_df_sorted,train,test)
    for train,test in skf:
#           print(datetime.datetime.now())
#           print("%s %s" % (train, test))
            clf.fit(x[train],y[train])
            prob=clf.predict_proba(x[test])
#           print(y[test],prob)
            fpr, tpr, thresholds = roc_curve(y[test], prob[:,-1])
            roc_auc = auc(fpr, tpr)
            precision,recall,thresholds =precision_recall_curve(y[test], prob[:,-1])
#        print('tpr:',tpr)
#        print('fpr:',fpr)
#        print('precision:',precision)
#        print('recall:',recall)
            prc_auc = auc(recall,precision)
            f1=0
            for i in range(len(precision)):
                f1c=2*precision[i]*recall[i]/(precision[i]+recall[i])
                if f1c>f1:
                    f1=f1c
#        print(roc_auc)
#        print('*',roc_auc,prc_auc,f1)
            print('mb3 cv: {0} fset: {1} shuf: {2} kern: {3} rseed: {4} auroc: {5} auprc: {6} f1: {7}'.format(icv,fset,shuf,kern,rseed,roc_auc,prc_auc,f1))
            avg_auroc=avg_auroc+roc_auc
            avg_f1=avg_f1+f1
            avg_auprc=avg_auprc+prc_auc
#            print(datetime.datetime.now())
         
    print('avg: auroc: {0} auprc: {1} f1: {2}'.format(avg_auroc/5,avg_auprc/5,avg_f1/5))


main()
