# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 19:23:45 2020
@author: Wangzihang
"""


# system modules
import os
import time
import sys
import pandas as pd

# numpy
import numpy,random,math
# classifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from imblearn.ensemble import BalanceCascade


#EPBoost
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_curve, precision_score, f1_score, recall_score, precision_recall_curve, auc, average_precision_score
from sklearn import metrics, svm
from catboost import CatBoostClassifier
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

kvalue = int(sys.argv[1])
test_file = '/EPBoost/EPBoost/dataset/'
test_select = 'TargetFinder/'
model_cellline = str(sys.argv[2])
cellline = str(sys.argv[3])
test_filepath = test_file+test_select
os.system("bedtools getfasta -fi ../hg19.fa -bed {}enhancers.bed -fo {}enhancers.fa".format(test_filepath,test_filepath))
os.system("bedtools getfasta -fi ../hg19.fa -bed {}promoters.bed -fo {}promoters.fa".format(test_filepath,test_filepath))
os.system("python3 ../seekr_py/src/kmer_counts.py {}enhancers.fa -o {}enhancers.txt -k {} -nb".format(test_filepath,test_filepath,kvalue))
os.system("python3 ../seekr_py/src/kmer_counts.py {}promoters.fa -o {}promoters.txt -k {} -nb".format(test_filepath,test_filepath,kvalue))
kmer = 4**kvalue
print(kmer)
enhancers_num=0
promoters_num=0
positive_num=0
negative_num=0
train_num = 0
test_num = 0



fin1 = open(test_filepath+'enhancers.bed','r')
fin2 = open(test_filepath+'promoters.bed','r')
enhancers = []
promoters = []
for line in fin1:
	data = line.strip().split()
	enhancers.append(data[3])
	enhancers_num = enhancers_num + 1
for line in fin2:
	data = line.strip().split()
	promoters.append(data[3])
	promoters_num = promoters_num + 1


#generate index
fin3 = open(test_filepath+'{}/train.csv'.format(cellline),'r')
fout1 = open(test_filepath+'{}/training.txt','w')


for line in fin3:
	if line[0] == 'b':
		continue
	else:
		data = line.strip().split(',')
		enhancer_index = enhancers.index(data[5])
		promoter_index = promoters.index(data[9])
		fout1.write(str(enhancer_index)+'\t'+str(promoter_index)+'\t'+data[10]+'\t'+data[11]+'\n')
		train_num = train_num + 1
fout1.close()
print(train_num)

#generate data matrix
arrays = numpy.zeros((train_num, kmer*2))
labels = numpy.zeros(train_num)
distance = numpy.zeros(train_num)

fin = open(test_filepath+'{}/training.txt','r')

fin1 = open(test_filepath+'{}/enhancers.txt','r')
fin2 = open(test_filepath+'{}/promoters.txt','r')
df1=[]
df2=[]


for m,line in enumerate(fin1):
	data1 = line.split(',')
	data1 = numpy.array(data1,dtype = float)
	df1.append(data1)
enhancer = df1
for n,line in enumerate(fin2):
	data2 = line.split(',')
	data2 = numpy.array(data2,dtype = float)
	df2.append(data2)
promoter = df2





for i,line in enumerate(fin):
	data = line.strip().split()
	enhancer_vec = enhancer[int(data[0])]
	promoter_vec = promoter[int(data[1])]
	enhancer_vec = enhancer_vec.reshape((1,kmer))
	promoter_vec = promoter_vec.reshape((1,kmer))
	arrays[i] = numpy.column_stack((enhancer_vec,promoter_vec))
	distance[i] = float(data[3])
	labels[i] = int(data[2])




X_train = numpy.column_stack((arrays,distance))
print(X_train.shape[0])
y_train = labels
print(numpy.sum(y_train))




cv = StratifiedKFold(n_splits = 10, shuffle = True, random_state = 0)
acc       = []
auc       = []
aupr      = []
auprc     = []
recall    = []
precision = []
f1        = []
mcc       = []
i         = 0

def plot_AUROC(fpr,tpr):
    plt.figure(1, figsize=(8.5,8.5))
    
    plt.plot(fpr,tpr,label = 'Fold'+str(i)+': AUC = '+str('%.3f'%auc[-1]),linewidth = 2)
    ax=plt.gca();
    ax.spines['bottom'].set_linewidth(3);
    ax.spines['left'].set_linewidth(3);
    ax.spines['right'].set_linewidth(3);
    ax.spines['top'].set_linewidth(3);
    plt.tick_params(labelsize=20)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('arial') for label in labels]

    ax.tick_params(axis='x',width=2,colors='black')

    ax.tick_params(axis='y',width=2,colors='black')
    

    font1 = {'family' : 'arial',
    'weight' : 'normal',
    'size'   : 10,
    }
    font2 = {'family' : 'arial',
    'weight' : 'normal',
    'size'   : 30,
    }
    plt.xlabel(u'False Positive Rate', font2)
    plt.ylabel(u'True Positive Rate', font2)
    plt.title('ROC Curve', font2)
    plt.legend(prop = font1)
    plt.savefig('AUROC{}.png'.format(kvalue),dpi = 300)

def plot_AUPRC(rec,prec):
    plt.figure(2, figsize=(8.5,8.5))
    
    plt.plot(rec,prec,label = 'Fold'+str(i)+': AUPR = '+str('%.3f'%aupr[-1]),linewidth = 2)
    ax=plt.gca();
    ax.spines['bottom'].set_linewidth(3);
    ax.spines['left'].set_linewidth(3);
    ax.spines['right'].set_linewidth(3);
    ax.spines['top'].set_linewidth(3);
    plt.tick_params(labelsize=20)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('arial') for label in labels]

    ax.tick_params(axis='x',width=2,colors='black')

    ax.tick_params(axis='y',width=2,colors='black')
    

    font1 = {'family' : 'arial',
    'weight' : 'normal',
    'size'   : 10,
    }
    font2 = {'family' : 'arial',
    'weight' : 'normal',
    'size'   : 30,
    }
    plt.title('Precision/Recall Curve', font2)
    plt.xlabel(u'Recall', font2) 
    plt.ylabel(u'Precision', font2)
    plt.legend(prop = font1)
    plt.savefig('AUPRC{}.png'.format(kvalue),dpi = 300)

max_num = 0

for train,test in cv.split(X_train, y_train):
        i+=1
        print('validation:', i)
        estimator = CatBoostClassifier(iterations = 1000,depth = 10,learning_rate = 0.1,logging_level = None,scale_pos_weight = 45)
        #estimator = svm.SVC(kernel = 'rbf',C = 10, gamma = 0.012)
        #estimator = lgb.LGBMClassifier(is_unbalance = True, learning_rate = 0.012)
        estimator.load_model('{}{}/best_model3'.format(test_filepath,model_cellline))

        y_pred = estimator.predict(X_train[test,:])
        y_proba_pred = estimator.predict_proba(X_train[test,:])[:,1]
        TP = numpy.sum(numpy.logical_and(numpy.equal(y_train[test],1),numpy.equal(y_pred,1)))
        FP = numpy.sum(numpy.logical_and(numpy.equal(y_train[test],0),numpy.equal(y_pred,1)))
        TN = numpy.sum(numpy.logical_and(numpy.equal(y_train[test],0),numpy.equal(y_pred,0)))
        FN = numpy.sum(numpy.logical_and(numpy.equal(y_train[test],1),numpy.equal(y_pred,0)))

        accuracy = (TP+TN)/(TP+FP+TN+FN)
        acc.append(accuracy)
        fpr, tpr, th = metrics.roc_curve(y_train[test],y_proba_pred ,pos_label=1)
        auc.append(metrics.auc(fpr, tpr))
        plot_AUROC(fpr,tpr)
        aupr.append(metrics.average_precision_score(y_train[test],y_proba_pred))
        prec, rec, thres = metrics.precision_recall_curve(y_train[test],y_proba_pred ,pos_label=1)
        auprc.append(metrics.auc(rec, prec))

        plot_AUPRC(rec,prec)
        recall.append(metrics.recall_score(y_train[test],y_pred))
        precision.append(metrics.precision_score(y_train[test],y_pred))
        f1.append(metrics.f1_score(y_train[test],y_pred))
        m_c_c = (TP*TN - FP*FN)/(math.sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP)))
        mcc.append(m_c_c)



print ('Epnet,kmer feature selection using EPBoost:')
print ('acc:',numpy.mean(acc),numpy.std(acc))
print ('auc',numpy.mean(auc),numpy.std(auc))
print ('aupr',numpy.mean(aupr),numpy.std(aupr))
print ('auprc',numpy.mean(auprc),numpy.std(auprc))
print ('recall:',numpy.mean(recall),numpy.std(recall))
print ('precision',numpy.mean(precision),numpy.std(precision))
print ('f1',numpy.mean(f1),numpy.std(f1))
print ('mcc:',numpy.mean(mcc),numpy.std(mcc))
