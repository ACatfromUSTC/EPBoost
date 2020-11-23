# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:13:35 2020
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
cellline = str(sys.argv[2])
enchrome = str(sys.argv[3])
enstart = str(sys.argv[4])
enend = str(sys.argv[5])
prchrome = str(sys.argv[6])
prstart = str(sys.argv[7])
prend = str(sys.argv[8])
test_filepath = test_file+test_select+cellline
enname = enchrome+':'+enstart+'-'+enend
prname = prchrome+':'+prstart+'-'+prend
os.system("bedtools getfasta -fi ../hg19.fa -bed enhancer.bed -fo enhancer.fa")
os.system("bedtools getfasta -fi ../hg19.fa -bed promoter.bed -fo promoter.fa")
os.system("python3 ../seekr_py/src/kmer_counts.py enhancer.fa -o enhancer.txt -k {} -nb")
os.system("python3 ../seekr_py/src/kmer_counts.py promoter.fa -o promoter.txt -k {} -nb")
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


estimator = CatBoostClassifier(iterations = 1000,depth = 10,learning_rate = 0.1,logging_level = None,scale_pos_weight = 45)
estimator.load_model('{}{}/best_model3'.format(test_filepath,model_cellline))

y_pred = estimator.predict(X_train[test,:])
y_proba_pred = estimator.predict_proba(X_train[test,:])[:,1]
