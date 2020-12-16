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
model_file = 'dataset/'
#model_select = 'TargetFinder/'
cellline = str(sys.argv[2])
if cellline in ['GM12878','HUVEC','HeLa-S3','IMR90','K562','NHEK']:
    model_select = 'Targetfinder/'
elif cellline in ['FoeT','Mon','nCD4','tB','tCD4','tCD8']:
    model_select = 'DeepTACT/'
else:
    print('The model of the input cell line hasnot prepared yet, you can choose a similar cell line to make a prediction.')
    sys.exit(0)
enchrome = str(sys.argv[3])
enstart = str(sys.argv[4])
enend = str(sys.argv[5])
prchrome = str(sys.argv[6])
prstart = str(sys.argv[7])
prend = str(sys.argv[8])
enmid =  (int(enstart)+int(enend))//2
newenstart = enmid - 1500
newenend = newenstart + 2999
prmid = (int(prstart)+int(prend))//2
newprstart = prmid-1000
newprend = newprstart + 1999
distance = abs(prmid-enmid)
dis =  '%.4f' % math.log((2000000/distance),10)
model_filepath = model_file+model_select
enoldname = cellline+'|'+enchrome+':'+str(enstart)+'-'+str(enend)
proldname = cellline+'|'+prchrome+':'+str(prstart)+'-'+str(prend)
enname = cellline+'|'+enchrome+':'+str(newenstart)+'-'+str(newenend)
prname = cellline+'|'+prchrome+':'+str(newprstart)+'-'+str(newprend)

kmer = 4**kvalue

train_num = 1



fin1 = open('enhancer.bed','w')
fin2 = open('promoter.bed','w')
for i in range(2): #In order to be suitable for using SEEKR.
    fin1.write(enchrome+'\t'+str(newenstart)+'\t'+str(newenend)+'\t'+enname+'\n')
    fin2.write(prchrome+'\t'+str(newprstart)+'\t'+str(newprend)+'\t'+prname+'\n')
fin1.close()
fin2.close()

os.system("bedtools getfasta -fi ../hg19.fa -bed enhancer.bed -fo enhancer.fa")
os.system("bedtools getfasta -fi ../hg19.fa -bed promoter.bed -fo promoter.fa")
os.system("python3 ../seekr_py/src/kmer_counts.py enhancer.fa -o enhancer.txt -k {} -nb".format(kvalue))
os.system("python3 ../seekr_py/src/kmer_counts.py promoter.fa -o promoter.txt -k {} -nb".format(kvalue))



#generate data matrix
arrays = numpy.zeros((train_num, kmer*2))
labels = numpy.zeros(train_num)
distance = numpy.zeros(train_num)


fin1 = open('enhancer.txt','r')
fin2 = open('promoter.txt','r')
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

enhancer_vec = enhancer[0]
promoter_vec = promoter[0]
enhancer_vec = enhancer_vec.reshape((1,kmer))
promoter_vec = promoter_vec.reshape((1,kmer))
arrays[0] = numpy.column_stack((enhancer_vec,promoter_vec))
distance[0] = float(dis)

X_train = numpy.column_stack((arrays,distance))
print(X_train.shape[0],X_train.shape[1])


estimator = CatBoostClassifier(iterations = 1000,depth = 10,learning_rate = 0.1,logging_level = None,scale_pos_weight = 45)
estimator.load_model('{}/{}/best_model{}'.format(model_filepath,cellline,kvalue))

y_pred = estimator.predict(X_train)
y_proba_pred = estimator.predict_proba(X_train)[:,1]
if enchrome != prchrome:
    print('The two elements are not in the same chrosome, please recheck your input!')
else:
    print('For Promoter '+enoldname+', Enhancer '+proldname+' in cell line '+ cellline+' :')
    if y_pred[0] == 0:
        print('The two elements are predicted not to be interacted by EPBoost, the interaction prediction score is %.4f.'%y_proba_pred[0])
    else:
        print('The two elements are predicted interacted by EPBoost, the interaction prediction score is %.4f.'%y_proba_pred[0])
