# EPBoost
A **simple** but **accurate** method to identify enhancer-promoter interactions using intrinsic features generated from genomic sequence.

## **Training**
### **STEP1** <br>
* `python DataPrepare.py`<br>
   In this process we will pad our enhancers into *3000-bp-long* and promoters into *2000-bp-long*:<br>
   only **one** file is needed: ***dataset/TargetFinder(or DeepTACT)/celllinename/pairs.csv***<br>
   and **three** files will be produced: ***enhancers.bed***,   ***promoters.bed***,   ***train.csv***<br>
### **STEP2** <br>
* `python EPBoost.py k`<br>
   This is the training program, the _k_ determines the length of the kmer which can be ranged from 3 to 7, the imbalance ratio in training set and test set are both 1:20.
* `python EPBoost2.py k`<br>
   This is the training program to compare with DeepTACT, the _k_ determines the length of the kmer which can be ranged from 3 to 7, the imbalance ratio in training set is 1:20 and in test set is 1:5.
### **Note** <br>
   In percentage-normalization process, we basically adapted the code in ***seer_py*** which is originally from https://github.com/CalabreseLab/seekr.
   We use hg19.fa file as a reference genome which can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/.

## **Test**
* `python EPBoost_Test.py k`<br>
