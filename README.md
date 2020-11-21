# EPBoost
  EPBoost is a **quick** and **accurate** method to identify enhancer-promoter interactions using intrinsic features generated from genomic sequences. It exploits the kmer content counts of the sequences as inputs and trains and predicts with a CatBoost model. To evaluate the performance of the model we extracted interaction data of 12 cell lines from **TargetFinder** https://github.com/shwhalen/targetfinder and **DeepTACT** https://github.com/liwenran/DeepTACT. 


## **Training**
### **STEP1** <br>
* `python DataPrepare.py`<br>
   In this process we will pad input enhancers into *3000bp long* and promoters into *2000bp long*:<br>
   only **one** file is needed: ***dataset/TargetFinder(or DeepTACT)/celllinename/pairs.csv***<br>
   and **three** files will be produced: ***enhancers.bed***,   ***promoters.bed***,   ***train.csv***<br>
### **STEP2** <br>
* `python EPBoost.py k`<br>
   This is the training program, the _k_ determines the length of the kmer which can be ranged from 3 to 7, the imbalance ratio in training set and test set are both 1:20.
* `python EPBoost2.py k`<br>
   This is the training program to compare with DeepTACT, the _k_ determines the length of the kmer which can be ranged from 3 to 7, the imbalance ratio in training set is 1:20 and in test set is 1:5.
### **Note** <br>
   In  processes counting and normalizing the kmer contents, we basically adapted the code in ***seer_py*** which is originally from https://github.com/CalabreseLab/seekr.
   We use hg19.fa file as a reference genome which can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz.

## **Test**
* `python EPBoost_Test.py k model_cell test_cell`<br>
   This is the test program, the _k_ determines the length of the kmer which can be ranged from 3 to 7, the _model_cell_ defines the trained model we use for predicting, the _test_cell_ refers to the cell line we would like to make a prediction.


## **Requirements**


## **License**
  EPBoost is licensed under the MIT License - details can be found in the LICENSE.md file
## **Download**
  Please download the related profiles from https://github.com/ACatfromUSTC/EPBoost/. 
