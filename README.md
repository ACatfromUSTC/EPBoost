# EPBoost
  EPBoost is a **quick** and **accurate** method to identify enhancer-promoter interactions using intrinsic features generated from genomic sequences. It exploits the kmer content counts of the sequences as inputs and trains and predicts with a CatBoost model.
  
## **Material**
 To evaluate the performance of the model we extracted interaction data of 12 cell lines from **TargetFinder** https://github.com/shwhalen/targetfinder and **DeepTACT** https://github.com/liwenran/DeepTACT. 
## **Usage** 
### **Training**
#### **STEP1** <br>
* `python DataPrepare.py`<br>
   In this process we will pad input enhancers into *3000bp long* and promoters into *2000bp long*:<br>
   only **one** file is needed: ***dataset/TargetFinder(or DeepTACT)/celllinename/pairs.csv***<br>
   and **three** files will be produced: ***enhancers.bed***,   ***promoters.bed***,   ***train.csv***<br>
#### **STEP2** <br>
* `python EPBoost_Train.py k`<br>
   This is the training program, a 10-fold validation is also included. The _k_ determines the length of the kmer which can be ranged from 3 to 7, the imbalance ratio in training set and test set are both 1:20. In the process, profiles of _enhancers.bed_, _promoters.bed_, _train.csv_ are needed and _enhancers.fa_, _promoters.fa_, _enhancers.txt_, _promoters.txt_, _training.txt_ are intermediate processing files. At last, a best_model will be generated and saved.
* `python EPBoost2_Train.py k`<br>
   This is the training program to compare with DeepTACT, a 10-fold validation is also included. The _k_ determines the length of the kmer which can be ranged from 3 to 7, the imbalance ratio in training set is 1:20 and in test set is 1:5. In the process, profiles of _enhancers.bed_, _promoters.bed_, _train.csv_ are needed and _enhancers.fa_, _promoters.fa_, _enhancers.txt_, _promoters.txt_, _training.txt_ are intermediate processing files. At last, a best_model will be generated and saved.
#### **Note** <br>
   In  processes counting and normalizing the kmer contents, we basically adapted the code in ***seer_py*** which is originally from https://github.com/CalabreseLab/seekr.
   We use hg19.fa file as a reference genome which can be downloaded by <br>
   `$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz`<br>
   `$ tar zxvf chromFa.tar.gz`<br>
   `$ cat chr*.fa > hg19.fa`<br>

### **Test**
* `python EPBoost_Test.py k model_cell_line test_cell_line`<br>
   This is the test program, the _k_ determines the length of the kmer which can be ranged from 3 to 7, the _model_cell_line_ defines the trained model we use for predicting, the _test_cell_line_ refers to the cell line we would like to make a prediction.

### **Predict**
* `python Predict.py k cell_line enchrome enstart enend prchrome prstart prend`<br>
   This is the predicting program, the _k_ determines the length of the kmer which can be ranged from 3 to 7 (here we provide a model with ), the _cell_line_ defines the trained model we use for predicting, the _enchrome_, _enstart_, _enend_, _prchrome_, _prstart_, _prend_ refer to the locations of the enhancer and promoter we would like to make a prediction, respectively.
## **Requirements**
* Python (run on 3.6.8)
* scikit-learn (run on 0.21.3)
* numpy (run on 1.16.2)
* bedtools (run on 2.28.0)
* catboost (run on 0.20)
* matplotlib (run on 3.1.2)

## **License**
  EPBoost is licensed under the MIT License - details can be found in the LICENSE.md file
## **Download**
  Please download the related profiles from https://github.com/ACatfromUSTC/EPBoost/. 
