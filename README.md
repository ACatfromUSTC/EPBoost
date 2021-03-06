# EPBoost
  EPBoost is a **quick** and **accurate** method to identify enhancer-promoter interactions using intrinsic features generated from genomic sequences. It exploits the kmer content counts of the sequences as inputs and trains and predicts with a CatBoost model.
  
## **Material**
 To evaluate the performance of the model we extracted interaction data of 12 cell lines from **TargetFinder** https://github.com/shwhalen/targetfinder and **DeepTACT** https://github.com/liwenran/DeepTACT. 
## **Usage** 

### **FilePath**
  When using EPBoost, the actual filepaths should be set properly. Take cell line NHEK as an example:<br>
  > **EPBoost** <br>
  >> **EPBoost_Test.py** <br>
  >> **Predict.py** <br>
  >> **DataPrepare.py** <br>
  >> **EPBoost_Train.py** <br>
  >> **EPBoost2_Train.py** <br>  
  >> **dataset** <br>
  >>> **DeepTACT** <br>
  >>> **TargetFinder** <br>
  >>>> **NHEK** <br>
      
  
### **Training**
#### **STEP1** <br>
* `$ python Dataprepare.py cell_line`<br>
   In this process we will pad input enhancers into *3000bp long* and promoters into *2000bp long* , the _cell_line_ refers to the name of the cell line:<br>
   only **one** file is needed: ***dataset/TargetFinder(or DeepTACT)/celllinename/pairs.csv***<br>
   and **three** files will be produced: ***enhancers.bed***,   ***promoters.bed***,   ***train.csv***<br>
##### **Example**
* `$ python Dataprepare.py NHEK` <br>
#### **STEP2** <br>
* `$ python EPBoost_Train.py k cell_line`<br>
   This is the training program, a 10-fold validation is also included. The _k_ determines the length of the kmer which can be ranged from 3 to 7, the _cell_line_ refers to the name of the cell line, the imbalance ratio in training set and test set are both 1:20. In the process, profiles of _enhancers.bed_, _promoters.bed_, _train.csv_ are needed and _enhancers.fa_, _promoters.fa_, _enhancers.txt_, _promoters.txt_, _training.txt_ are intermediate processing files. At last, a best_model will be generated and saved.
* `$ python EPBoost2_Train.py k cell_line`<br>
   This is the training program to compare with DeepTACT, a 10-fold validation is also included. The _k_ determines the length of the kmer which can be ranged from 3 to 7,  the _cell_line_ refers to the name of the cell line, the imbalance ratio in training set is 1:20 and in test set is 1:5. In the process, profiles of _enhancers.bed_, _promoters.bed_, _train.csv_ are needed and _enhancers.fa_, _promoters.fa_, _enhancers.txt_, _promoters.txt_, _training.txt_ are intermediate processing files. At last, a best_model will be generated and saved.
##### **Example**
* `$ python EPBoost_Train.py 3 NHEK` <br>
#### **Note** <br>
   In  processes counting and normalizing the kmer contents, we basically adapted the code in ***seer_py*** which is originally from https://github.com/CalabreseLab/seekr.
   We use hg19.fa file as a reference genome which can be downloaded by <br>
* `$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz`<br>
 *  `$ gzip -d hg19.fa.gz`<br>

### **Test**
* `$ python EPBoost_Test.py k model_cell_line test_cell_line`<br>
   This is the test program, the _k_ determines the length of the kmer which can be ranged from 3 to 7, the _model_cell_line_ defines the trained model we use for predicting, the _test_cell_line_ refers to the cell line we would like to make a prediction. eg:`$ python EPBoost_Test.py 3 GM12878 NHEK`
##### **Example**
* `$ python EPBoost_Test.py 3 NHEK GM12878` <br>
### **Predict**
* `$ python Predict.py k cell_line enchrome enstart enend prchrome prstart prend`<br>
   This is the predicting program, the _k_ determines the length of the kmer which can be ranged from 3 to 7 (here we provide a model setting k at 3), the _cell_line_ defines the trained model we use for predicting, the _enchrome_, _enstart_, _enend_, _prchrome_, _prstart_, _prend_ refer to the locations of the enhancer and promoter we would like to make a prediction, respectively.
#### **Example**
* `$ python3 Predict.py 3 NHEK chr1 3399800 3400600 chr1 3541200 3542000` <br>
**Output:** For Promoter NHEK|chr1:3399800-3400600, Enhancer NHEK|chr1:3541200-3542000 in cell line NHEK :
The two elements are predicted interacted by EPBoost, the interaction prediction score is 0.99766.
* `$ python3 Predict.py 3 NHEK chr1 10000000 10001000 chr1 10004000 10005000` <br>
**Output:** For Promoter NHEK|chr1:10000000-10001000, Enhancer NHEK|chr1:10004000-10005000 in cell line NHEK :
The two elements are predicted not to be interacted by EPBoost, the interaction prediction score is 0.0001.
* `$ python3 Predict.py 3 NHEK chr1 10000000 10001000 chr2 10004000 10005000` <br>
**Output:** The two elements are not in the same chrosome, please recheck your input!

## **Requirements**
* Python (run on 3.6.8)
* scikit-learn (run on 0.21.3)
* numpy (run on 1.16.2)
* bedtools (run on 2.28.0)
* catboost (run on 0.20)
* matplotlib (run on 3.1.2)
* tqdm (run on 4.38.0)

## **License**
  EPBoost is licensed under the MIT License - details can be found in the LICENSE.md file
## **Download**
  Please download the related profiles from https://github.com/ACatfromUSTC/EPBoost/. 
