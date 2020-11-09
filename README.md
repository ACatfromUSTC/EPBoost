# EPBoost
A **simple** but **accurate** method to identify enhancer-promoter interactions using intrinsic features generated from genomic sequence.

## **Training**
**STEP1** <br>
* `python DataPrepare.py`<br>
   In this process we will pad our enhancers into *3000-bp-long* and promoters into *2000-bp-long*:<br>
   only **one** file is needed: ***dataset/TargetFinder(or DeepTACT)/celllinename/pairs.csv***<br>
   and **three** files will be produced: ***enhancers.bed***,   ***promoters.bed***,   ***train.csv***<br>
**STEP2** <br>
* `python epnetfeature.py k`<br>

## **Test**
