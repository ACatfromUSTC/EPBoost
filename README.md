# EPBoost
A simple but accurate method to identify enhancer-promoter interactions using intrinsic features generated from genomic sequence.

#**Training
1. python DataPrepare.py
   In this process we will pad our enhancers into 3000-bp-long and promoters into 2000-bp-long, 
   1 file is needed in dataset/TargetFinder(or DeepTACT)/celllinename/pairs.csv
   3 files will be produced.
      **file1: enhancers.bed
      **file2: promoters.bed
      **file3: train.csv
#Test
