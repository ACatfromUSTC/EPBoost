# -*- coding: utf-8 -*-
"""
Created on Tue Dec 3 23:15:03 2019

@author: Wangzihang
"""

import math
names = ['GM12878','HUVEC','HeLa-S3','IMR90','K562','NHEK','FoeT','Mon','nCD4','tB','tCD4','tCD8']
cmd = '/EPBoost/EPBoost/dataset/TargetFinder/' # or '/EPBoost/EPBoost/dataset/DeepTACT/'

#can be 'GM12878','HUVEC','HeLa-S3','IMR90','K562','NHEK','FoeT','Mon','nCD4','tB','tCD4','tCD8'
cellline = names[5]
cellline_dir = names[5]+'/'
fin = open(cmd+cellline_dir+'pairs.csv','r')
fout = open(cmd+cellline_dir+'train.csv','w')
fout1 = open(cmd+cellline_dir+'enhancers.bed','w',newline = '')
fout2 = open(cmd+cellline_dir+'promoters.bed','w',newline = '')
fout.write('bin'+','+'enhancer_chrom'+','+'enhancer_start'+','+'enhancer_end'+','+'enhancer_name'+','+
           'promoter_chrom'+','+'promoter_start'+','+'promoter_end'+','+'promoter_name'+','+'label'+','+'distance'+'\n')




enhancers = []
promoters = []
#generate enhances
def generen(chrom,start,newstart,label):
    enlist   = []
    enstart = newstart
    enend   = newstart + 2999
    en      = str(chrom)+':'+str(enstart)+'-'+str(enend)
    enlist.append(en)
    return(enlist)

#generate promoters
def generpr(chrom,start,newstart,label):
    
    prlist = []
    prstart = newstart
    prend   = newstart + 1999
    pr      = str(chrom)+':'+str(prstart)+'-'+str(prend)
    prlist.append(pr)        
    return(prlist)
        
#for cell lines in TargetFinder
for line in fin:
    if line[0] == 'b':
        continue
    else:
        data  = line.strip().split(',')
        label = data[7]
        distancebin = data[0]+','+data[1].strip()
        distance = int(data[3])
        logdis = '%.4f' % math.log((2000000/distance),10)
        #enhancer sequnce information
        en_chrom = data[2]
        en_start = int(data[6])
        en_end   = int(data[4])
        en_len = en_end - en_start
        en_mid   = math.ceil(int((en_start+en_end)/2))
        en_newend = en_mid+1500
        en_newstart = en_mid-1500 
        
        #promoter sequence information
        pr_chrom = data[8]
        pr_start = int(data[11])
        pr_end   = int(data[9])
        pr_len = pr_end - pr_start
        pr_mid = math.ceil(int((pr_start+pr_end)/2))
        pr_newend = pr_mid+1000
        pr_newstart = pr_mid-1000
        
        #generate new enhancers
        en = generen(en_chrom,en_start,en_newstart,label)
        #generate new promoters
        pr = generpr(pr_chrom,pr_start,pr_newstart,label)
        
        #generate pairs
        for item1 in en:
            data1 = item1.split(':')
            data11 = data1[1].split('-')
            for item2 in pr:
                data2 = item2.split(':')
                data22 = data2[1].split('-')
                fout.write(distancebin+','+en_chrom+','+str(data11[0])+','+str(data11[1])+','+cellline+'|'+item1+','+pr_chrom+','+str(data22[0])+','+str(data22[1])+','+cellline+'|'+item2+','+label+','+str(logdis)+'\n')


        #generate beds
        for item1 in en:
            if item1 in enhancers:
                pass
            else:
                enhancers.append(item1)
                data1 = item1.split(':')
                data11 = data1[1].split('-')
                fout1.write(str(data1[0])+'\t'+str(data11[0])+'\t'+str(data11[1])+'\t'+cellline+'|'+item1+'\n')


        for item2 in pr:
            if item2 in promoters:
                pass
            else:
                promoters.append(item2)
                data2 = item2.split(':')
                data22 = data2[1].split('-')
                fout2.write(str(data2[0])+'\t'+str(data22[0])+'\t'+str(data22[1])+'\t'+cellline+'|'+item2+'\n')

fout.close()
fout1.close()
fout2.close()


#for cell lines in DeepTACT

'''for line in fin:
    if line[0] == 'e':
        continue
    else:
        data  = line.strip().split(',')
        label = data[8]
        #distancebin = data[0]+','+data[1].strip()
        corr = '1'
        mid1 = int((int(data[1])+int(data[2]))/2)
        mid2 = int((int(data[5])+int(data[6]))/2)
        distance = int(abs(mid1-mid2))
        logdis = '%.4f' % math.log((2000000/distance),10)
        #enhancer sequnce information
        en_chrom = data[0]
        en_start = int(data[1])
        en_end   = int(data[2])
        en_len = en_end - en_start
        en_mid   = math.ceil(int((en_start+en_end)/2))
        en_newend = en_mid+1499
        en_newstart = en_mid-1500 
        
        #promoter sequence information
        pr_chrom = data[4]
        pr_start = int(data[5])
        pr_end   = int(data[6])
        pr_len = pr_end - pr_start
        pr_mid = math.ceil(int((pr_start+pr_end)/2))
        pr_newend = pr_mid+999
        pr_newstart = pr_mid-1000
        
        #generate new enhancers
        en = generen(en_chrom,en_start,en_newstart,label)
        #generate new promoters
        pr = generpr(pr_chrom,pr_start,pr_newstart,label)
        
        #generate pairs
        for item1 in en:
            data1 = item1.split(':')
            data11 = data1[1].split('-')
            for item2 in pr:
                data2 = item2.split(':')
                data22 = data2[1].split('-')
                fout.write(corr+','+en_chrom+','+str(data11[0])+','+str(data11[1])+','+cellline+'|'+item1+','+pr_chrom+','+str(data22[0])+','+str(data22[1])+','+cellline+'|'+item2+','+label+','+str(logdis)+'\n')


        #generate beds
        for item1 in en:
            if item1 in enhancers:
                pass
            else:
                enhancers.append(item1)
                data1 = item1.split(':')
                data11 = data1[1].split('-')
                fout1.write(str(data1[0])+'\t'+str(data11[0])+'\t'+str(data11[1])+'\t'+cellline+'|'+item1+'\n')


        for item2 in pr:
            if item2 in promoters:
                pass
            else:
                promoters.append(item2)
                data2 = item2.split(':')
                data22 = data2[1].split('-')
                fout2.write(str(data2[0])+'\t'+str(data22[0])+'\t'+str(data22[1])+'\t'+cellline+'|'+item2+'\n')

fout.close()
fout1.close()
fout2.close()'''
                
        
        
