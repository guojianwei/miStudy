# -*- coding: utf-8 -*-
"""
Created on Thu Dec 08 21:13:00 2016

@author: G2
#refGene 来自ucsc 可以定位mRNA外显子位置
"""

import os 
import pandas as pd
'''
选择hg18用作测试程序的正确性
< D:\data\webdata\refGene38\refGene, data\sample_mrna_seq.txt
>r'data\prodata\hg38\mrna_exon_locations\\
'''
refgene_file = r'D:\data\webdata\refGene38\refGene.txt'
out_folder = r'data\prodata\hg38'
df = pd.read_csv(r'data\sample_mrna_seq.txt')
#df = pd.read_csv(r'D:\Projects\cellsandmachines-avishkar-5f1d4a1ac1d2\data\human\hg18_mrna_seq\sample',names=['mrna_gi_no','seq'])

refgene = pd.read_csv(refgene_file,sep='\t',header=None)#names=['bin','name','chrom']
#refGene 种有 ref名字相同的列，无法确定，去掉重复的列
refs = set(refgene[1])
drefs = {}
for x in refgene[1]:
    if drefs.has_key(x):
        drefs[x] = 1
    else:
        drefs[x] = 0
isdp = refgene[1].apply(lambda x: drefs[x])
refgene = refgene[isdp==0]    


def convert_NM2exonloc(nm):
    # > chr|+x,y
    rec = refgene[refgene[1]==nm]
    n = len(rec)
    if(n==0):
        return ''
    elif(n>1):
        return ''
    else:
        x = map(int,rec.iloc[0,9].rstrip(',').split(','))
        y = map(int,rec.iloc[0,10].rstrip(',').split(','))
        x.sort()
        y.sort()
        x = map(str,x)
        y = map(str,y)
        loc = rec.iloc[0,2] +'|' + rec.iloc[0,3] +','.join(x+y)
        return loc 

df = df[['mrna_gi_no']]
df['RNA_nucleotide_accession'] = df['mrna_gi_no'].apply(lambda x:x.split('.')[0]) #只保留accession不要版本号
df['locinfo'] = map(convert_NM2exonloc,df['RNA_nucleotide_accession'])
df = df[df['locinfo']!='']
df['chr'] = df['locinfo'].apply(lambda x:x.split('|')[0])
df['locations'] = df['locinfo'].apply(lambda x:x.split('|')[1])
#locations: 开始位置...,结束位置...

df = df[['mrna_gi_no','chr','locations']]
df.to_csv(os.path.join(out_folder,'mrna_exon_locations.txt'),index=False)
grouped = df.groupby(by=['chr'])
for chrom,group in grouped:
    group.to_csv(os.path.join(out_folder,'mrna_exon_locations\\'+str(chrom)+'_mrna_exon_locations.txt'),index=False)
