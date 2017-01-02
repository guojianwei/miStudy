# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 16:49:17 2016

@author: G2
"""

#from pyspark import SparkContext, SparkConf
import pandas as pd
from data import Clashdata,RNAdata
from models import RNA
from rnahybrid import *
from pentamers import pentamers_discover
import INN_HB
#序列碱基含量，1,2？
def composite_count( seg ):
    '''
        < seq
        > 各种碱基含量占seq长度的比例
        对于AA，AA和AAA都算1个 AAAA算两个
    '''
    seg = seg.upper().replace('T','U')
    n = float(len(seg))
    patterns = ['A','U','G','C','AA','UU','GG','CC',
                'AU','AG','AC','UG','UC','GC',
                'UA','GA','CA','GU','CU','CG']
    comp = [seg.count(x)/n for x in patterns]
    return comp
#line_index1 = 0
def gc_content( seg ):
 #   global line_index1
   # line_index1 += 1
    seg = seg.upper().replace('T','U')
    n = float(len(seg))
    gcn = seg.count('G')+seg.count('C')
    if n == 0:
        #print "sp_feature gc_content 0 error",line_index1,line_index1/5
        return 0
    return gcn/n

#line_index = -1
def position_pre(mrna):
#    global line_index
#    line_index += 1
    mrna = mrna.upper().replace('T','U')
    pospre = []
    positions = [0,1,9,10,11,12,13,14,15,16]
    patterns = ['A','G','C','U']
    for pat in patterns:
        for pos in positions:
         #   pospre.append(mrna[pos]==pat)
            try:
                pospre.append(mrna[pos]==pat)
            except IndexError:
                pospre.append(False)
    for pos in positions:
        try:
            pospre.append(mrna[pos]=='A' or mrna[pos]=='U')
        except IndexError:
            pospre.append(False)
    pospre = map(int, pospre)
    return pospre


# 计算deldelg需要用到delg:因为 deldelg = delg - mfe_open
def deldelg_seq(rna_fold, mrna_seq, delg, site_loc):
    output = rna_fold(mrna_seq[site_loc[0]:site_loc[1]])
    try:
        mfe_open = parse_rna_fold_output(output)
    except:
        print "Caught exception while processing mRNA seq:", mrna_seq
        raise
    deldelg = delg - mfe_open
    return deldelg

#read pentamers with significant enrichment
f = open(r'../data/pentamers.txt','r')
pentamers = f.readlines()
pentamers = [p.rstrip() for p in pentamers]
f.close()
inn_hb = INN_HB.INN_HB()
index2 = -1
def compute_features(record):
    global index2
    index2 += 1
    (Mi,M,label,seed_loc,sitestart,siteend,energe,seed,rnapos1,rnafeatures)= record
    features = [Mi._name,M._name]
    # 碱基含量
    loc = (int(sitestart),int(siteend))
    seg = M._seq[loc[0]-1:loc[1]]
    com_features = composite_count(seg)
    features += com_features
    # GC含量，上下游 #上下游是指？30,36
    upstreamloc1,upstreamloc2 = loc[0]-30-1 , loc[0]-46-1
    upstreamloc1 = upstreamloc1 if upstreamloc1>0 else 0
    upstreamloc2 = upstreamloc2 if upstreamloc2>0 else 0
    gc_feature = map(gc_content,[ M._seq[upstreamloc1: loc[0]-1],
                                M._seq[upstreamloc2: loc[0]-1],
                                seg,
                                M._seq[loc[1]:loc[1]+30],
                                M._seq[loc[1]:loc[1]+46] ])
    features += gc_feature
    # 特定位置碱基类型
    pos0 = rnapos1+1
    pos20 = pos0-20
    pos20 = pos20 if pos20>=1 else 1
    seg2 = M._seq[pos20-1:pos0][::-1]
    posp_features = position_pre(seg2)
    features += posp_features
    
    mrna_gi_no,strand,seq_len,mrna_cds_start,mrna_cds_end,mrna_conscores = rnafeatures
    #print rnafeatures
    isUTR, smallt200, bigt800 = 0,0,0
    if loc[1]>mrna_cds_end:
        isUTR = 1
    if abs(loc[1]-mrna_cds_end)<200 or abs(loc[1]-seq_len)<200:
        smallt200 = 1
    if abs(loc[1]-mrna_cds_end)>800 and abs(loc[1]-seq_len)>800:
        bigt800 == 1
    # UTR位置，长度，保守性分数
    features += [isUTR, smallt200, bigt800, seq_len-mrna_cds_end]
    seed_scores = mrna_conscores.split(',')[rnapos1-7-1:rnapos1-1]
   # print "11111111:",seed_scores
    try:
        seed_scores_sum = reduce(lambda x,y: x+y, map(float,seed_scores))
    except TypeError:
        print index2,seed_scores
    features.append(seed_scores_sum/len(seed_scores))
    #print seed_scores_sum
    # Free energy of seed sequence binding
    free_energy_seed = inn_hb.base_pairing(Mi._seq[1:8])
    features.append(free_energy_seed)
#    features.append(energe)
    '''
        Free energy of seed sequence binding 7
        Pentamer motif match 8
    '''
    features = map(str, features)
    return features
    
    
if __name__ == '__main__':
    
    mirna_seq = pd.read_csv('../data/sample_mirna_seq.txt')
    mrna_seq = pd.read_csv('../data/sample_mrna_seq.txt')
    mrna_features = pd.read_csv('../data/sample_mrna_features.txt')
    dmirnas = RNAdata(mirna_seq)
    dmrnas = RNAdata(mrna_seq)
    drnafeatures = dict(map(lambda x:(x[0],x),mrna_features.values.tolist()))
    df = pd.read_csv('../data/candidate_noverlap.txt')
    df.drop(labels='structure',axis=1,inplace=1)
    #df = df[49655:]
    #print df.loc[49655,]
    values = df.values.tolist()
    trainset = map(lambda x:(RNA(x[0],dmirnas.get_rnaseq(x[0])),
                                        RNA(x[1],dmrnas.get_rnaseq(x[1])),
                                        x[2], x[3].split('-'),x[4],x[5],x[6],x[7],x[8],
                                        drnafeatures[x[1]]),
                                values)
    del df,values
    res = map(compute_features,trainset)
    with open(r'../data/trainset_features.txt','w') as f:
        for line in res:
            f.write(','.join(line)+'\n')
        
    '''
    sc = SparkContext()
    train_set = sc.parallelize(values)
    train_set = train_set.map(lambda x:(RNA(x[0],dmirnas.get_rnaseq(x[0])),
                                        RNA(x[1],dmrnas.get_rnaseq(x[1])),
                                        x[2], x[3].split('-'),x[4],x[5],x[6],
                                        drnafeatures[x[1]]))
    train_set = train_set.map(compute_features)
    train_set.saveAsTextFile("./g3data/featured")
    '''                                  
    
