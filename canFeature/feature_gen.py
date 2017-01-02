#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 15:47:48 2016

@author: g3
"""
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna,generic_dna
from rnafold import *
from data import RNAdata
import pandas as pd
#序列碱基含量，1,2？
def composite_seq( seq ):
    '''
        < seq
        > 各种含量占seq长度的比例
        对于AA，AA和AAA都算1个 AAAA算两个
    '''
    temp = Seq(seq,alphabet=generic_dna)
    temp = temp.transcribe().upper()
    n = float(len(seq))
    patterns = ['A','U','G','C']
    comp = [temp.count(x) for x in patterns]
    comp2mer = []
    for i in range(4):
        for j in range(i+1,4):
            comp2mer.append((comp[i]+comp[j])/n)
    return comp2mer
    
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

def compute_features(records,dmirnas,dmrnas):
    
    for each in records:
         can = each.split('\t')
         mi = dmirnas.get_rnaseq(can[0])
         m = dmrnas.get_rnaseq(can[1])
         start = int(can[3].split(',')[0].split('-')[0])
         end = int(can[3].split(',')[0].split('-')[1])
         composite = []
         composite += ( composite_seq(m[start-1:end]) )
         for offset in [5,10,15,20,25,30,46]:             
             composite += ( composite_seq(m[end:end+offset]) )
             start0 = start-offset
             if start0<0:
                 start0 = 0
             composite += ( composite_seq(m[start0:start]) )
         print composite
    
if __name__ == '__main__':
    mirna_seq = pd.read_csv('../data/sample_mirna_seq.txt')
    mrna_seq = pd.read_csv('../data/sample_mrna_seq.txt')
    dmirnas = RNAdata(mirna_seq)
    dmrnas = RNAdata(mrna_seq)
    f = open('../data/mi_m_siteloc_candica111.txt','rb')
    lines = f.readlines()
    f.close()
    compute_features(lines[:2],dmirnas,dmrnas)
    
    

