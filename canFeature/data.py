# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 17:26:04 2016

@author: G2
"""

class Clashdata:
    def __init__(self,df_clash):
        temp = map(lambda x: (x[0]+x[1],[[x[2],x[3]]]) ,df_clash.values)
        self.clash_dict = {}        
        for name_loc in temp:
            if not self.clash_dict.has_key(name_loc[0]):
                self.clash_dict[name_loc[0]] = name_loc[1]
            else:
                self.clash_dict[name_loc[0]] += name_loc[1]
    def find_loc(self,mi,m):
        '''
        <mi,m
        >loc or None
        '''
        try:
            return self.clash_dict[mi+m]
        except:
            return None
    def is_clash(self,mi,m,loc):
        cloc = self.find_loc(mi+m)
        res = False
        if(cloc!=None): 
            for db_loc in cloc:
                if((db_loc[0]-20)<=loc[0] and loc[1]<=db_loc[1]):
                    res = True 
        else:
            res = False
        return res

class RNAdata:
    def __init__(self,df_rna_seq):
        self.rna_seqs = dict(map(lambda x: (x[0],x[1]) ,df_rna_seq.values))
        
    def get_rnaseq(self,rna_name):
        '''
            < rna name
            > rna seq or None
        '''
        try:
            return self.rna_seqs[rna_name]
        except:
            return None
    def get_allseqs(self):
        '''
            >all rna seqs
        '''
        return self.rna_seqs.items()
