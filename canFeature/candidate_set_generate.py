# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 14:26:56 2016

@author: G2
"""
import pandas as pd
from data import Clashdata,RNAdata
from models import RNA
from rnahybrid import *
import random
#def generate():

clash_loc = pd.read_csv('../data/sample_clash_locations.txt')
mirna_seq = pd.read_csv('../data/sample_mirna_seq.txt')
mrna_seq = pd.read_csv('../data/sample_mrna_seq.txt')
dclash = Clashdata(clash_loc)
dmirnas = RNAdata(mirna_seq)
dmrnas = RNAdata(mrna_seq)

'''
现根据hybrid生成所有candicate，再做特征提取
找到target的类型，找到site的起始位置，结合能
target结合能
'''
def get_labeled_target_locations(rna_hybrid,
                                 micro_rna, mrna, clipl_locations):
    '''
        Gets potential target locations in @mrna due to @micro_rna
        If the locations overlap/is contained in a clip location
        in clip_locations then
        return 1 for that location otherwise return 0 for that location
    '''
    # TODO: Move this to config
    output = rna_hybrid(mrna._seq, micro_rna._seq, energy_threshold=-15.0)
    loc_attrs = parse_rna_hybrid_output(output)
    constraints = [(1, 6), (2, 7), (3, 8)]
    for cons in constraints:
        output = rna_hybrid(mrna._seq, micro_rna._seq,
                            energy_threshold=-1.0,
                            helix_constraint=cons)
        loc_attrs += parse_rna_hybrid_output(output, enforce_seed_match=True)
    loc_attrs = resolve_ovarlapping_locations(loc_attrs)

    if loc_attrs is None:
        return None

    candidate_set = []
    for loc_attr in loc_attrs:
        label = 0
        if clipl_locations is not None:
            for other_loc in clipl_locations:
                if contains(loc_attr._loc, other_loc):
                    label = 1
                    break
        candidate_set.append((label, loc_attr))
    return (micro_rna,mrna,candidate_set)
    
def process_generate(rna_hybrid, itr):
    """
        > [(@Mirna,@Mrna,siteloc,candicasite)]
        candicasite: lable locattr
    """
    
    respg = []
    i = 0
    for (MiRNA,MRNA,clipl_locs) in itr:
        i += 1
        print i, MiRNA._name,MRNA._name
        respg.append( get_labeled_target_locations(
                                rna_hybrid,MiRNA,MRNA,clipl_locs))
    return respg
    
def rna_gen(milist,mlist):
    '''
    <milist,mlist
    >[[@miRNA,@mrna,locs]...]
    '''
    res = []
    for mi in milist:
        mi_seq = dmirnas.get_rnaseq(mi)
        MiRNA = RNA(mi,mi_seq)
        if mi_seq == None:
            continue
        for m in mlist:
            m_seq = dmrnas.get_rnaseq(m)
            if m_seq == None:
                continue
            MRNA = RNA(m,m_seq)
            clipl_locs = dclash.find_loc(mi,m)
            res.append((MiRNA,MRNA,clipl_locs))
            print '%s-%s'%(mi,m)
    return res
def can_to_file(candicateset):
    """
        <[(@Mirna,@Mrna,siteloc,candicasite)]
        >
    """
    tofile = '../data/mi_m_siteloc_candica111.txt'
    with open(tofile,'wb') as f:
        for (Mi,M,cans) in candicateset:
            for can in cans:
                f.write('\t'.join([Mi._name,M._name,str(can[0]),str(can[1])])+'\n')
                
def canset_gen(mirnalist,mrnalist):
    """
        > [(@Mirna,@Mrna,siteloc,candicasite)]
    """
    print 'generating fullsets...'
    rh = RNAHybrid()
    fullset = rna_gen(mirnalist,mrnalist)
    print 'getting lables...'
    res = process_generate(rh,fullset) 
    return res

    

print 'fetching rnas ...'
hsa_l7ab = clash_loc[:2]
dfmi = hsa_l7ab[['microRNA_name']].drop_duplicates()
dfm = hsa_l7ab[['mrna_gi_no']].drop_duplicates()
mirnalist = list(dfmi['microRNA_name'].values)
mrnalist = list(dfm['mrna_gi_no'].values)
res = canset_gen(mirnalist,mrnalist)
can_to_file(res)
'''
if __name__ == '__main__':   
    print 'fetching rnas ...'
    dfmi = clash_loc[['microRNA_name']].drop_duplicates()
    dfm = clash_loc[['mrna_gi_no']].drop_duplicates()
    mirnalist = list(dfmi['microRNA_name'].values)
    mrnalist = list(dfm['mrna_gi_no'].values)
    mirnalist = mirnalist[:30] 
    # TODO: 选择表达最丰富的10种miRNA家族,大约40-50个miRNA
    #mrnalist = random.sample(mrnalist,1000) #随机选择1000种mrna
    mrnalist = mrnalist[:100] #选择1000种mrna
    res = canset_gen(mirnalist,mrnalist)
    #can_to_file(res)
'''    
#    
#    temp_can = rna_gen(dfmi['microRNA_name'].values,dfm['mrna_gi_no'].values,mirnas,mrnas)
#    labled_can = []
#    for mi,miseq,m,mseq,loc in temp_can:
#        labled_can.append([mi,m,get_labeled_target_locations(rh,miseq,mseq,loc),loc])
'''
    df_can = pd.DataFrame(temp_can,columns=['microRNA_name','miseq','mrna_gi_no','mseq'])
    df_can = pd.merge(df_can,hsa_l7ab,how='left',on =['microRNA_name','mrna_gi_no'])
    df_can.drop(labels=['miseq','mseq'],axis=1).to_csv('../data/temp.csv',index=False)
'''