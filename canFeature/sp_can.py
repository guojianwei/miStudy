#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 20:28:51 2016

@author: g3
"""
from pyspark import SparkContext, SparkConf
import pandas as pd
from data import Clashdata,RNAdata
from models import RNA
from rnahybrid import *
#import random
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
        clip_locations也应该坐标从1开始？，rna_hybrid从1开始
    '''
    # TODO: Move this to config
   # output = rna_hybrid(mrna._seq, micro_rna._seq, energy_threshold=-15.0)
   # loc_attrs = parse_rna_hybrid_output(output)
    loc_attrs = []
    constraints = [(2, 8)] #只关注2-8
    for cons in constraints:
        output = rna_hybrid(mrna._seq, micro_rna._seq,
                            energy_threshold=-1.0,
                            helix_constraint=cons)
        loc_attrs += parse_rna_hybrid_output(output, enforce_seed_match=cons) # G3 改动
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
        else:
            label = 2
        candidate_set.append((label, loc_attr))
    return (micro_rna,mrna,candidate_set)
    
def process_generate(rna_hybrid, itr):
    """
        > [(@Mirna,@Mrna,siteloc,candicasite)]
        candicasite: lable locattr
    """
    print 'Begin'
    def do_gen(x):
        MiRNA,MRNA,clipl_locs = x[0],x[1],x[2]
        return get_labeled_target_locations( rna_hybrid,MiRNA,MRNA,clipl_locs )

    rdditr = itr
    rddres = rdditr.map(do_gen)
    def to_local(candicate):
        resofcans = []
        (Mi,M,cans) = candicate
        for can in cans:
            resofcans.append(','.join([Mi._name,M._name,str(can[0]),str(can[1])])+'\n')
        return resofcans
    rddres = rddres.map(to_local)
    rddres = rddres.collect()
    #rddres.saveAsTextFile("./g3data/f")
    with open(r'/home/hadoop/gthird_party/output82/f.txt','w') as f:
        f.write(','.join(['mi_name','m_name','lable','loc','sitestart','siteend','energe','seed','rnapos1','structure'])+'\n')
        for each in rddres:
            for site in each:
                f.write(site)
    
def rna_gen(milist,mlist,mysc):
    '''
    <milist,mlist
    >[[@miRNA,@mrna,locs]...]
    '''
    Milist = []
    Mlist = []
    for mi in milist:
        mi_seq = dmirnas.get_rnaseq(mi)
        if mi_seq == None:
            continue
        Milist.append( RNA(mi,mi_seq) )
    for m in mlist:
        m_seq = dmrnas.get_rnaseq(m)
        if m_seq == None:
                continue
        Mlist.append( RNA(m,m_seq) )
    print 'Ms',len(Mlist)
    print 'Mis',len(Milist)
    Mlist = mysc.parallelize(Mlist,480)
    n1 = Mlist.count()
    print 'Mlist',n1
    res = Mlist.flatMap(lambda MRNA:
        map(lambda MiRNA:(MiRNA,MRNA,dclash.find_loc(MiRNA._name,MRNA._name)),Milist))
    n2 = res.count()
    print n2
   # print '%s-%s'%(mi,m)
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
                
def canset_gen(mirnalist,mrnalist,mysc):
    """
        > [(@Mirna,@Mrna,siteloc,candicasite)]
    """
    print 'generating fullsets...'
    rh = RNAHybrid()
    fullset = rna_gen(mirnalist,mrnalist,mysc)
    print 'getting lables...'
    process_generate(rh,fullset) 

    

print 'fetching rnas ...'
f = open(r'../data/mirna_list.txt','r')
lines = f.readlines()
f.close()
mirnalist = map(lambda x:x.rstrip(),lines) 
f = open(r'../data/mrna_list.txt','r')
lines = f.readlines()
f.close()
mrnalist = map(lambda x:x.rstrip(),lines)
conf = SparkConf()
conf.setAppName('G3test')
sc = SparkContext(conf = conf)    
canset_gen(mirnalist,mrnalist,sc)
sc.stop()
#can_to_file(res)
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
