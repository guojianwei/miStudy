# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 11:07:44 2016

@author: G2
miRDB:
The potential role of non-seed-based target binding was also
investigated. In this analysis, all pentamers from the non-seed region
of a miRNA were screened against the 50-end of the target-binding
site to identify any sequence match. In this way, significant enrichment
of pentamer-matching motifs were discovered in the target
sites (P ¼ 2.0E-15 with v2 test).
"""

def pentamers_discover(struct):
    '''
        在9:发现长度>=5的匹配不包括gu
        返回长度超过5的匹配模式列表
    '''
    mrna_mismatch, mrna_match,mirna_match, mirna_mismatch = struct
    max_align_len = [mrna_mismatch, mrna_match, mirna_mismatch, mirna_match]
    max_align_len = max(map(lambda x: len(x), max_align_len))
    mrna_mismatch += " " * (max_align_len - len(mrna_mismatch))
    mrna_match += " " * (max_align_len - len(mrna_match))
    mirna_mismatch += " " * (max_align_len - len(mirna_mismatch))
    mirna_match += " " * (max_align_len - len(mirna_match))
    #把这4个补到相同的长度（最长）
    pos = 0
    mipos = 0
    for i in range(1,max_align_len):
        pos += 1
        if not (mirna_mismatch[-i] == " " and mirna_match[-i] == " "):
            mipos += 1
        if mipos == 9:
            break
    start_pos = pos
    i = start_pos
    long_base_pair=[]
    while i<max_align_len+1:  
        start_pos = i
        end_pos = i
        for j in range(i,max_align_len+1):
            if ( mirna_match[-j] != " " and \
                not( (mirna_match[-j] == "G" and mrna_match[-j] == "U") or\
                (mirna_match[-j] == "U" and mrna_match[-j] == "G")) ):
                start_pos = j
                end_pos = j
                break
        for j in range(start_pos,max_align_len+1):
            if ( mirna_match[-j] == " " or \
                (mirna_match[-j] == "G" and mrna_match[-j] == "U") or\
                (mirna_match[-j] == "U" and mrna_match[-j] == "G") ):
                end_pos = j
                break
        if end_pos - start_pos>=5:
            long_base_pair.append(mirna_match[-end_pos+1:-start_pos+1])
        i = end_pos+1
    long_base_pair =[x[::-1] for x in long_base_pair]
    return long_base_pair
    
def chisquare(values):
    import scipy.stats
    '''
    <values are the 4 numbers in the table
    >return chisquare and p value
    '''
    a,b,c,d = values
    n = float(a+b+c+d)
    ab,cd,ac,bd = float(a+b),float(c+d),float(a+c),float(b+d)
    ea,eb,ec,ed = ab*ac/n,ab*bd/n, cd*ac/n,cd*bd/n
    return scipy.stats.chisquare([a,b,c,d],[ea,eb,ec,ed])
    
if __name__ == '__main__':
    '''
    a = '        UU     C        C'
    b = '   AAUAC  GCCUA CUACCUCA '
    c = '   UUAUG  UGGAU GAUGGAGU '
    d = 'UUG     U'
    print pentamers_discover([a,b,c,d])
    '''
    import pandas as pd
    d_pentamers = {}
    df = pd.read_csv(r'../data/candidate_noverlap.txt')
    npositive = len(df[df['lable']==1])
    nnegtive = len(df)-npositive
    for i in range(len(df)):
        label = 1 if df['lable'][i]==1 else 0
        structure = df['structure'][i]
        pentamers = pentamers_discover(structure.split(':'))
        for p in pentamers:
            for j in range(0,len(p)-4):
                if d_pentamers.has_key(p[j:j+5]):
                    d_pentamers[p[j:j+5]][label] += 1
                else:
                    d_pentamers[p[j:j+5]] = [0,0]
                    d_pentamers[p[j:j+5]][label] += 1
    
    pchi = []
    for key,value in d_pentamers.items():
        a = [key,value,chisquare([value[0],value[1],nnegtive-value[0],npositive-value[1]])]
        pchi.append(a)
    li = []
    for x in pchi:
        if x[2][1]<2.0e-15:
            li.append(x[0])
    print li
    #print chisquare([19,24,10,10])