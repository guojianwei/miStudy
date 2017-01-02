# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 15:29:26 2016

@author: G2
"""
'''
conserv_score.py
1.mrna名字 -> 外显子坐标 -> 保守性分数
2.负链与正链反向互补
phastCons100way 保守性得分也是从ucsc数据库中下载所以也是0-based，与refgene中坐标一致！
所以只要refgene外显子坐标正确（hg38）phastCons100way（hg38）正确，结果一定正确，
mrna末端会有多个A并不在外显子坐标内所以mrnaseq会比cons多出几个（这些A没有分数）
'''
import pandas as pd
import os
def get_step(line):
    line = line.split()
    return line[3].split('=')[1]

def process_conserv(lines):
    def get_location(line):
        line = line.split() #line[1].split('=')[1],染色体组
        return line[2].split('=')[1]

    lines = lines
    index={}
    for i,v in enumerate(lines):
        if v[0]>'9': #不是数字，而是字母
            index[get_location(v)] = i+1
    return index
    
def find(locs,lines,index,fix_locations):
    def find_near(loc):
        fi = 0
        for i,v in enumerate(fix_locations):
            if loc<v:
                fi = i-1
                break
        return fix_locations[fi]
    start = int(locs[0])
    end = int(locs[1])
    near_fixlocation = find_near(start) #距离最近的fix location
    near_fixlocation2 = find_near(end) #距离最近的fix location
    if(near_fixlocation!=near_fixlocation2):#防止夸两段
        return ''
    line_i = index[str(near_fixlocation)] # 该 fix 在lines中的下标
    start = line_i+start-near_fixlocation    
    end = line_i+end-near_fixlocation
    for i in range(line_i,end): # 防止从中间有 fix location行
        if lines[i][0]>'9':
            return ''
    res = [x.rstrip('\n') for x in lines[start:end]]
    res = [str(float(x)) for x in res]
    return ','.join(res)

def lcos_to_scores(locations,lines,index,fix_locations):
    '''
    #locations: 开始位置...,结束位置...
    '''
    strand = locations[0]
    locations = locations[1:].split(',')
    n = len(locations)/2
    scores = []
    for i in range(n):
        score = find((locations[i],locations[i+n]),lines,index,fix_locations) 
        if score != '':
            scores.append(score)
        else:
            return ''
    return strand+'|'.join(scores)
       
def main(chrom):
    folder = r'data\prodata\hg38'
    print 'Reading and processing phastcons ...'
    f = open(r'D:\data\MicroRNA\conservationscores\hg38\chr'+str(chrom)+'.phastCons100way.wigFix','rb') #from ucsc
    lines = f.readlines()
    f.close()
    index = process_conserv(lines) # {染色体位置 : 文件位置}
    fix_locations = map(int,index.keys())
    fix_locations.sort() # 染色体起始位置

    print 'chr%s ... '%(str(chrom))
    dfloc = pd.read_csv(os.path.join(folder,'mrna_exon_locations\\chr'
                                        +str(chrom)+'_mrna_exon_locations.txt'))
    print 'Finding scores for mrnaexon in chr%s ... '%(str(chrom))
    dfloc['scores'] =[ lcos_to_scores(x,lines,index,fix_locations) for x in dfloc['locations']]    
    #dfloc = dfloc[dfloc['scores']!='']
    dfloc[['mrna_gi_no','scores']].to_csv(os.path.join(folder,'scores_mrna_exon\\chr'+str(chrom)+'_scores.txt'),index=False)
    del index,fix_locations,lines
    
    
        
# 1512143   1534687 13049808 13053659
# 12920889  12890130 12475704 12506463
# NM_199454
# NM_017818
'''
特殊情况：肯定横跨两段
开始239629072 239327323 215648222
结束239629539 239629518 215944382
start = 13232420    # NM_001010889 13232420 13198962 12674661 12708119
start = 134212702 # NM_001195025 mm9 30 vertebrate
'''

def scores_data_output():
    """
    把hg38的各个染色体的mrna分数合并成一个文件
    """
    folder = r'data\prodata\hg38\scores_mrna_exon'
    df = pd.DataFrame(columns=['mrna_gi_no','scores'])
    for item in os.listdir(folder):
        tempdf = pd.read_csv(os.path.join(folder, item))
        df = pd.concat([df,tempdf],ignore_index=True)
    df = df.dropna(axis=0)
    df = df[df['scores']!='']
    df['strand']=df['scores'].apply(lambda x: str(x)[0])
    def get_sores(signed_scores):
        x = str(signed_scores)
        strand = x[0]
        scores = x[1:].replace('|',',')
        if strand == '+':
            return scores
        elif strand == '-':
            scores = scores.split(',')
            scores.reverse()
            return ','.join(scores)
    df['mrna_conscores']=df['scores'].apply(get_sores)
    df = df[['mrna_gi_no','strand','mrna_conscores']]
    df.to_csv('data/sample_clash_mrnascores2_negreversed.txt',index=False)
    return df
if __name__ == '__main__':
    #stage 1
#    for i in range(1,23):
#        main(str(i))
#    main('X')
#    main('Y')
#    print "Scores fetched. Merging them..."
    #stage 2
    scores_data_output()
