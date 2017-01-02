# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 09:28:00 2016

@author: G2
"""
"""
处理spark生成的candidate提取特征？
"""
import os
import pandas as pd
def merge_hd(folder):
    '''
    简单的把每个文件合并成一个文件。
    '''
    all_records = []
    for item in os.listdir(folder):
        if os.path.isfile(os.path.join(folder, item)):
            with open( os.path.join(folder, item),'r' ) as f:
                lines = f.readlines()
                print item,len(lines)
                all_records += lines
    print 'Writing file...'
    with open (os.path.join(folder,r'merged\all_records.txt'),'w') as f:
        for records in all_records:
            f.write(records)
    
#merge_hd(r'D:\Projects\scripts\datasetgen\data\output')

"""
all_records格式
51600行43*1200，每行是一个mirna和一个mrna的所有候选位点
['hsa-let-7a\tNM_014431.2\t0\t35-62,-20.100,None,C  GG   UC   CUC  GUG      C: AC   GC  CCU   GC   GCCUCG : UG   UG  GGA   UG   UGGAGU :
U  AUA  UU        A', 'hsa-let-7a...']\n
"""
def save_ascsvhd(folder):
    f = open(os.path.join(folder,r'merged\all_records.txt'),'r')
    lines = f.readlines()
    f.close()
    print 'records...'
    data = []
    for records in lines:
        records = records.rstrip()[1:-1].split(',') #去掉每一行的换行和前后中括号
        n = len(records)
        i = 0
        while i < n:
            if(len(records[i])==0):
                break
            each = []
            each += records[i].lstrip()[1:].split('\\t') #去掉第一行的hsa前的空格和引号
            each.append(records[i+1]) #第二行是自由能
            each.append(records[i+2]) #位点类型
            each.append(records[i+3])
            each.append(records[i+4][:-1]) #双链结构，去掉最后的引号
            data.append(each)
            i += 5
    df = pd.DataFrame(data,columns=['mi_name','m_name','lable','loc','energe','seed','rnapos1','structure'])
    df.to_csv(os.path.join(folder,r'merged\candidate.txt'),index=False)

def delete_overlap(folder):
    '''
    我们真的需要去重吗？——不去重前的数据：
        hsa-let-7a,NM_014431.2,1,4469-4491,-27.2,1-8,4492
        hsa-let-7b,NM_014431.2,0,4469-4491,-27.2,1-8,4492
        hsa-let-7c,NM_014431.2,2,4469-4491,-27.2,1-8,4492
        hsa-let-7d,NM_014431.2,2,4471-4490,-24.3,2-8,4493
        hsa-let-7e,NM_014431.2,2,4469-4491,-29.4,1-9,4492
        hsa-let-7f,NM_014431.2,2,4469-4491,-21.6,1-8,4492
    从mirtarget数据库来看这样生成训练集是合适的。并没有实验验证的数据说明hsa-let-7b,NM_014431.2的结合
    说明直接确定类标是合理的。
    '''
    #microRNA_name,mrna_gi_no,mRNA_start,mRNA_end_extended
    df_clash = pd.read_csv(r'data/sample_clash_locations.txt')
    df_fulset = pd.read_csv(os.path.join(folder,r'f.txt'))
    dict_site_locs = {}        
    for x in df_clash.values.tolist():
        if dict_site_locs.has_key(x[1]):
            dict_site_locs[x[1]].append( (x[2],x[3]) )
        else:
            dict_site_locs[x[1]] = [ (x[2],x[3]) ]
    def get_overlap(m_name,label,loc):
        loc = map(int,loc.split('-'))
        if label == 1:
            return False
        else:
            try:
                cliplocs = dict_site_locs[m_name]
            except:
                cliplocs = []
            for site in cliplocs:
                if site[0]<=loc[0] and loc[1]<=site[1]:
                    return True
            return False
    
    df_fulset['overlap'] =map(get_overlap, df_fulset['m_name'],df_fulset['lable'],df_fulset['loc'])
    df_fulset = df_fulset[df_fulset['overlap']==False]
    df_fulset.drop(labels='overlap',axis=1,inplace=1)
    df_fulset.to_csv(os.path.join(folder,r'candidate_noverlap.txt'),index=False)
    
def static(folder):
    df = pd.read_csv(os.path.join(folder,r'candidate_noverlap.txt'))
   # df.drop(labels='structure',axis=1,inplace=1)
    return df

#merge_hd(folder)
#save_ascsvhd(folder)

folder = r'D:\Projects\scripts\datasetgen\data\output82'
#merge_hd(folder)
#save_ascsvhd(folder)
delete_overlap(folder)
df = static(folder)
#df.to_csv(r'data/sample_final_set.txt',index=False)
#df2 = df[df['lable']!=2]
#df2 = df[df['lable']!=2]
#df2.to_csv(r'D:\tmp\p28.txt',index=False)