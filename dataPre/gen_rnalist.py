# -*- coding: utf-8 -*-
"""
算法输入需要的miRNA list 和 mrna list
miRNA选择丰富表达的家族的miRNA ps miRNA表达最丰富的家族是<<MicroRNA target prediction using
    thermodynamic and sequence curves>>中实验数据的最丰富，但本实验使用的CLASH数据集肯定不同，
    所以最丰富的种类的确定应以本实验数据集为准
mrna选择特征比较全的mrna：+链，有scores的mrna
"""
import pandas as pd

def getrnalist():
    """
    只有features的mrna才可能加入训练集
    """
    df = pd.read_csv(r'data/sample_mrna_features.txt')
    #df = df[df['strand']=='+'] #+链有确定的scores，目前不确定-链的scores是否正确
    df = df[['mrna_gi_no']] #[:1200] #随机选择1200种
    df = df.drop_duplicates()
    df.to_csv(r'data/mrna_list.txt',header=None,index=False)

def getmilist_convert():
    """
    <mirnas0_sum_ver21_seq_family10
    >mirna list to convert
    mmc1 mirna ver15
    """
    df = pd.read_csv(r'mirnas0_sum_ver21_seq_family10.csv') 
    df = df[['Mirnas/Samples']]
    df.to_csv(r'data/minamewillto15.txt',header=None,index=False)

def getmilist():
    #df = pd.read_csv(r'data/Conversion Summary Reportmi10most-15.csv')
    #df = df[['VER_15']]
    #df = df.dropna()
    dfm = pd.read_csv(r'data/sample_mrna_features.txt')
    dfm = list(dfm['mrna_gi_no'].values)
    dfloc = pd.read_csv(r'data/sample_clash_locations.txt')
    dfloc['flag'] = dfloc['mrna_gi_no'].apply(lambda x: True if str(x) in dfm else False)
    dfloc = dfloc[dfloc['flag']==True]
    df = dfloc[['microRNA_name']]
    df = df.drop_duplicates()
    df.to_csv(r'data/mirna_list.txt',header=None,index=False)
    
if __name__ == '__main__':
    getrnalist()
    getmilist()