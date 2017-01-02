# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 20:25:18 2016

@author: G2
"""
import pandas as pd
from Bio import SeqIO
import os

def read_mRNAs_forcds(fname):
    #把给定文件的（genban,）的NCBI的索引号gi NM
    #rna序列数据.fa.gbk

    records = {}
    for x in SeqIO.parse(fname, "genbank"):
        accession = x.id
        if records.has_key(accession):
            records[accession] = ''
        else:
            records[accession] = x

    return records
#msg = read_mRNAs_forcds(os.path.join(folder,'h_rna.gbk'))
'''
并不能确定features中那个是CDS信息
genbank中cds遵循 0 base start 1 base end 例如[0:529] 长529,第一个下标0，最后一个下标528
cds也是同理
'''
'''
i = 0
for x in SeqIO.parse(os.path.join(folder,'h_rna.gbk'), "genbank"):
    try: 
        f=x.features
    except :
        continue
    try:
        print '---'
        print f[0].type
        print f[1].type
        print f[2].type
        print f[3].type
        print f[4].type
        print f[5].type
        i += 1
    except:
        pass
    if i>100:
        break
'''
#x.features[cds].location.start.position  .location.end.position
def get_cds_mrna():
    """
        根据refGene.txt计算出sample_mrna_seq.txt中mrna的所有外显子在ucsc中chrome的坐标，
        输出到_mrna_exon_locations文件中
    """
    clash_data = r'data\train_data_mirnaseq15_mrnaseqncbi_target_1based.csv'
    df = pd.read_csv(clash_data)
    df_mrna_cds = df[['mrna_gi_no']]
    del df
    df_mrna_cds = df_mrna_cds.drop_duplicates()
    folder = 'D:\data'
    dmg = read_mRNAs_forcds(os.path.join(folder,'h_rna_108_38.gbk')) 
    def getcds(mrna_acc):
        try:
            x = dmg[mrna_acc]
        except:
            return ''
        try:
            f = x.features
        except:
            return ''
        for fi in f:
            if fi.type == 'CDS':
                return str(len(x.seq))+'_'+str(fi.location.start.position)+'_'+str(fi.location.end.position)
        return ''
    df_mrna_cds['mrna_cds'] = df_mrna_cds['mrna_gi_no'].apply(getcds)
    df_mrna_cds = df_mrna_cds[df_mrna_cds['mrna_cds']!='']
    #df_mrna_cds.to_csv(r'data\sample_mrna_cdsloc.txt',index=False)
    with open(r'data\sample_mrna_cdsloc.txt','w') as f:
        f.write('mrna_gi_no,seq_len,mrna_cds_start,mrna_cds_end\n')
        for i,v in enumerate(df_mrna_cds.values):
            line = ','.join([v[0]]+v[1].split('_'))+'\n'
            f.write(line)

get_cds_mrna()
df1 = pd.read_csv(r'data\sample_clash_mrnascores2_negreversed.txt')
df2 = pd.read_csv(r'data\sample_mrna_cdsloc.txt')
df = pd.merge(df1,df2,how='inner',left_on ='mrna_gi_no',right_on='mrna_gi_no')
df = df[['mrna_gi_no', 'strand', 'seq_len',
             'mrna_cds_start', 'mrna_cds_end','mrna_conscores']]
df[['seq_len','mrna_cds_start', 'mrna_cds_end']] = df[['seq_len','mrna_cds_start', 'mrna_cds_end']].astype(int)
df.to_csv(r'data\sample_mrna_features.txt',index=False)