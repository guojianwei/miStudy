# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 15:00:27 2016

@author: G2
数据生成：1、mmc1中生成mirna:target 的位置 mrna序列 具体见train_data_generate 注释
2、保守性分数
3、cds信息 genbank中cds信息
这段代码：
train_data_generate 函数从mmc1中生成用于训练集的mirna:target 的位置(1based)
依赖文件：gene2ensembl 从Ensembl(86)的rna号到 ncbi的 rna reqseq号RNA_nucleotide_accession.version
        mmc1
        Results-Homo_sapiens_Tools_IDMapper_.csv ensembl的不同ID版本历史转换60->86
        h_rna_108_38.fa

"""
import pandas as pd
import RNApar as myrnapar
def readascsvmmc1():
    mmc1 = "data/mmc1.txt"
    with open(mmc1,'r') as f:
        lines = f.readlines()
        #30行是表头，31行是第一行数据，18544是最后一行
        
        head = lines[30]
        lines = lines[31:]   
        head = head.split()
        data = []
        for line in lines:
            line = line.split()
            data.append(line)
        df = pd.DataFrame(data,columns=head)
        return df

def readasdictg2e():
    '''
    从Ensembl(86)的rna号到 ncbi的 rna reqseq号RNA_nucleotide_accession.version
    '''
    g2e = r'D:\data\webdata\gene2ensembl'
    df = pd.read_csv(g2e,sep='\t')
    df = df.drop_duplicates(['RNA_nucleotide_accession.version','Ensembl_rna_identifier']) #602133
    df.index = range(len(df))
    map_dict = {}
    for i,v in enumerate(df['Ensembl_rna_identifier']):
        map_dict[v] = df['RNA_nucleotide_accession.version'][i]
    return map_dict
        
def read_ID_His_Converter():
    '''
    只要86版的
    '''
    def get86(releaseno):
        s = releaseno.split(':')
        if(s[0] == '86'):
            return s[1][:-2]
        else:
            return ''
    id_convert =  pd.read_csv('D:\data\webdata\Results-Homo_sapiens_Tools_IDMapper_.csv') #是否应该所有的est都需要重新映射？
    id_convert['est86'] = id_convert['Releases'].apply(get86)
    id_convert = id_convert[id_convert['est86']!='']
    id_convert = id_convert[['Requested ID','est86']]
    return id_convert

def train_data_generate():
    '''
    测序数据mmc1使用ensemble 60和mirbase15，但ensemble 60已经找不到下载mrna,所以使用ncbi的ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/RNA/ GRCh38.p7 NCBI Homo sapiens Annotation Release 108
    1.ensembl的http://asia.ensembl.org/Homo_sapiens/Tools/IDMapper?db=core可以进行不同ID版本历史转换
    2.转换ID后,ncbi有从gi映射到enseble86的文件直接映射，从而找到NCBI的ACCESSION.version
    3.mmc1中的位置可能有误差，重新找到位置信息
    ？？？rna片段直接使用blast映射？
    # train_data_mirnaseq15_mrnaseqncbi_target_1based 产生中间文件
    '''
    e2g = readasdictg2e() # gene2ensembl 从Ensembl(86)的rna号到 ncbi的 rna reqseq号RNA_nucleotide_accession.version
    mmc1 = readascsvmmc1()
    mmc1['estold'] = mmc1['mRNA_name'].apply(lambda x: x.split('_')[1])
    id_convert = read_ID_His_Converter()    # ensemble 60 -> 86
    df = pd.merge(mmc1,id_convert,how='left',left_on='estold',right_on='Requested ID')
    df['estold'][df['Requested ID'].notnull()] = df[df['Requested ID'].notnull()]['est86']
    def get_mrnagi(no):
        try:
            return e2g[no]
        except:
            return ''
    mrna_gi_no = df['estold'].apply(get_mrnagi)
    df['mrna_gi_no'] = mrna_gi_no
    df = df[df['mrna_gi_no']!='']
    #df.to_csv('data/train_data_mirna_mrna_target.csv',index=False)
    df = df[['microRNA_name','miRNA_end','miRNA_seq','mRNA_start','mRNA_end_extended','mRNA_seq_extended','estold','mrna_gi_no']]
    #df = df[['seq_ID','microRNA_name','miRNA_end','miRNA_seq','mRNA_start','mRNA_end_extended','mRNA_seq_extended','estold','mrna_gi_no']]
    mrna_dict = myrnapar.read_mRNAs('D:/data/h_rna_108_38.fa')
    def get_mrnaseq(gi):
        try:
            return mrna_dict[gi]
        except:
            return ''
    df['mrna_seq_ncbi'] = df['mrna_gi_no'].apply(get_mrnaseq)
    df = df[df['mrna_seq_ncbi']!='']
    def simple_blast(seq,seg):
        """
        此处是否不严谨，是否有多段序列一样，是否应该是1-based
        """
        if seq.count(seg)>1: #解决有多段序列一样
            return ''
        index = seq.find(seg)   #seq[1800:1800+len(seg)]
        if(index<0):
            return ''
        else:
            return str(index+1)+'_'+str(index+len(seg))
            
    df['location'] = map(simple_blast,df['mrna_seq_ncbi'],df['mRNA_seq_extended'])
    df = df[df['location']!='']
    df['mRNA_start']=df['location'].apply(lambda x: x.split('_')[0])
    df['mRNA_end_extended']=df['location'].apply(lambda x: x.split('_')[1])
    df.to_csv('data/train_data_mirnaseq15_mrnaseqncbi_target_1based.csv',index=False)        
    #df.to_csv('data/train_data_mirnaseq15_mrnaseqncbi_target_1based_seqID.csv',index=False)    
    return df

def data_set_output():
    '''
    根据clash的mmc1数据经过mrna的从ensemble 60到86的映射使用86的测序数据，并且找到mrnaseqextend片段位置，所产生的记录
    依赖文件：clash_mmc1_to_dataset.py产生的train_data_mirnaseq15_mrnaseqncbi_target_1based
    产生：1.mirna 测序 mirbase15
        2. mrna 测序数据 ensemble86
        3. clash 位点，形式是mirna，mrna ，start，end(1based)
    '''
    clash_data = 'data/train_data_mirnaseq15_mrnaseqncbi_target_1based.csv'
    df = pd.read_csv(clash_data)
    df['microRNA_name'] = df['microRNA_name'].apply(lambda x:'hsa-'+x.split('_')[2])
    
    df_mrna_seq = df[['mrna_gi_no','mrna_seq_ncbi']]
    df_mrna_seq = df_mrna_seq.drop_duplicates()
    df_mrna_seq.to_csv('data/sample_mrna_seq.txt',index=False)
    del df_mrna_seq

    df_mirna_seq = df[['microRNA_name','miRNA_seq']]    
    df_mirna_seq = df_mirna_seq.drop_duplicates()
    df_mirna_seq.to_csv('data/sample_mirna_seq.txt',index=False)
    del df_mirna_seq

    df_clash_locations = df[['microRNA_name','mrna_gi_no','mRNA_start','mRNA_end_extended']]
    df_clash_locations = df_clash_locations.drop_duplicates() 
    df_clash_locations.to_csv('data/sample_clash_locations.txt',index=False)
    del df_clash_locations
    return df
    
if __name__ == '__main__':
   #mmc1 = train_data_generate()
    df = train_data_generate()
    #df = data_set_output()