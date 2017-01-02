# -*- coding: utf-8 -*-
"""
Created on Sat Nov 05 09:32:32 2016

@author: G2
"""
import os
from Bio import SeqIO
"""
rna序列数据.fa.gbk
ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/RNA/对应ANNOTATION RELEASE NAME:	NCBI Homo sapiens Annotation Release 108
ASSEMBLY NAME:	GRCh38.p7
ASSEMBLY ACCESSION:	GCF_000001405.33
详细信息：ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/README_CURRENT_RELEASE
"""

def read_mRNAs(fname):
    #把给定文件的（h_rna.fa,）的NCBI的索引号accession NM
    records = {}
    for x in SeqIO.parse(fname, "fasta"):
        accession = x.name.split('|')[3]
        if records.has_key(accession):
            records[accession] = ''
        else:
            records[accession] = str(x.seq)
    return records
    
def read_miRNAs(fname):
    #把从miRBase中的（mature15.fa）转换为miRNA名字: miRNA序列的字典形式
    records = {}
    for x in SeqIO.parse(fname, "fasta"):
        records[x.id] = str(x.seq)
    return records



if __name__ == '__main__':
    #res = read_mRNAs(os.path.join(folder,'h_rna.fa'))
    #mi_dict = read_miRNAs(os.path.join(folder,'mature15.fa'))
    #msf = SeqIO.parse(os.path.join(folder,'h_rna.fa'), "fasta")
    pass