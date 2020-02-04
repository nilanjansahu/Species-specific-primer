from Bio import Entrez
from Bio import SeqIO
import os
import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from multiprocessing import Pool
import multiprocessing
from Bio.SeqUtils import MeltingTemp as mt

def list_of_rna(d):
    with open('16S', "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id.find(d)!=-1:
                #print(record.seq)
                return str(record.seq).upper()

def find_primer(choto,list_of_rna_to_be_compared):
    n=0
    for sequence in list_of_rna_to_be_compared:
        #print(sequence)
        if sequence.find(choto)!=-1:
            n+=1
            return None
    if n==0:
        return choto
def helper_find_primer(args):
    return find_primer(*args)
    
def possible_forward_seq(seq_d,n,seq_length,mt_min,mt_max,na, tris, mg, dntps,saltcor):
    a=seq_d[n:n+seq_length]
    if (a[-1]=='C' or a[-1]=='G') and len(a)==seq_length:
        if mt_min < mt.Tm_NN(a, Na=na, Tris=tris, Mg=mg, dNTPs=dntps,saltcorr=saltcor) and mt.Tm_NN(a, Na=na, Tris=tris, Mg=mg, dNTPs=dntps,saltcorr=saltcor) < mt_max:
            return str(a)
    
def helper_possible_forward_seq(args):
    return possible_forward_seq(*args)

def possible_reverse_seq(seq_d,n,seq_length,mt_min,mt_max,na, tris, mg, dntps,saltcor):
    a=seq_d[n:n+seq_length]
    if (a[0]=='C' or a[0]=='G') and len(a)==seq_length:
        if mt_min < mt.Tm_NN(a, Na=na, Tris=tris, Mg=mg, dNTPs=dntps,saltcorr=saltcor) and mt.Tm_NN(a, Na=na, Tris=tris, Mg=mg, dNTPs=dntps,saltcorr=saltcor) < mt_max:
            return str(a)
    
def helper_possible_reverse_seq(args):
    return possible_reverse_seq(*args)

def length_of_pcr_product(seq_d,possible_forward_primers,possible_reverse_primers,min_pcr_len,max_pcr_len):
    if seq_d.count(possible_forward_primers)==1 :
        if seq_d.count(possible_reverse_primers)==1 :
            a=seq_d.find(possible_reverse_primers)-seq_d.find(possible_forward_primers)+len(possible_reverse_primers)
            if min_pcr_len < a and a < max_pcr_len:
                return [possible_forward_primers,possible_reverse_primers]
    return None

def helper_length_of_pcr_product(args):
    return length_of_pcr_product(*args)


if __name__ == '__main__':
    possible_forward_prmr=[]
    possible_reverse_prmr=[]
    list_of_rna_to_be_compared=[]
    Entrez.email = "nilanjansahu@gmail.com"

    micro_oranism='Lentibacillus populi strain WD4L-1'
    #possible seq length of primers you want
    max_seq_length=24
    min_seq_length=20
    #from the list blastn
    num_of_16s_seq=100
    #pcr product length 
    min_pcr_len=190
    max_pcr_len=200
    #melting temperature
    mt_min=58
    mt_max=62
    #sodium concentration
    na=50
    #tris conc
    tris=10
    #Mg conc
    mg=1.5
    #dNTPs
    dntps=1
    #salt correction
    saltcor=7

    handle = Entrez.esearch(db="nucleotide", retmax=30, term=micro_oranism + ' 16s ribosomal RNA[Title]', idtype="acc")
    record = Entrez.read(handle)
    d = record['IdList'][0]
    seq_d=''
    records = SeqIO.parse(Entrez.efetch(db="nucleotide", id=str(d), rettype="gb"), "gb")
    for i,record in enumerate(records):
        for feature in record.features:    
            if feature.type == "rRNA": 
                if feature.location.extract(record).seq[0]=="A" or feature.location.extract(record).seq[0]=="T" or feature.location.extract(record).seq[0]=="G" or feature.location.extract(record).seq[0]=="C" or feature.location.extract(record).seq[0]=="U":
                    if len(feature.location.extract(record).seq)>1000:
                        writerr = open('seq_d.fas','w+')
                        print('>'+str(d)+' '+Entrez.read(Entrez.esummary(db="nucleotide", id=str(d)))[0]['Title'])
                        writerr.write('>'+str(d)+' '+Entrez.read(Entrez.esummary(db="nucleotide", id=str(d)))[0]['Title']+'\n')
                        seq_d=feature.location.extract(record).seq
                        print(seq_d)
                        writerr.write(str(seq_d))
                        writerr.close()
  
    num_processs=32
    p = Pool(num_processs)
    possible_forward_prmr=np.unique([x for x in p.map(helper_possible_forward_seq,[(seq_d,n,seq_length,mt_min,mt_max,na, tris, mg, dntps,saltcor) for n in range(len(seq_d)) for seq_length in range(min_seq_length,max_seq_length)]) if x is not None])
    possible_reverse_prmr=np.unique([x for x in p.map(helper_possible_reverse_seq,[(seq_d,n,seq_length,mt_min,mt_max,na, tris, mg, dntps,saltcor) for n in range(len(seq_d)) for seq_length in range(min_seq_length,max_seq_length)]) if x is not None])
    list_of_rna_to_be_compared=p.map(list_of_rna,pd.read_csv('seq_d.csv').iloc[:,1].values[1:num_of_16s_seq])
    possible_forward_primers = np.unique([x for x in p.map(helper_find_primer,[(x, list_of_rna_to_be_compared) for x in possible_forward_prmr]) if x is not None])
    possible_reverse_primers = np.unique([x for x in p.map(helper_find_primer,[(x, list_of_rna_to_be_compared) for x in possible_reverse_prmr]) if x is not None])
    primers=np.asarray([x for x in p.map(helper_length_of_pcr_product,[(seq_d,possible_forward_primer,possible_reverse_primer,min_pcr_len,max_pcr_len) for possible_forward_primer in possible_forward_primers for possible_reverse_primer in possible_reverse_primers]) if x is not None])
    p.close()
    print('forward primer')
    forward_primer=np.unique(primers[:,0])
    print(forward_primer)
    print('reverse primer')
    reverse_primer=[]
    for pr in np.unique(primers[:,1]):
        reverse_primer.append(str(Seq(pr, IUPAC.unambiguous_dna).complement()))
    print(np.unique(reverse_primer))

    
    
    









