__author__ = 'Delasa'

import pandas as pd
import numpy as np
import math
import time
import datetime
from Bio import SeqIO
from load_data import get_config
from Bio.Seq import Seq
from progressbar import ProgressBar

import os

config = get_config()
R1_length = config['settings']['R1_length']
R2_length = config['settings']['R2_length']
R3_length = config['settings']['R3_length']
R4_length = config['settings']['R4_length']

genome_seq = config['data']['genome_seq']
gtf = config['data']['gtf']
output_path = config['outputs']['path']

def print_timestamp():
    import time
    import datetime
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print(st)

def get_transcript_id(string):
   try:
        strAr = string.split(';')
        if strAr[1].find('transcript_id') != -1:
            transcript_id = strAr[1].strip('transcript_id ').strip('"')
            return transcript_id
   except KeyError:
       transcript_id = ''
       return transcript_id



def get_exon_id(string):
    strAr = string.split(';')
    exon_number = strAr[2].strip('exon_number ').strip('"')
    return exon_number


def get_oId(string):
    strAr = string.split(';')
    oId = strAr[4].strip('oId ').strip('"')
    return oId


def get_geneName(string):
    strAr = string.split(';')
    gene_name = strAr[3].strip('gene_name ').strip('"')
    return gene_name



def seqfinderold(start_pos, end_pos, fa_lines):
        start_line = int(math.floor(start_pos/60)) + 1
        end_line = int(math.floor(end_pos/60)) + 1

        rem1 = start_pos%60
        rem2 = end_pos%60
        if start_line == end_line:
                bps = fa_lines[end_line][rem1-1: rem2]
        else:
                bps = ''
                bps_first = fa_lines[start_line][rem1 - 1:]
                bps = bps_first.strip('\n')
                bps_end = fa_lines[end_line][:rem2].strip('\n')
                for i in range(start_line + 1, end_line):
                        line = fa_lines[i]
                        bps = bps + line.strip('\n')
                bps = bps + bps_end
        return bps

def preprocess_gtf():
    print('started preprocessing gtf file')
    merged_df = pd.read_csv(gtf, sep='\t', header=None)
    transcripts = merged_df[8].copy()

    transcript_ids = transcripts.apply(get_transcript_id)
    exons_ids = transcripts.apply(get_exon_id)
    oId = transcripts.apply(get_oId)
    geneName = transcripts.apply(get_geneName)

    merged_df['transcript_id'] = transcript_ids
    merged_df['exon_number'] = exons_ids
    merged_df['oId'] = oId
    merged_df['geneName'] = geneName

    merged_df.rename(columns={0:'chr', 1:'source',2:'type',3:'start',4:'end',5:'dot',6:'strand',7:'dot2',8:'properties'},inplace=True)
    merged_df.set_index(['transcript_id','exon_number'], inplace=True)
    print('processing gtf file is complete')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    return merged_df


def open_genome():
    handle = open(genome_seq, "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    #print record_dict.keys()
    return record_dict
######### seq finder
def seqfinder(chr, start_pos, end_pos,strand, genome):
    #genome = open_genome()
    seqs = genome[str(chr)]

    extracted_seq = seqs[start_pos:end_pos]
    if strand == '-':
        return extracted_seq.reverse_complement()
    else:
        return extracted_seq

def get_seqs(status='R1', tcons_id='TCONS_00000001', exons_number='1', gtf_processed='', genome=''):
    chr = gtf_processed.loc[tcons_id,exons_number]['chr']
    start = gtf_processed.loc[tcons_id,exons_number]['start']
    end = gtf_processed.loc[tcons_id,exons_number]['end']
    strand = gtf_processed.loc[tcons_id,exons_number]['strand']
    exon_length = abs(end - start + 1)

    seqs = genome[str(chr)]

    if status == 'R1':
        if exon_length < R1_length:
            return False
        else:
            seq = seqfinder(chr, start - R1_length , start,strand, genome)
    elif status == 'R2':
        if exon_length < R2_length:
            return False
        else:
            seq = seqfinder(chr, start , start + R2_length,strand, genome)
    elif status == 'R3':
        seq = seqfinder(chr, end - R3_length + 1, end + 1,strand, genome)
    elif status == 'R4':
        seq = seqfinder(chr, end + 1, end + R4_length + 1,strand, genome)
    else:
        print 'invalid R'
        return False
    seq.id = tcons_id + '_' +exons_number
    return seq


def extract_and_write_df(status, df, filename='output', SF_name='SF', genome=''):
    thislist = []
    print('started extracting sequences for %s in %s' %(SF_name,status))
    #gtf_processed = preprocess_gtf()
    pbar = ProgressBar()
    for tcons,exon in pbar(df.index):
       # print(tcons, exon)
        seq = get_seqs(status=status,tcons_id=tcons, exons_number=exon, gtf_processed=df, genome=genome)
        if seq != False:
            thislist.append(seq)
    if SF_name not in os.listdir(output_path):
        os.makedirs(output_path + SF_name)

    R_handle = open(output_path + SF_name + '/' + status + '_' + filename + ".fasta", "w")
    SeqIO.write(thislist, R_handle, "fasta")
    R_handle.close()
    print '%s is done for %s' %(status, SF_name)

















