__author__ = 'Delasa'

import os
import yaml
import pandas as pd
#from Bio import SeqIO
#from Bio.Seq import Seq

def get_config():
   # print(os.getcwd())
    f = open('config.yaml')
    config = yaml.safe_load(f)
    f.close()
    return config

config = get_config()
data_path = config['data']['path']
SFs = config['data']['SFs']
expressions = config['data']['expressions']
output_path = config['outputs']['path']

input_type = config['settings']['input_data_type']
sep = config['settings']['sep']

def read_data(filename):
    print 'started reading %s' %filename
    if input_type == 'excel':
        df = pd.read_excel(filename)
    else:
        df = pd.read_csv(filename, sep=sep)
    print('%s loaded '%filename)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    return df

def get_expressed_SFs(SF_df):
    print('started finsing expressed splicing factors ...')
    expressions_df = read_data(expressions)

    SF_df['SF_gene'] = SF_df['gene_id'].apply(lambda x: x.upper())
    #print(SF_df.head())
    expressions_df['gene_id'] = expressions_df['gene_id'].apply(lambda x: x.upper())
    expressions_df.set_index('transcript_id')
    #print(expressions_df.shape)

    users_SF = pd.merge(SF_df, expressions_df, left_on='SF_gene', right_on='gene_id', how='inner')
    users_SF.to_csv(output_path + 'expressed_SFs.csv')

    print('expressed SFs: ', users_SF.shape)
    print ('%s expressed SFs are found and saved in the output folder' %users_SF.shape[0])
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    return users_SF

def load_data():
    expression_df = read_data(expressions)
    expression_df.set_index('transcript_id', inplace=True)
    expression_values_df = expression_df._get_numeric_data()

    SF_df = read_data(SFs)
    expressed_SFs = get_expressed_SFs(SF_df)
    SF_ids = expressed_SFs.transcript_id.values

    return SF_ids, expression_values_df