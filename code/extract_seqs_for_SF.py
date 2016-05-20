__author__ = 'Delasa'
from extract_seqs import extract_and_write_df, preprocess_gtf, get_config, open_genome
import pandas as pd

config = get_config()
data_path = config['data']['path']
output_path = config['outputs']['path']
run_mode = config['settings']['mode']



def extract_seqs_for_SFs():
    gtf_processed = preprocess_gtf()
    gtf_processed.reset_index(inplace=True)
    gtf_processed.set_index('transcript_id', inplace=True)

    coexpression_file = output_path + 'coexpression_with_SFs_%s.txt'%run_mode
    coexpression_df = pd.read_csv(coexpression_file, sep='\t')
    SFs = set(coexpression_df.SF_id.values)

    for f in SFs:
        print '\nstarted processing %s' %f
        selected_sf_FNs = coexpression_df[coexpression_df.SF_id == f]
        selected_sf_FNs.set_index('transcript_id', inplace=True)
        concat = pd.merge(gtf_processed,selected_sf_FNs, left_index=True, right_index=True , how='inner')

        concat.index.name = 'transcript_id'
        concat.reset_index(inplace=True)
        concat.set_index(['transcript_id', 'exon_number'], inplace=True)

        genome = open_genome()
        extract_and_write_df(status='R1', df=concat, filename=f.strip('.csv'), SF_name=f, genome=genome)
        extract_and_write_df(status='R2', df=concat, filename=f.strip('.csv'), SF_name=f, genome=genome)
        extract_and_write_df(status='R3', df=concat, filename=f.strip('.csv'), SF_name=f, genome=genome)
        extract_and_write_df(status='R4', df=concat, filename=f.strip('.csv'), SF_name=f, genome=genome)

        print ('\nfinished processing file %s' %f)
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#extract_seqs_for_SFs()