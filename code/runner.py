__author__ = 'Delasa'

from load_data import get_config, load_data
from make_coexpression_network import make_coexpression_network
from extract_seqs_for_SF import extract_seqs_for_SFs
from run_MEME import run_meme
from make_cosplicing_network import find_cosplied_transcipts


config = get_config()

data_path = config['data']['path']
SFs = config['data']['SFs']
output_path = config['outputs']['path']
input_type = config['settings']['input_data_type']
sep = config['settings']['sep']
run_mode = config['settings']['mode']

correlation_method = config['settings']['correlation_method']
expressions = config['data']['expressions']
corr_cutoff = config['settings']['corr_cutoff']

if __name__ == '__main__':
    print('started brewing a co-splicing network fory you ...\n')
    SF_ids, expression_values_df = load_data() #loading data

    make_coexpression_network(SF_ids=SF_ids, expression_values_df=expression_values_df) #making coexpression networks

    extract_seqs_for_SFs() #extracting R1, R2, R3, and R4 sequences for running denovo motif discovery

    run_meme() #meme will be run automatically. but the commands will be saved in the meme_script.sh in the output directory for reproducing the results

    cosplicing_df = find_cosplied_transcipts() #finding transcripts that are co-expressed with a SF and have at least one signfiicant motif in their splice juctions
    cosplicing_df.to_csv(output_path + 'cosplicing_network_%s.txt'%run_mode, sep='\t', index=False)
    print('cosplicing network is ready.')
    print('########################## DONE :)  ######################################')


