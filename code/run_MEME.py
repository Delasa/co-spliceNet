__author__ = 'Delasa'
import os

from load_data import get_config
import pandas as pd

def run_meme():
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('started running MEME ...')
    config = get_config()
    output_path = config['outputs']['path']
    run_mode = config['settings']['mode']

    meme_mode = config['settings']['meme_mode']
    motif_min_width = config['settings']['motif_min_width']
    motif_max_width = config['settings']['motif_max_width']
    motif_evalue = config['settings']['motif_evalue']
    number_of_motifs = config['settings']['number_of_motifs']

    if motif_evalue == 'default':
        evalue = ''
    else:
        evalue = '-evt %s' %motif_evalue

    coexpression_file = output_path + 'coexpression_with_SFs_%s.txt'%run_mode
    coexpression_df = pd.read_csv(coexpression_file, sep='\t')
    SFs = list(set(coexpression_df.SF_id.values))

    scripts_outlines = []
    Rs = ['R1','R2','R3','R4']

    for sf in SFs:
        for R in Rs:
            dir_name = 'meme_' + R + '_' + sf
            input_seq_path = output_path + sf + '/' + R + '_' + sf + '.fasta'
            output_meme_path = output_path + sf + '/' + dir_name
            #if dir_name in os.listdir(output_path + sf):
            #    print('%s is already processed' %sf)
            #else:
            print('~~~~~~~~~~')
            cmd = 'meme %s -o %s -dna -minw %s -maxw %s -maxsize 200000 -mod %s -nmotifs %s  %s ' %(input_seq_path,output_meme_path,motif_min_width,motif_max_width,meme_mode,number_of_motifs,evalue)
            scripts_outlines.append(cmd + '\n')
            print cmd
            os.system(cmd)

        print('processing %s is done ' %sf)
        print('---------------------------------\n')
    scripts_file = open(output_path + 'meme_script.sh','w')
    scripts_file.writelines(scripts_outlines)
    scripts_file.close()

#run_meme()