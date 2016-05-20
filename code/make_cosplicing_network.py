__author__ = 'Delasa'
import pandas as pd
import os
import os.path as path
import numpy as np
from load_data import get_config

config = get_config()
data_path = config['data']['path']
output_path = config['outputs']['path']
run_mode = config['settings']['mode']


def find_cosplied_transcipts():
    coexpression_file = output_path + 'coexpression_with_SFs_%s.txt'%run_mode
    coexpression_df = pd.read_csv(coexpression_file, sep='\t')
    SFs = list(set(coexpression_df.SF_id.values))
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print('started making cosplicing netwoks ...')
    results_dfs = []
    for sf in SFs:
        meme_output_path = output_path + sf + '/'
        meme_results = os.listdir(meme_output_path)
        for meme_res in meme_results:
            if path.isdir(meme_output_path + meme_res):
                #print 'processing ',meme_res, ' started'
                status = None
                if meme_res.lstrip('meme_').find('R1') == 0:
                    status = 'R1'
                elif meme_res.lstrip('meme_').find('R2') == 0:
                    status = 'R2'
                elif meme_res.lstrip('meme_').find('R3') == 0:
                    status = 'R3'
                elif meme_res.lstrip('meme_').find('R4') == 0:
                    status = 'R4'
                else:
                    print('invalid meme output')
                    return False
                SF_name = meme_res.lstrip(status + '_')#region_strip(mo)
                this_df = find_enriched_transcripts(meme_output_path + meme_res, status, sf)
                if type(this_df) != bool:
                    results_dfs.append(this_df)
        print 'processing ',meme_res, ' finished'
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    final_df = pd.concat(results_dfs, axis=0)
    return final_df

def find_enriched_transcripts(path, status, SF_id):
    print (path, status, SF_id)
    test_file = open(path + '/meme.txt','r')
    test_lines = test_file.readlines()
    n = 0
    for i in range(len(test_lines)):
        line = test_lines[i]
        if line.find('SUMMARY OF MOTIFS') != -1:
            print i, line
            n = i
            break
    if n == 0:
        print(status, SF_id)
        return False

    filtered_lines = test_lines[n+8:len(test_lines)-12] #skip the last 13 lines
    dic = {}
    dic_exons = {}
    for thisTranscript in filtered_lines:
        resAr = thisTranscript.split()
        if len(resAr) > 1:
            trans_exon = resAr[0].strip()
            hyphen = trans_exon.find('_',6)
            tcons_id = trans_exon[:hyphen]
            exon_number = trans_exon[hyphen+1:]
            p_value = float(resAr[1].strip())

            if p_value < 0.05: #added this part to get the exons with significant motifs only
                if tcons_id in dic.keys():
                    dic[tcons_id].append(p_value)
                    dic_exons[tcons_id].append(exon_number)
                else:
                    dic[tcons_id] = [p_value]
                    dic_exons[tcons_id] = [exon_number]
    df = pd.DataFrame(data=None,  columns=['SF_id','target', 'interaction','avg_p_value','log_pvalue']) #log_pvalue can be used to weight the edges in the co-splicing network
    df.index.name = 'target_id'
    for transcript in dic.keys():
        p_value = np.mean(dic[transcript])
        if SF_id.strip() != transcript.strip():
            df.ix[SF_id.strip() + ':' + transcript.strip()] = SF_id.strip(),transcript.strip(),status, p_value, -1 * np.log2(p_value)
    test_file.close()
    return df

#cosplicing_df = find_cosplied_transcipts()
#cosplicing_df.to_csv(output_path + 'cosplicing_network_%s.txt'%run_mode, sep='\t', index=False)