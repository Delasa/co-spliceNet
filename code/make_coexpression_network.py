__author__ = 'Delasa'
from scipy.stats import spearmanr,pearsonr
from load_data import get_config
from progressbar import ProgressBar
pbar = ProgressBar()

config = get_config()
correlation_method = config['settings']['correlation_method']
expressions = config['data']['expressions']
corr_cutoff = config['settings']['corr_cutoff']
run_mode = config['settings']['mode']
SFs = config['data']['SFs']
output_path = config['outputs']['path']

'''
expression_df = read_data(expressions)
expression_df.set_index('transcript_id', inplace=True)
expression_values_df = expression_df._get_numeric_data()

SF_df = read_data(SFs)
expressed_SFs = get_expressed_SFs(SF_df)
SF_ids = expressed_SFs.transcript_id.values
'''
#SF_ids, expression_values_df = load_data()

def perform_correlation(thisSF, expression_values_df, SF_ids):
    results = []
    thisSf_exp = list(expression_values_df.ix[thisSF].values)
    for thisTranscript in expression_values_df.index:
        thisgene_exp = list(expression_values_df.ix[thisTranscript].values)

        if thisTranscript not in SF_ids: #do not perform corr analysis on SFs
            if correlation_method == 'pearson':
                corr = pearsonr(thisSf_exp, thisgene_exp)
            else:
                corr = spearmanr(thisSf_exp, thisgene_exp)
        if corr[0] >= corr_cutoff:
            results.append(thisTranscript + '\t' + thisSF + '\t' + str(corr[0]) + '\t' + str(corr[1])+ '\n')
        #print(thisTranscript, thisSF, corr)
    return results

def make_coexpression_network(SF_ids, expression_values_df):
    outlines = ['transcript_id\tSF_id\tcorrelation_coefficient\tp-value\n']
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('started performing correlation analysis for %s splicing related proteins' %len(SF_ids))
    if run_mode == 'demo': #make co-expression and co-splicing network for one-SF only
        print('mode is demo. correlation analysis will be performed for one SF only.')
        thisSF = 'TCONS_00015312' #just do this for IRE1A
        thisSF_corrs = perform_correlation(thisSF, expression_values_df , SF_ids)
        outlines = outlines + thisSF_corrs
    else:
        print('mode is full. The coexpression network will be made for all the SFs present in the expression data')
        for thisSF in pbar(SF_ids):
            thisSF_corrs = perform_correlation(thisSF, expression_values_df , SF_ids)
            outlines = outlines + thisSF_corrs

    outfile = open(output_path + 'coexpression_with_SFs_%s.txt'%run_mode,'w')
    print('finished correlation analysis. coexpression_with_SFs.txt is saved in the output directory')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    outfile.writelines(outlines)
    outfile.close()

#make_coexpression_network()