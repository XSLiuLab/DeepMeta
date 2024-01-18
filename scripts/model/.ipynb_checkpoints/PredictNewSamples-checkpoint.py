# Predict TCGA (or other new) samples using a trained model
#print("\n\nStarting to run PredictNewSamples.py with a demo example of 10 TCGA tumors...")

import numpy as np
import pandas as pd
from keras import models
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--exp", help="expression")
parser.add_argument("-f", "--fingerprint", help="gene fingerprint feature")
parser.add_argument("-o", "--output", help="prediction results")
args = parser.parse_args()
exp_file = args.exp
f_file = args.fingerprint
save_file = args.output

def load_data(filename):
    data = []
    gene_names = []
    data_labels = []
    lines = open(filename).readlines()
    sample_names = lines[0].replace('\n', '').split('\t')[1:]
    dx = 1

    for line in lines[dx:]:
        values = line.replace('\n', '').split('\t')
        gene = str.upper(values[0])
        gene_names.append(gene)
        data.append(values[1:])
    data = np.array(data, dtype='float32')
    data = np.transpose(data)

    return data, data_labels, sample_names, gene_names

if __name__ == '__main__':
    model_saved = models.load_model("/home/data/sda/wt/DeepDEP/model_exp_paper.h5")
    # model_paper is the full 4-omics DeepDEP model used in the paper
    # user can choose from single-omics, 2-omics, or full DeepDEP models from the
    # /data/full_results_models_paper/models/ directory
    
    # load TCGA genomics data and gene fingerprints
    # data_mut_tcga, data_labels_mut_tcga, sample_names_mut_tcga, gene_names_mut_tcga = load_data(
    #     "/data/tcga_mut_data_paired_with_ccl.txt")
    data_exp_tcga, data_labels_exp_tcga, sample_names_exp_tcga, gene_names_exp_tcga = load_data(exp_file)
    # data_cna_tcga, data_labels_cna_tcga, sample_names_cna_tcga, gene_names_cna_tcga = load_data(
    #     "/data/tcga_cna_data_paired_with_ccl.txt")
    # data_meth_tcga, data_labels_meth_tcga, sample_names_meth_tcga, gene_names_meth_tcga = load_data(
    #     "/data/tcga_meth_data_paired_with_ccl.txt")
    data_fprint_1298DepOIs, data_labels_fprint, gene_names_fprint, function_names_fprint = load_data(f_file)
    print("\n\nDatasets successfully loaded.\n\n")

    batch_size = 500
    first_to_predict = data_exp_tcga.shape[0]
    # predict the first 10 samples for DEMO ONLY, for all samples please substitute 10 by data_mut_tcga.shape[0]
    # prediction results of all 8238 TCGA samples can be found in /data/full_results_models_paper/predictions/

    t = time.time()
    data_pred = np.zeros((first_to_predict, data_fprint_1298DepOIs.shape[0]))
    for z in np.arange(0, first_to_predict):
        data_pred_tmp = model_saved.predict([data_exp_tcga[np.repeat(z, data_fprint_1298DepOIs.shape[0])],
                                             data_fprint_1298DepOIs], batch_size=batch_size, verbose=0)
        data_pred[z] = np.transpose(data_pred_tmp)
        print("TCGA sample %d predicted..." % z)
        
    # write prediction results to txt
    data_pred_df = pd.DataFrame(data=np.transpose(data_pred), index=gene_names_fprint,
                                columns=sample_names_exp_tcga[0:first_to_predict])
    pd.DataFrame.to_csv(data_pred_df, path_or_buf=save_file, sep='\t', index_label='CRISPR_GENE', float_format='%.4f')
    
