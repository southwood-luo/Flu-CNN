import os
import argparse
import sys
import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc
import pickle
from Bio import SeqIO
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='UMAP',
                                 formatter_class=argparse.RawTextHelpFormatter)
args = parser.parse_args()


def plot_umap(adata, n_neighbors, arg, namespace='6interlayer_fc1in'):
    my_palette1 = {
        "Avian": "#B22222",
        "Human": "#6495ED",
    }
    n = str(n_neighbors)
    sc.pl.umap(adata, color='Host', size=20,# palette=my_palette1,
               save=f'figures/interlayer/{arg.seg}_{namespace}_{arg.data_type}_species{n}.png')
    sc.pl.umap(adata, color='HA', size=20,
               save=f'figures/interlayer/{arg.seg}_{namespace}_{arg.data_type}_ha{n}.png')


def load_dict(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def save_dict(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    segMaxLength = {'PB2': 2304, 'PB1': 2280, 'PA': 2163, 'HA': 1743,
                    'NP': 1515, 'NA': 1422, 'MP': 759, 'NS': 714}
    seg = 'NA'
    args.seg = seg
    args.mafft = 'mafft'
    args.data_type = 'protein'
    seg_max_length_dna = {'PB2': 2364, 'PB1': 2364, 'PA': 2229, 'HA': 2040,
                          'NP': 1554, 'NA': 1608, 'MP': 771, 'NS': 744}
    seg_max_length_protein = {'PB2': 798, 'PB1': 798, 'PA': 744, 'HA': 690,
                              'NP': 528, 'NA': 528, 'MP': 258, 'NS': 258}
    args.alphabet_path = 'alphabet_protein.json'
    args.num_features = 27
    args.alphabet = 'atcg'

    features_in_hook = load_dict(f'results/internal_layer/features_in_hook_fc1_{seg}')
    features_out_hook = load_dict(f'results/internal_layer/features_out_hook_fc1_{seg}')
    obs = load_dict('results/internal_layer/obs')

    ha_list = list()
    na_list = list()
    host_group_list = list()

    tt = list()
    for i in features_in_hook:
        p = i[0]
        tt.append(list(p.detach().numpy()[0]))
    adata = AnnData(pd.DataFrame(tt), dtype='float64')

    for key in ['Host','HA','NA']:
        adata.obs[key] = obs[key]

    n_neighbors = 10
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X')
    sc.tl.louvain(adata, resolution=0.2)
    sc.set_figure_params(dpi_save=500)
    sc.tl.umap(adata, min_dist=1)

    plot_umap(adata, n_neighbors, arg=args)
    print('++++++++++++++++++')
