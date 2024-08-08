import os
import argparse
import numpy as np
from model import CharCNN
import torch
from torch.autograd import Variable
import torch.nn.functional as F
import pandas as pd
from anndata import AnnData
import scanpy as sc
import matplotlib.pyplot as plt
import pickle
import torch.nn as nn
import math
import copy
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Character level CNN text classifier testing',
                                 formatter_class=argparse.RawTextHelpFormatter)
# model
parser.add_argument('--model_path', default=None,
                    help='Path to pre-trained acouctics model created by DeepSpeech training')
parser.add_argument('--dropout', type=float, default=0.5, help='the probability for dropout [default: 0.5]')
parser.add_argument('--l0', type=int, default=1014, help='maximum length of input sequence to CNNs [default: 1014]')
parser.add_argument('--kernel-num', type=int, default=100, help='number of each kind of kernel')
parser.add_argument('--kernel-sizes', type=str, default='3,4,5',
                    help='comma-separated kernel size to use for convolution')
# data
parser.add_argument('--test-path', metavar='DIR',
                    help='path to testing data csv', default='data/ag_news_csv/test.csv')
parser.add_argument('--batch-size', type=int, default=20, help='batch size for testing [default: 128]')
parser.add_argument('--alphabet-path', default='alphabet.json', help='Contains all characters for prediction')
# device
parser.add_argument('--num-workers', default=4, type=int, help='Number of workers used in data-loading')
parser.add_argument('--cuda', action='store_true', default=False, help='enable the gpu')
# logging options
parser.add_argument('--save-folder', default='Results/', help='Location to save epoch models')
parser.add_argument('--class_type', type=str, default='host',
                    help='which type data to analysis [default:host]')
parser.add_argument('--data_type', type=str, default='dna',
                    help='which type data to analysis [default:dna]')
args = parser.parse_args()


def char2Index(character, arg):
    alphabet = arg.alphabet
    return alphabet.find(character)


def oneHotEncode(sequence, arg):
    X = torch.zeros(arg.num_features, arg.l0)
    for index_char, char in enumerate(sequence[::-1]):
        if char2Index(char, arg) != -1:
            X[char2Index(char, arg)][index_char] = 1.0
    return X


def seqProcess(seq, arg):
    inputs = oneHotEncode(seq, arg)
    inputs = inputs.unsqueeze(0)
    return inputs


def load_dict(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def predictor(inp):
    inp = Variable(inp)
    logit = model(inp)
    # log_softmax = nn.LogSoftmax(dim=1)
    probs = torch.sigmoid(logit)[0]
    # probs = F.softmax(logit, dim=1)[0]
    return probs, logit


if __name__ == '__main__':
    segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    dict_reassort_prob = dict()
    # segs = ['HA']
    for seg in segs:
        args.seg = seg
        args.mafft = 'unmafft'
        args.data_type = 'protein'
        print('\n开始测试：' + seg)
        if args.data_type == 'dna':
            seg_max_length_dna = {'PB2': 2364, 'PB1': 2364, 'PA': 2229, 'HA': 2040,
                                  'NP': 1554, 'NA': 1608, 'MP': 771, 'NS': 744}
            args.l0 = seg_max_length_dna[seg]
            args.alphabet_path = 'alphabet_dna.json'
            args.num_features = 5
            args.alphabet = 'atcg'
        elif args.data_type == 'protein':
            seg_max_length_protein = {'PB2': 798, 'PB1': 798, 'PA': 771, 'HA': 609,
                                      'NP': 528, 'NA': 501, 'MP': 285, 'NS': 258, 'all8': 4443}
            # seg_max_length_protein = {'PB2': 1041, 'PB1': 906, 'PA': 1068, 'HA': 852,
            #                           'NP': 582, 'NA': 717, 'MP': 582, 'NS': 339}
            args.l0 = seg_max_length_protein[seg]
            args.alphabet_path = 'alphabet_protein.json'
            args.alphabet = 'abcdefghijklmnopqrstuvwxyz'
            args.num_features = 27

        if args.class_type == 'host':
            args.num_class = 2
            args.path = f'data/other_try/zoonotic_strain.csv'
            if args.data_type == 'protein':
                # args.model_path = f'models_CharCNN/3parts/CharCNN_host_protein_{seg}_{args.mafft}_best_nll_3parts.pth.tar'
                args.model_path = f'models_CharCNN/3parts_smallsub/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll_3parts_smallsub.pth.tar'
                # args.model_path = f'models_CharCNN/111and114_unmafft/CharCNN_host_protein_{seg}_{args.mafft}_best_nll.pth.tar'
                # args.model_path = f'models_CharCNN/CharCNN_host_protein_{seg}_{args.mafft}_best_nll.pth.tar'

        # load testing data
        print("\nLoading testing data...")
        print(args.path)
        model = CharCNN(args)
        print("=> loading weights from '{}'".format(args.model_path))
        assert os.path.isfile(args.model_path), "=> no checkpoint found at '{}'".format(args.model_path)
        checkpoint = torch.load(args.model_path)
        model.load_state_dict(checkpoint['state_dict'])
        model.eval()
        dict_genomeID_seq = load_dict('dict/dict_genomeID_seq')
        dict_genomeID_seg_seq = load_dict('dict/dict_genomeID_seg_seq')
        dict_name_genomeID = load_dict('dict/dict_name_genomeID')
        dict_name_seg_seq_lack = load_dict('data/other_try/dict_name_seg_seq_lack')

        fread = open(args.path, 'r')
        for rline in fread:
            lrline = rline.rstrip().split(',')
            name = lrline[1]
            strain_name = name.replace(' ', '_').replace('-', '_').lower()
            strain_name = strain_name[::-1].replace(' ', '', 1)[::-1]
            if strain_name == "strain":
                continue
            if 'zoonotic' in lrline[6]:
                continue
            if strain_name not in dict_name_genomeID:
                continue
            genomeID = dict_name_genomeID[strain_name]
            if genomeID not in dict_genomeID_seg_seq:
                continue

            if seg in dict_genomeID_seg_seq[genomeID]:
                seq = dict_genomeID_seg_seq[genomeID][seg]
            else:
                # seq = dict_name_seg_seq_lack[name][seg]
                continue

            inputs = seqProcess(seq, args)
            prob_origin, logits_ori = predictor(inputs)
            # logits_ori = logits_ori.detach().numpy()[0]
            prob_float = float(prob_origin[0])
            print(seg, strain_name, prob_float)
            if genomeID not in dict_reassort_prob:
                dict_reassort_prob[genomeID] = {'PB2': '', 'PB1': '', 'PA': '', 'HA': '',
                                                'NP': '', 'NA': '', 'MP': '', 'NS': '', 'info': lrline[6]}
            dict_reassort_prob[genomeID][seg] = str(prob_float)
    dict_reassort_prob = pd.DataFrame.from_dict(dict_reassort_prob)
    dict_reassort_prob = dict_reassort_prob.T
    dict_reassort_prob.to_csv(f'results/zoonotic_strain/zoonotic_strain_3parts_smallsub2_not.csv')
    print('\n++++++++++++++++++')
