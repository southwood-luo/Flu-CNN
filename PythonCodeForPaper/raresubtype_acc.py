import argparse
import time
import sys

import numpy as np
from model import CharCNN
import torch
from torch.autograd import Variable
import torch.nn.functional as F
import pandas as pd
import copy
from Bio import SeqIO
import os

os.environ["CUDA_VISIBLE_DEVICES"] = '1'

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


aaLigal = 'ACDEFGHIKLMNPQRSTVWY'
ntLigal = 'ATCG'

dictCodon = dict()
fread = open('data/codon.txt', 'r')
for rline in fread:
    lrline = rline.strip().split('\t')
    aa = lrline[1]
    codons = lrline[0].split(',')
    for codon in codons:
        if len(codon) != 3:
            print('error codon!')
        else:
            dictCodon[codon] = aa

fread.close()


def translateCDS(seq):
    length = len(seq)
    if length % 3 != 0:
        print('wrong codon number, cannot divided by 3!')
        sys.exit()
    pSeq = ''
    for i in range(length):
        start = i * 3
        codon = seq[start:start + 3]
        if codon in dictCodon:
            aa = dictCodon[codon]
            if aa not in aaLigal or aa == 'STOP':
                break
            else:
                pSeq = pSeq + aa
        else:
            break
    return pSeq

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


def predictor(inputs):
    inputs = Variable(inputs)
    logit = model(inputs)
    probs = F.softmax(logit, dim=1)[0]
    return probs


if __name__ == '__main__':

    segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    # segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    # segs = [, 'all8']
    subs = ['h5n1', 'h7n9', 'h9n2']
    args.data_type = 'protein'
    args.mafft = 'unmafft'
    for sub in subs:
        print(sub)
        dict_acc_prob = dict()
        for seg in segs:
            args.seg = seg
            print('\n节段：' + seg)
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
                args.l0 = seg_max_length_protein[seg]
                args.alphabet_path = 'alphabet_protein.json'
                args.num_features = 27
                # args.alphabet = 'merikldsqtvhaygnpwfc'
                args.alphabet = 'abcdefghijklmnopqrstuvwxyz'

            if args.class_type == 'host':
                args.num_class = 2
                args.val_interval = 100
                # args.model_path = f'models_CharCNN/111and114_unmafft/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll.pth.tar'
                # args.model_path = f'models_CharCNN/3parts_smallsub/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll_3parts_smallsub.pth.tar'
                # args.model_path = f'models_CharCNN/3parts/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll_3parts.pth.tar'
                args.model_path = f'models_CharCNN/0103_unmafft_model/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll.pth.tar'
                args.test_path = f'data/rare_subtype_fas/{sub}/{sub}_unmafft_dna_{seg}.fas'

            # load testing data
            print("\nLoading testing data...")
            print(args.test_path)
            model = CharCNN(args)
            print("=> loading weights from '{}'".format(args.model_path))
            assert os.path.isfile(args.model_path), "=> no checkpoint found at '{}'".format(args.model_path)
            checkpoint = torch.load(args.model_path)
            model.load_state_dict(checkpoint['state_dict'])
            model.eval()

            # fw = open(f'results/3parts/3parts_test_fa_smallsub/3parts_unmafft_protein_{seg}_test_result_smallsub.txt', 'w')

            for record in SeqIO.parse(args.test_path, 'fasta'):
                strain_info = str(record.description)
                seqfas = str(record.seq).upper()
                seq = translateCDS(seqfas).lower()
                genomeID = str(strain_info.split('|')[0])
                # print(genomeID)
                host_group = strain_info.split('|')[3]
                strain_name = strain_info.split('|')[1]
                subtype = strain_info.split('|')[2]
                inputs = seqProcess(seq, args)
                prob_origin = predictor(inputs)
                prob_float = float(prob_origin[0])
                if genomeID not in dict_acc_prob:
                    dict_acc_prob[genomeID] = {'PB2': '', 'PB1': '', 'PA': '', 'HA': '',
                                               'NP': '', 'NA': '', 'MP': '', 'NS': '',
                                               'host_group': host_group, 'strain_name': strain_name,
                                               'subtype': subtype}
                dict_acc_prob[genomeID][seg] = str(prob_float)
        dict_acc_prob = pd.DataFrame.from_dict(dict_acc_prob)
        dict_acc_prob = dict_acc_prob.T
        # dict_acc_prob.to_csv(f'results/3parts/3parts_raresub_smallsub/3parts_rare_result_smallsub_{sub}.csv')
        dict_acc_prob.to_csv(f'results/3parts/0103_rare_result_{sub}.csv')
        #     fw.write(genomeID + '\t' + strain_name + '\t' + host_group + '\t' + subtype + '\t' + str(prob_float) + '\n')
        # fw.close()

    print('\n++++++++++++++++++')
