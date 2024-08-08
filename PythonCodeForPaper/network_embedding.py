import argparse
import pickle
from model import CharCNN
import torch
from torch.autograd import Variable
import torch.nn.functional as F
from Bio import SeqIO
import os
import sys
import psutil
import gc
import objgraph

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


def load_dict(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def save_dict(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


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


# 做个钩子，抓取中间层

def hook(module, fea_in, fea_out):
    features_in_hook.append(fea_in)
    features_out_hook.append(fea_out)
    fea_in = None
    fea_out = None
    return None


def predictor(inputs):
    inputs = Variable(inputs)
    logit = model(inputs)
    # for (name, module) in model.named_modules():
    #     if name == layer_name:
    #         module.register_forward_hook(hook=hook)
    # model.remove()
    probs = F.softmax(logit, dim=1)[0]
    return probs


if __name__ == '__main__':
    # segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS','all8']
    # segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    segs = ['NA']
    args.data_type = 'protein'
    args.mafft = 'unmafft'
    for seg in segs:
        if os.path.exists(f'results/internal_layer/features_in_hook_fc1_{seg}.pkl'):
            features_in_hook = load_dict(f'results/internal_layer/features_in_hook_fc1_{seg}')
            features_out_hook = load_dict(f'results/internal_layer/features_out_hook_fc1_{seg}')
        else:
            features_in_hook = []
            features_out_hook = []
        args.seg = seg
        print('\n看节段的宿主类型：' + seg)

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
            args.model_path = f'models_CharCNN/3parts_smallsub/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll_3parts_smallsub.pth.tar'
            # args.model_path = f'models_CharCNN/3parts/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll_3parts.pth.tar'
            # args.model_path = f'models_CharCNN/0103_unmafft_model/CharCNN_host_{args.data_type}_{seg}_{args.mafft}_best_nll.pth.tar'
            args.test_path = f'data/3parts/3parts_test_fa/3parts_8seg_unmafft_protein_{seg}_test.fa'

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
        n = 0
        out_num = len(features_in_hook)
        conti_flag = False
        for record in SeqIO.parse(args.test_path, 'fasta'):
            n += 1
            if n <= out_num:
                continue
            if n % 1000 == 0 and conti_flag:
                # print(objgraph.show_growth())
                # objgraph.show_most_common_types(limit=10)
                # tensor_objects = [obj for obj in gc.get_objects() if isinstance(obj, torch.Tensor)]
                #
                # for i in range(len(tensor_objects)):
                #     del tensor_objects[0]
                # torch.cuda.empty_cache()
                # gc.collect()
                # tensor_objects = [obj for obj in gc.get_objects() if isinstance(obj, torch.Tensor)]

                # print(len(tensor_objects))
                mem = psutil.virtual_memory()
                available = mem.available
                free = mem.free
                total = mem.total
                print(f"可用内存：{available / 1024 / 1024:.2f} MB")
                print(f"空闲内存：{free / 1024 / 1024:.2f} MB")
                print(f"总内存：{total / 1024 / 1024:.2f} MB")
                # print(f"handle内存：{sys.getsizeof(handle) / 1024 / 1024} MB")
                print('##')
                save_dict(features_in_hook, f'results/internal_layer/features_in_hook_fc1_{seg}')
                save_dict(features_out_hook, f'results/internal_layer/features_out_hook_fc1_{seg}')
            conti_flag = True
            handle = model.fc1.register_forward_hook(hook)
            strain_info = str(record.description)
            genomeID = str(strain_info.split('|')[0])

            seq = str(record.seq).lower()
            host_group = strain_info.split('|')[2]
            strain_name = strain_info.split('|')[1]
            subtype = strain_info.split('|')[3]
            inputs = seqProcess(seq, args)
            prob_origin = predictor(inputs)
            # print(n)
            handle.remove()
            gc.collect()
            prob_float = float(prob_origin[0])
            # fw.write(genomeID + '\t' + strain_name + '\t' + host_group + '\t' + subtype + '\t' + str(prob_float) + '\n')
        # fw.close()
        save_dict(features_in_hook, f'results/internal_layer/features_in_hook_fc1_{seg}')
        save_dict(features_out_hook, f'results/internal_layer/features_out_hook_fc1_{seg}')



    print('\n++++++++++++++++++')

