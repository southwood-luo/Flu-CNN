import pandas as pd
import pickle
import numpy as np


def load_dict(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


dict_name_hostgroup = load_dict('../dict/dict_name_hostgroup')
dict_genomeID_name = load_dict('../dict/dict_genomeID_name')
dict_seg_subtype_time_genomeID_sample = load_dict('../dict/rare_subtype_sample/dict_seg_subtype_time_genomeID_sample')
segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
segs2 = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA']
subtypes = ['h9n2', 'h5n1', 'h7n9']

for subtype in subtypes:
    with open(f'hostdetection_res/hostdetection_mothods_re_{subtype}.csv', 'a') as fw:
        fw.write('Sample Name' + ',' + 'Total Sequences' + ',' + 'Gene' + ',' + 'Detection Model' + ',' + 'Accuracy Rate' + '\n')

# vidhop
for subtype in subtypes:
    for seg in segs:
        fn = f'vidhop_raresub/vidhop_pred_{subtype}_{seg}.csv'
        data = pd.read_csv(fn, index_col=0)
        for i in range(1, 21, 1):
            sample_genomeID = dict_seg_subtype_time_genomeID_sample[seg][subtype][str(i)]
            sample_genomeID = [int(i) for i in sample_genomeID]
            sample_data = data.loc[sample_genomeID]

            equal_count = (sample_data['host_true'] == sample_data['host_pred']).sum()
            percentage = equal_count / len(sample_data)
            with open(f'hostdetection_res/hostdetection_mothods_re_{subtype}.csv', 'a') as fw:
                fw.write(subtype.upper()+'_'+seg+'_resample_'+str(i) + ',' + str(len(sample_data)) + ',' + seg + ',' +
                         'VIDHOP' + ',' + str(percentage) + '\n')


# dnts
for subtype in subtypes:
    for seg in segs2:
        fn = f'dnts_raresub/results_{subtype}_{seg}.csv'
        data = pd.read_csv(fn, index_col=0)
        data['new_column_name'] = data.index.str.split('|').str.get(0)
        data['host_true'] = data.index.str.split('|').str.get(2)
        data.index = data.new_column_name

        for i in range(1, 21, 1):
            sample_genomeID = dict_seg_subtype_time_genomeID_sample[seg][subtype][str(i)]
            # sample_genomeID = [int(i) for i in sample_genomeID]
            sample_data = data.loc[sample_genomeID]

            equal_count = (sample_data['host_true'] == sample_data['host_pred']).sum()
            percentage = equal_count / len(sample_data)
            with open(f'hostdetection_res/hostdetection_mothods_re_{subtype}.csv', 'a') as fw:
                fw.write(subtype.upper()+'_'+seg+'_resample_'+str(i) + ',' + str(len(sample_data)) + ',' + seg + ',' +
                         'ML-DNTs' + ',' + str(percentage) + '\n')

# earlybird
for subtype in subtypes:
    for seg in segs:
        fn = f'earlybird_web/{subtype}/earlybird_{subtype}_{seg}.csv'
        data = pd.read_csv(fn, index_col=0)
        data['genomeID'] = data.index
        data['host_true'] = data['genomeID'].apply(lambda x: dict_name_hostgroup[dict_genomeID_name[str(x)]])
        data['host_pred'] = data['host_pred'].apply(lambda x: x.lower())

        for i in range(1, 21, 1):
            sample_genomeID = dict_seg_subtype_time_genomeID_sample[seg][subtype][str(i)]
            sample_genomeID = [int(i) for i in sample_genomeID]
            sample_data = data.loc[sample_genomeID]

            equal_count = (sample_data['host_true'] == sample_data['host_pred']).sum()
            percentage = equal_count / len(sample_data)
            with open(f'hostdetection_res/hostdetection_mothods_re_{subtype}.csv', 'a') as fw:
                fw.write(subtype.upper()+'_'+seg+'_resample_'+str(i) + ',' + str(len(sample_data)) + ',' + seg + ',' +
                         'FluPhenotype' + ',' + str(percentage) + '\n')

# phynogenetic
for subtype in subtypes:
    fn = f'Sampletest_all/{subtype}_results.csv'
    data = pd.read_csv(fn)
    data.replace(np.nan, 'NA', inplace=True)

    for data_row in data.iterrows():
        with open(f'hostdetection_res/hostdetection_mothods_re_{subtype}.csv', 'a') as fw:
            fw.write(data_row[1][0] + ',' + str(data_row[1][1]) + ',' + data_row[1][2] + ',' +
                     'Phylogenic' + ',' + str(data_row[1][4]) + '\n')

# char-cnn
for subtype in subtypes:
    for seg in segs:
        # fn = f'char-cnn_raresub/0103_rare_result_{subtype}.csv'
        fn = f'char-cnn_raresub/3parts_rare_result_smallsub_{subtype}.csv'
        data = pd.read_csv(fn, index_col=0)
        data = data[[seg, 'host_group']]
        data = data.rename(columns={'host_group': 'host_true', seg: 'prob'})
        data.dropna()
        data['host_pred'] = data['prob'].apply(lambda x: 'human' if float(x) >= 0.5 else 'avian')

        for i in range(1, 21, 1):
            sample_genomeID = dict_seg_subtype_time_genomeID_sample[seg][subtype][str(i)]
            sample_genomeID = [int(i) for i in sample_genomeID]
            sample_data = data.loc[sample_genomeID]

            equal_count = (sample_data['host_true'] == sample_data['host_pred']).sum()
            percentage = equal_count / len(sample_data)
            with open(f'hostdetection_res/hostdetection_mothods_re_{subtype}.csv', 'a') as fw:
                fw.write(subtype + '_' + seg + 'resample_' + str(i) + ',' + str(len(sample_data)) + ',' + seg + ',' +
                         'Char-CNN' + ',' + str(percentage) + '\n')

# for subtype in subtypes:
#     for seg in segs:
#         fn = f'char-cnn_raresub/0103_rare_result_{subtype}.csv'
#         data = pd.read_csv(fn, index_col=0)
#         data = data[[seg, 'host_group']]
#         data = data.rename(columns={'host_group': 'host_true', seg: 'prob'})
#         data.dropna()
#         data['host_pred'] = data['prob'].apply(lambda x: 'human' if float(x) >= 0.5 else 'avian')
#
#         for i in range(1, 21, 1):
#             sample_genomeID = dict_seg_subtype_time_genomeID_sample[seg][subtype][str(i)]
#             sample_genomeID = [int(i) for i in sample_genomeID]
#             sample_data = data.loc[sample_genomeID]
#
#             equal_count = (sample_data['host_true'] == sample_data['host_pred']).sum()
#             percentage = equal_count / len(sample_data)
#             with open(f'hostdetection_res/hostdetection_mothods_re_{subtype}.csv', 'a') as fw:
#                 fw.write(subtype + '_' + seg + 'resample_' + str(i) + ',' + str(len(sample_data)) + ',' + seg + ',' +
#                          'Char-CNN' + ',' + str(percentage) + '\n')
#
print('**'*20)
