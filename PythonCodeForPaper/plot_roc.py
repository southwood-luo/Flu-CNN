import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def plot_roc(y_test, y_score_all, name):
    plt.figure()
    for i in range(len(y_score_all)):
        fpr, tpr, thread = roc_curve(y_test, y_score_all[i], pos_label=1)
        roc_auc = auc(fpr, tpr)
        print('roc_auc:', roc_auc)
        if segs[i] == 'all8':
            plt.plot(fpr, tpr, color='darkorange',
                     lw=2, label=f'Mix ROC curve (area = %0.5f)' % roc_auc)
        else:
            plt.plot(fpr, tpr, color='darkorange',
                     lw=2, label=f'{segs[i]} ROC curve (area = %0.5f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic for test dataset')
    plt.legend(loc="lower right")
    plt.savefig(f'{name}_roc.png', )
    plt.show()


segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS', 'all8']
data = pd.DataFrame()
host_true = []

prob_list_all = []
for seg in segs:
    print(seg)
    file = f'../results/unmafft_protein_{seg}_test_result.txt'
    prob_list = []
    with open(file, 'r') as fr:
        lines = fr.readlines()
        for line in lines:
            info = line.rstrip().split('\t')
            prob = info[4]
            prob_list.append(float(prob))
            if seg == 'PB2':
                host_true.append(1 if info[2] == 'human' else 2)
    print('##')
    prob_list_all.append(prob_list)
plot_roc(host_true, prob_list_all, f'figures/roc/')

print('######')
