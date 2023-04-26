import pandas as pd
from sklearn.metrics import classification_report

df_vidhop = pd.read_csv('data/vidhop.csv', index_col=1).assign(Method='VIDHOP')
df_vidhop = df_vidhop.drop(df_vidhop.columns[0], axis=1)
df_vidhop = df_vidhop.rename(columns={'host_pred': 'Mix'})
df_dnts = pd.read_csv('data/dnts.csv', index_col=0).assign(Method='ML-DNTs')
df_dnts['new_column_name'] = df_dnts.index.str.split('|').str.get(0)
df_dnts.index = df_dnts.new_column_name
df_dnts = df_dnts.rename(columns={'host_pred': 'Mix'})
df_cnn = pd.read_csv('data/flucnn_results.csv', index_col=1).assign(Method='Flu-CNN')
df_cnn = df_cnn.drop(df_cnn.columns[0], axis=1)
df_cnn = df_cnn.rename(columns={'host_pred': 'Mix'})
df_result = pd.DataFrame(columns=['Method', 'Segment', 'Precision', 'Recall', 'F1-score', 'Accuracy'])

df = {'VIDHOP':df_vidhop, 'ML-DNTs':df_dnts, 'Flu-CNN': df_cnn}
segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS', 'Mix']

for method_name in df:
    for seg in segs:
        method_df = df[method_name]
        y_true = method_df['host_true']
        if seg not in method_df:
            continue
        y_pred = method_df[seg]
        report = classification_report(y_true, y_pred, output_dict=True)
        row = {'Segment': seg,
               'Method': method_name,
               'Precision': report['weighted avg']['precision'],
               'Recall': report['weighted avg']['recall'],
               'F1-score': report['weighted avg']['f1-score'],
               'Accuracy': report['accuracy']}
        df_result = df_result.append(row, ignore_index=True)  # 将行添加到数据框中

df_result = df_result.iloc[:, [1, 0] + list(range(2, len(df_result.columns)))]
df_result['Segment'] = pd.Categorical(df_result['Segment'], categories=segs)
df_result = df_result.sort_values(by=['Segment', 'Method'])
df_result = df_result.T
df_result.to_csv('test_metrics.csv', header=False)
