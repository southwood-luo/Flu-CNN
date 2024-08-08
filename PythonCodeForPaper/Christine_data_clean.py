import pandas as pd
import pickle


def load_dict(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


dict_name_genomeID = load_dict('../dict/dict_name_genomeID')
fn = f'../data/other_try/zoonotic_strain.csv'
# meta0 = pd.read_csv(fn)
fread = open(fn, 'r')
for rline in fread:
    lrline = rline.rstrip().split(',')
    name = lrline[1]
    strain_name = name.replace(' ', '_').replace('-', '_').lower()
    strain_name = strain_name[::-1].replace(' ', '', 1)[::-1]
    if strain_name == "strain":
        lrline.insert(1, 'genomeID')
        with open('../data/other_try/zoonotic_strain2.csv','a') as fw:
            fw.write(','.join(lrline)+'\n')
        continue
    if strain_name not in dict_name_genomeID:
        continue
    genomeID = dict_name_genomeID[strain_name]
    lrline.insert(1, genomeID)
    with open('../data/other_try/zoonotic_strain2.csv', 'a') as fw:
        fw.write(','.join(lrline)+'\n')
