#%%
'''
Takes in alignment and alignment index.
Plots distribution of lengths for each alignment group.
'''
import csv, sys, re
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline


#%%
#Defining some helpful functions
def load_pid_data_from_txt(txt_location):
    '''Loads text data into dictionary. Returns it and a list of all unique sequence names.'''
    uniq_set = list()
    out_dict = dict()
    txt_file = open(txt_location)
    txt_list = txt_file.readlines()
    species_pairs = [re.sub(r'\/.*?\n', '', z) for z in
                    [re.sub(r'\/.*?,', ',', y) for y in
                    [re.sub(' vs ', ',', x) for x in
                    [re.sub('Results for ', '', w) for w in txt_list[1::7]]]]]
    perc_identity = [re.sub('  Percent identity: ','', y) for y in 
                    [re.sub('\n','', x) for x in txt_list[5::7]]]
    i = 0
    for pair in species_pairs:
        uniq_set.append(pair.split(',')[0])
        uniq_set.append(pair.split(',')[1])
        out_dict[(pair.split(',')[0],pair.split(',')[1])] = round(float(perc_identity[i]))
        i+=1
    return sorted(list(set(uniq_set)), key=lambda x: re.sub('_.*_','',x)), out_dict


def load_pid_data(csv_location):
    '''Loads csv data into dictionary. Returns it and a list of all unique sequence names.'''
    first_line = True
    uniq_set = list()
    out_dict = dict()
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            if first_line:
                first_line = False
                continue
            uniq_set.append(row[0])
            uniq_set.append(row[1])
            out_dict[(row[0],row[1])] = round(float(row[3]))
    return sorted(list(set(uniq_set)), key=lambda x: re.sub('_.*_','',x)), out_dict

def generate_matrix(uniq_set, pid_dict):
    pid_matrix = list()
    for seq1 in uniq_set:
        temp_list = list()
        for seq2 in uniq_set:
            if (seq1,seq2) not in pid_dict.keys() and (seq2,seq1) not in pid_dict.keys():
                temp_list.append(100)
                continue
            elif (seq1,seq2) not in pid_dict.keys():
                temp_list.append(pid_dict[(seq2,seq1)])
                continue
            elif (seq2,seq1) not in pid_dict.keys():
                temp_list.append(pid_dict[(seq1,seq2)])
                continue
        pid_matrix.append(temp_list)
    return pid_matrix

#%%
#Set to True if running for the eukaryotic-only alignment
eukarya = True

#%%

#uniq_set, pid_dict = load_pid_data('./pids/ES39_Euk.csv')
#uniq_set, pid_dict = load_pid_data('./pids/ES39_Asg.csv')
#uniq_set, pid_dict = load_pid_data('./pids/ES39_SAR+EXC_Asg_ClustalO.csv')
#uniq_set, pid_dict = load_pid_data('./pids/ES39_SAR+EXC_Asg_Probcons.csv')
#uniq_set, pid_dict = load_pid_data_from_txt('./txts/ES39_Euk_manual.txt')
#uniq_set, pid_dict = load_pid_data_from_txt('./txts/ES39_SAR+EXC_Asg_ClustalO_manual_v2.txt')
uniq_set, pid_dict = load_pid_data_from_txt('./txts/CntrProt_Euk.txt')

#%%
#Special ordering for the Eukarya-only alignment
if eukarya: 
    uniq_set = ['LSUe_9606_HUMAN','LSUe_9598_PANTR','LSUe_10090_MOUSE','LSUe_10116_RAT','LSUe_13616_MONDO',
                'LSUe_9031_CHICK','LSUe_28377_ANOCA','LSUe_8355_XENLA','LSUe_7955_DANRE','LSUe_7897_LATCH',
                'LSUe_6238_CAEBR','LSUe_6239_CAEEL','LSUe_7227_DROME','LSUe_7160_AEDAL','LSUe_559292_YEAST',
                'LSUe_33169_EREGO','LSUe_284812_SCHPO','LSUe_284591_YARLI','LSUe_45157_CYAME','LSUe_3702_ARATH',
                'LSUe_39947_ORYSJ','LSUe_55529_GUITH','LSUe_36329_PLAF7','LSUe_353151_CRYHT','LSUe_554453_PLAMI',
                'LSUe_753081_BIGNC','LSUe_1928728_PAMIC','LSUe_3055_CHLRE','LSUe_5911_TETTH','LSUe_857967_ICHMG',
                'LSUe_1333499_NAOCE','LSUe_2044563_ANTWI','LSUe_1315715_PIJUD','LSUe_238746_MICPO','LSUe_296543_THAPS',
                'LSUe_35144_PRYPA','LSUe_420245_LEIBR','LSUe_185431_TRYB2','LSUe_352472_DICDI','LSUe_104782_ADIVA']

#%%

mx = generate_matrix(uniq_set, pid_dict)

#%%

plt.figure(figsize=(25, 25))

g = sns.heatmap(mx, cmap="viridis", annot=True, annot_kws={"size": 14}, fmt='d', xticklabels=uniq_set, yticklabels=uniq_set, linewidths=.1)

plt.rcParams["axes.labelsize"] = 15
g.set_ylim(len(mx), 0)
g.set_xticklabels(g.get_xticklabels(), rotation=45, ha='right', fontsize=14)
g.set_yticklabels(g.get_yticklabels(), fontsize=14)

g.figure.axes[-1].set_ylabel('% identities', size=20)

# %%
