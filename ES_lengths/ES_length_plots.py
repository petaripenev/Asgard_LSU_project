#%%
'''
Takes in alignment and alignment index.
Plots distribution of lengths for each alignment group.
'''
import csv, sys, re
import numpy as np
import matplotlib
import pandas as pd
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqUtils import GC
import seaborn as sns
sns.set_style('whitegrid')
%matplotlib inline
#%%
alignment = AlignIO.read(open('../data/Supp_Data_S1_global_LSU_MSA_v5.fa'), "fasta")

#%%

def hs_index_to_aln_index(alignment_obj, index_range):
    '''Given a sequence range of HS and an alignment object,
    returns a truncated alignment by the HS index.
    '''
    index_start=0
    index_end=0
    for aln in alignment_obj:
        if aln.id == 'LSUe_9606_HUMAN/1-5227':
            aln_anchor_map=dict()
            aln_i = 0
            seq_i = 0
            for letter in aln:
                if letter != '-':
                    seq_i+=1
                    aln_anchor_map[seq_i] = aln_i
                aln_i +=1
            index_start=int(str(aln.seq[:aln_anchor_map[index_range[0]]]).count('-')) + index_range[0] - 1
            index_end=int(str(aln.seq[:aln_anchor_map[index_range[1]]]).count('-')) + index_range[1]
    return alignment[:, index_start:index_end]

#Construct truncated alignment objects covering the ES regions
es39_alignment = hs_index_to_aln_index(alignment, (4891,5123))
es9_alignment = hs_index_to_aln_index(alignment, (1578,1621))
es27_alignment = hs_index_to_aln_index(alignment, (3052,3763))
es7_alignment = hs_index_to_aln_index(alignment, (597,1464))

#%%
#Check that ES39 alignment is OK
print(es39_alignment)

#%%
#Quick GC content
for seq in alignment:
    if re.search(r'LSUASG', seq.id):
        print(GC(str(seq.seq).replace('-','')), seq.id)

#%%
def count_length(seq):
    allowed_bases = ['A','T','G','C','U']
    seq = seq.upper()
    total_dna_bases = 0
    for base in allowed_bases:
        total_dna_bases = total_dna_bases + seq.count(base.upper())
    return total_dna_bases

#%%
#Loki names
lokiarchaeota_names = [
    'LSUASG_SHMU01000207.1:1-2097/1-1912',
    'LSUASG_1538547_LOKI/1-3298',
    'LSUASG_F3H4_G_18380/2268-5563',
    'LSUASG_NJBH01000007.1:89956-93415/1-3275',
    'LSUASG_SDNK01000007.1:46432-48523_SDNK01000056.1:1-1270/1-3259',
    'LSUASG_SDNY01000005.1:53412-56906/1-3294',
    'LSUASG_QMYW01000148.1:199-3284/1-3085',
    'LSUASG_CP042905.1:656088-659193/1-3106']

#%%
#For renaming purposes (to look nicer on plot)
nice_name_dict = {
    'LSUb' : 'Bacteria',
    'LSUa' : 'Archaea',
    'LSUASG' : 'Asgard',
    'LSUe' : 'Eukarya',
}

#Initiate lists for populating pandas
group_order = []
es7_length_order = []
es9_length_order = []
es27_length_order = []
es39_length_order = []
rRNA_23S = []
name_order = []
lokiarchaeota = []

#Populate the lists with filtered data
i = 0
for aln in alignment:
    name_order.append(aln.name)
    group_order.append(nice_name_dict[aln.name.split('_')[0]])
    if aln.name in lokiarchaeota_names:
        lokiarchaeota.append(True)
    else:
        lokiarchaeota.append(False)
    if ((count_length(alignment[i].seq) == 1073) or     #Filter for sequences covering only ES9
        (count_length(alignment[i].seq) == 981) or
        (count_length(alignment[i].seq) == 959) or
        (count_length(alignment[i].seq) == 1227)):
        es7_length_order.append(None)
        es9_length_order.append(count_length(es9_alignment[i].seq))
        es27_length_order.append(None)
        rRNA_23S.append(None)
        es39_length_order.append(None)
        i+=1
        continue
    if count_length(alignment[i].seq) <= 2100:          #Filter for sequences covering only ES39
        print(alignment[i].id)
        es7_length_order.append(None)
        es9_length_order.append(None)
        es27_length_order.append(None)
        rRNA_23S.append(None)
        es39_length_order.append(count_length(es39_alignment[i].seq))
    if count_length(alignment[i].seq) > 2100:           #Full length sequences filtered here
        es7_length_order.append(count_length(es7_alignment[i].seq))
        es9_length_order.append(count_length(es9_alignment[i].seq))
        es27_length_order.append(count_length(es27_alignment[i].seq))
        rRNA_23S.append(count_length(alignment[i].seq))
        es39_length_order.append(count_length(es39_alignment[i].seq))
    i+=1

#Organize lists in dictionary of lists, pass them to a pandas dataframe
df = pd.DataFrame({
    'Phylogenetic Group': group_order,
    'Name': name_order,
    'ES7 Nucleotide Length': es7_length_order,
    'ES9 Nucleotide Length': es9_length_order,
    'ES27 Nucleotide Length': es27_length_order,
    'ES39 Nucleotide Length': es39_length_order,
    'Full LSU rRNA length': rRNA_23S,
    'Lokiarcheota' : lokiarchaeota,
})

#Check if everything looks in order
df

#%%
#Plot figures
plt.figure(figsize=(16,12))
g = sns.boxplot(data=df, x='ES39 Nucleotide Length', y='Phylogenetic Group')
g = sns.swarmplot(data=df, x='ES39 Nucleotide Length', y='Phylogenetic Group', color=".25")

#%%
#Output tables with organized data
pd.set_option('display.max_colwidth', -1)
pd.set_option('display.max_rows', 500)
df[df['Name'].str.match('LSUa')].sort_values(by=['ES39 Nucleotide Length'], ascending=False)
#print(df['Name'])

#%%
#Ugly barplot
#g = sns.barplot(data=df[df.Name.str.contains('LSUe|LSUASG')], x='ES9 Nucleotide Length', y='Name',hue='Phylogenetic Group')

#%%
#Plot ES39
#palette = ["#808080", "#499331", "#7f0089", '#000000']
palette = sns.color_palette("Set2", n_colors=3)
plt.figure(figsize=(16,12))
g = sns.boxplot(data=df, x='ES39 Nucleotide Length', y='Phylogenetic Group')
g = sns.swarmplot(data=df, x='ES39 Nucleotide Length', y='Phylogenetic Group', hue='Lokiarcheota', palette=sns.color_palette(palette), size=6)

#%%
palette = ["#808080", "#7f0089"]
plt.figure(figsize=(16,12))
#g = sns.boxplot(data=df, x='ES9 Nucleotide Length', y='Phylogenetic Group')
g = sns.swarmplot(data=df, x='ES9 Nucleotide Length', y='Phylogenetic Group', hue='Lokiarcheota', palette=sns.color_palette(palette), size=10)

#%%
plt.figure(figsize=(16,12))
#g = sns.boxplot(data=df, x='Full LSU rRNA length', y='Phylogenetic Group')
g = sns.boxplot(data=df[df.Name.str.contains('LSUa|LSUb')], x='Full LSU rRNA length', y='Phylogenetic Group')


#%%
#LSUe data
#df[df['Name'].str.match('LSUe')][['ES7 Nucleotide Length','ES9 Nucleotide Length','ES27 Nucleotide Length', 'ES39 Nucleotide Length','Name', 'Full LSU rRNA length']]
#df[df['Name'].str.match('LSUe')][['ES7 Nucleotide Length','ES9 Nucleotide Length','ES27 Nucleotide Length', 'ES39 Nucleotide Length','Name', 'Full LSU rRNA length']].sort_values(by=['Full LSU rRNA length'], ascending=False).to_csv(r'./ES_lengths.csv')

#%%
#All data

df[['ES7 Nucleotide Length','ES9 Nucleotide Length','ES27 Nucleotide Length', 'ES39 Nucleotide Length','Name', 'Full LSU rRNA length']].sort_values(by=['Full LSU rRNA length'], ascending=False)
#df[['Name', 'Full LSU rRNA length','ES7 Nucleotide Length','ES9 Nucleotide Length','ES27 Nucleotide Length', 'ES39 Nucleotide Length']].to_csv(r'./ES_lengths.csv')

#%%
palette = ["#000000", "#7f0089"]
plt.figure(figsize=(16,12))
#g = sns.boxplot(data=df, x='Full LSU rRNA length', y='Phylogenetic Group')
#g = sns.swarmplot(data=df, x='Full LSU rRNA length', y='Phylogenetic Group', hue='Lokiarcheota', palette=sns.color_palette(palette), size=6)
g = sns.stripplot(data=df, x='ES9 Nucleotide Length', y='Phylogenetic Group', hue='Lokiarcheota', size=20, jitter=0)

# %%
