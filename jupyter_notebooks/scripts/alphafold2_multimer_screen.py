#!/usr/bin/env python


import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import seaborn as sns
sns.set_style("ticks")

from matplotlib.pyplot import *
 
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')

from datetime import datetime, timedelta
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as multi



from statannotations.Annotator import Annotator


def manhattan_plot(df_all, gene):
    
    data_all = df_all.copy()    
    data_all = data_all[data_all.A == gene]
    data_all.index = range(len(data_all))        
        
    data_all, data_best = max_iptm(data_all, 0.75)
        
    data_grouped_all = data_best.groupby('species')['identifier'].median()
        
    plt.figure(figsize=(30, 25))

    fig = sns.relplot(data=data_all, x='identifier', y='iptm_ptm', aspect=3.7, hue='species', s=20, color='grey', alpha=0.3, legend=False).set(ylim=(0, 1))
    sns.scatterplot(data=data_best.rename(columns={'A':'gene_id'}), x='identifier', y='iptm_ptm', hue='species', style='gene_id', s=40, alpha=0.55, edgecolor='black', linewidth=1, legend='brief')

    plt.title('{}'.format(gene), fontsize=15)

    fig.ax.set_xlabel('Pathogen_secretome')
    fig.ax.set_xticks(data_grouped_all)
    fig.ax.set_xticklabels(data_grouped_all.index)
    plt.axhline(0.75, linestyle='--')
    plt.xticks(rotation=45)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    #plt.tight_layout()
    plt.show()
    
    
def rename_df(df, column, old, new):
    data = df.copy()
    data[column] = [x.replace(old, new) for x in data[column]]
    return(data)


def small_secreted_candidates(data):
    
    df = data.copy()

    # molecular weight filter at 35 kDa

    df = df.rename(columns={'mol_weight[kDa]':'mol_weight_kDa'})
    df = df[df.mol_weight_kDa < 35]

    
    # signalP and targetP selection
    
    df = df[(df.signalp_pred != "OTHER") | (df.targetp_pred == "SP")]
    df_apop = df[df.apop_pred == "Apoplastic"]
    
    return(df, df_apop)

def spec_df(df, name):
    
    new = []
    
    tmp = [x.replace("_", "") for x in df.index.values]
    new = pd.DataFrame(tmp, columns=['prot_id'])
    new['species'] = name
    
    return(new)

def max_iptm(df, max_value):
    filtered = []
    
    complexes = np.unique(df.prot_id.values)
    
    for comp in complexes:
        tmp = df[df.prot_id == comp]
        best = np.max(tmp.iptm_ptm.values)
        
        if best > 0:
            A = np.unique(tmp.A.values)[0]
            B = np.unique(tmp.B.values)[0]
            model = tmp[tmp.iptm_ptm == best]['model'].values[0]
            species = np.unique(tmp.species.values)[0]
            
            filtered.append([comp, A, B, model, best, species])
            
    filtered = pd.DataFrame(filtered, columns=['prot_id', 'A', 'B', 'model', 'iptm_ptm', 'species']).sort_values(by='species')   
    
    filtered.index = range(len(filtered))
    filtered['identifier'] = filtered.index
    
    to_screen = filtered.copy()
    
    for i in to_screen.index:
        
        if to_screen.loc[i, 'iptm_ptm'] < max_value:
            to_screen.drop([i], inplace=True)
            
    
            
    return(filtered, to_screen)



def avg_score_thresh(df, column, thresh):
    
    high = df[df[column] > thresh][column]
    return(pd.DataFrame(high, columns=[column]))



def screen_complexes(df, thresh):
    
    filtered = []
    
    complexes = np.unique(df.prot_id.values)
    for counter, comp in enumerate(complexes):
        tmp = df[df.prot_id == comp]
        
        max_score = np.max(tmp.iptm_ptm.values)
        
        if max_score > thresh:
            avg_score = np.mean(tmp.iptm_ptm.values)
            std_score = np.std(tmp.iptm_ptm.values)
            A = np.unique(tmp.A.values)[0]
            B = np.unique(tmp.B.values)[0]
            species = np.unique(tmp.species.values)[0]
            
            filtered.append([comp, A, B, avg_score, std_score, max_score, species])   
            
    filtered = pd.DataFrame(filtered, columns=['prot_id', 'A', 'B', 'avg_score', 'std_score', 'best_model_score', 'species']).sort_values(by='species')
    

        
    return(filtered)



def remove_pkls(df, thresh):
    
    filtered = []
    
    complexes = np.unique(df.prot_id.values)
    for counter, comp in enumerate(complexes):
        tmp = df[df.prot_id == comp]
        
        max_score = np.max(tmp.iptm_ptm.values)
        
        if max_score < thresh:
            avg_score = np.mean(tmp.iptm_ptm.values)
            std_score = np.std(tmp.iptm_ptm.values)
            A = np.unique(tmp.A.values)[0]
            B = np.unique(tmp.B.values)[0]
            species = np.unique(tmp.species.values)[0]
            
            filtered.append([comp, A, B, avg_score, std_score, max_score, species])   
            
    filtered = pd.DataFrame(filtered, columns=['prot_id', 'A', 'B', 'avg_score', 'std_score', 'best_model_score', 'species']).sort_values(by='species')
    
    return(filtered)


def manhattan_plot(df_all, gene):
    
    data_all = df_all.copy()    
    data_all = data_all[data_all.A == gene]
    data_all.index = range(len(data_all))        
        
    data_all, data_best = max_iptm(data_all, 0.75)
        
    data_grouped_all = data_best.groupby('species')['identifier'].median()
        
    plt.figure(figsize=(30, 25))

    fig = sns.relplot(data=data_all, x='identifier', y='iptm_ptm', aspect=3.7, hue='species', s=20, color='grey', alpha=0.3, legend=False).set(ylim=(0, 1))
    sns.scatterplot(data=data_best.rename(columns={'A':'gene_id'}), x='identifier', y='iptm_ptm', hue='species', style='gene_id', s=40, alpha=0.55, edgecolor='black', linewidth=1, legend='brief')

    plt.title('{}'.format(gene), fontsize=15)

    fig.ax.set_xlabel('Pathogen_secretome')
    fig.ax.set_xticks(data_grouped_all)
    fig.ax.set_xticklabels(data_grouped_all.index)
    plt.axhline(0.75, linestyle='--')
    plt.xticks(rotation=45)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    #plt.tight_layout()
    plt.show()
    
    
def subset_candidates(cand, pathogen, df):
    new = pd.DataFrame()
    for i in cand[cand.species == pathogen]['prot_ids'].to_numpy():
        new = new.append(df[df.index.str.contains(i)])
        
        
    return(new)
        