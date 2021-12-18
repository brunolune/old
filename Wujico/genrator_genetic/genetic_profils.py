#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:57:09 2019

@author: bruno
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import euclidean_distances

from genetic_data import profil_min,profil_max,columns_nutri

profil_min_max = (profil_min + profil_max)
profil_min_max /= 2

#%%
def plot_dendrogram(model, index, **kwargs):
     # Children of hierarchical clustering
     children = model.children_
     # Distances between each pair of children (Since we don't have this information, we can use a uniform one for plotting)
     distance = np.arange(children.shape[0])
     # The number of observations contained in each cluster level
     no_of_observations = np.arange(2, children.shape[0]+2)
     # Create linkage matrix and then plot the dendrogram
     linkage_matrix = np.column_stack([children, distance,no_of_observations]).astype(float)
     # Plot the corresponding dendrogram
     fig, axes = plt.subplots(figsize=(20, 15))
     dendrogram(linkage_matrix, **kwargs, labels=index )
     #print(model.labels_)
#%%
# Calculate the mean between min and max
def det_grouped_profils(profils_target,dist_threshold):

    profil_selected,poids_demo = zip(*profils_target)
    #load selected_profiles to resume from a previous optimization
    #with open('log_GA/2019_10_09/Poitiers/saved_populations/profil_selected_2019_10_09_12_11_08_pop_300_pc_0.95_pm_0.4_cn_0.6_ni_251_iter_250_pn_0_testplotconv.pkl','rb') as f:
    #        profil_selected= pickle.load(f)
    profil_min_max_selected = profil_min_max.loc[ profil_selected , : ]
    profil_min_selected = profil_min.loc[ profil_selected , : ]
    profil_max_selected = profil_max.loc[ profil_selected , : ]
    
    #distance matrix
    distance_matrix=euclidean_distances(profil_min_max_selected)
    print(distance_matrix)
    x,y=np.meshgrid(poids_demo,poids_demo)
    mat_poids=x*y
    print(mat_poids)
    # Compute the agglomerative clustering on the mean between min and max
    clt = AgglomerativeClustering(linkage='complete',affinity='precomputed',
                                   distance_threshold = dist_threshold*np.mean(mat_poids),
                                   n_clusters = None,
                                   compute_full_tree = True)
    
    model = clt.fit(distance_matrix*mat_poids)
    #plot_dendrogram(model, profil_min_max_selected.index)
    
    # Assign the newly determined clusters
    profil_min_selected["cluster"] = model.labels_
    profil_max_selected["cluster"] = model.labels_
    #profil_min_selected=profil_min_selected.iloc[0:2,:]
    #profil_max_selected=profil_min_selected.iloc[0:2,:]
    # Cree des profils moyens a partir des clusters
    profil_min_grouped = profil_min_selected.groupby(['cluster']).mean()
    profil_max_grouped = profil_max_selected.groupby(['cluster']).mean()
    #list les noms de profils regroupés
    gp=profil_min_selected.reset_index().groupby(['cluster'])
    profil_grouped_profil_names=gp['nom'].apply(list)
    profil_grouped_profil_nbs=[gp.groups[i].to_list() for i in range(len(gp))]
    profil_grouped_names=pd.DataFrame({'names':profil_grouped_profil_names,'profils':pd.Series(profil_grouped_profil_nbs)})
    return profil_min_grouped,profil_max_grouped,profil_grouped_names
#%%
#   fix seed for reproducibility
#   np.random.seed(123)
    
#   nb_profil=10
profil_selected =list(profil_min.index[[0,36,20,46,9,28,51,42,7]]) #random.sample(list(profil_min.index), nb_profil)
poids_demo=[1,1,1,1,1,1,1,1,1]#[30,30,5,5,2,3,3,5,5] #np.random.randint(1,30,nb_profil)
profils_target=zip(profil_selected,poids_demo)
dist_threshold=650
profil_min_grouped,profil_max_grouped,profil_grouped_names=det_grouped_profils(profils_target,dist_threshold)
ratio_min_max= profil_max_grouped[columns_nutri].div(profil_min_grouped[columns_nutri])