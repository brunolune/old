#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:22:16 2019

@author: bruno
"""
import gym
from gym import spaces
import numpy as np
from gym.utils import seeding
import pandas as pd
from collections import defaultdict
import itertools
import os

class RLlibAIRecipesEnv(gym.Env):
    """ A recipes generator with Gym and RLlib"""
    metadata = {'render.modes': ['human']}
    
    def __init__(self, df_ingredients_used,df_profilmin,df_profilmax,pairs_values):
        
        self.df_ingredients_used = df_ingredients_used
        self.df_profilmin=df_profilmin
        self.df_profilmax=df_profilmax
        self.pairs_values=pairs_values
        
        self.nb_profils=1
        self.nb_indic_holi_jour=6 #SAIN,LIM,FF,PRAL,ORAC,IG,nutriscore
        
        #choix du repas
        self.repas='dej'
        #limitation apport nutritionnel au repas
        self.percent_nutrition={'pdej':0.25,'dej':0.3,'col':0.15,'din':0.3}
        self.df_profilmin=self.df_profilmin*self.percent_nutrition[self.repas] #30% pour dej
        self.df_profilmax=self.df_profilmax*self.percent_nutrition[self.repas] #30% pour dej
        #ratio
        self.df_profil_ratio=self.df_profilmax.div(self.df_profilmin)
        
        self.en_col_ind=self.df_ingredients_used.columns.get_loc("Energie, Règlement UE N° 1169/2011 (kcal)")
        
        
        
        
        
        #nb_nutriments si cols nutriments en fin de df_ingredients_used
        self.nb_nutriments=len(self.df_ingredients_used.columns[self.en_col_ind:])
        
        self.nb_cat_repas={'pdej':{'feculents':0,'legumineuses':1,'cereales':1,'legumes':0,
                        'vop':0,'fruits':2,'fruits_coque':1,'prod_laitiers':1,
                        'mat_grasses':1,'supplements':2},
                            'dej':{'feculents':1,'legumineuses':1,'cereales':1,'legumes':3,
                        'vop':1,'fruits':0,'fruits_coque':1,'prod_laitiers':1,
                        'mat_grasses':1,'supplements':0},
                            'col':{'feculents':1,'legumineuses':0,'cereales':0,'legumes':1,
                        'vop':0,'fruits':1,'fruits_coque':0,'prod_laitiers':1,
                        'mat_grasses':1,'supplements':2},
                            'din':{'feculents':0,'legumineuses':1,'cereales':1,'legumes':3,
                        'vop':1,'fruits':0,'fruits_coque':1,'prod_laitiers':0,
                        'mat_grasses':1,'supplements':2}}
                            
        #list of categories in dishes                    
        self.cat_dej=[[x]*self.nb_cat_repas['dej'][x] for x in self.nb_cat_repas['dej'].keys() if self.nb_cat_repas['dej'][x]!=0]
        self.cat_dej = [y for x in self.cat_dej for y in x] #to flatten list
        #reverse dictionnary of nb_cat_repas
        self.nb_repas_cat = defaultdict(dict)
        for key, val in self.nb_cat_repas.items():
            for subkey, subval in val.items():
                self.nb_repas_cat[subkey][key] = subval
        #nb of ingredients composing dishes by categories
        self.nb_feculents_jour=sum(self.nb_repas_cat['feculents'].values()) #self.nb_ingr_repas['pdej']['feculents']
        self.nb_legumineuses_jour=sum(self.nb_repas_cat['legumineuses'].values())
        self.nb_cereales_jour=sum(self.nb_repas_cat['cereales'].values())
        self.nb_legumes_jour=sum(self.nb_repas_cat['legumes'].values())
        self.nb_vop_jour=sum(self.nb_repas_cat['vop'].values())
        self.nb_fruits_jour=sum(self.nb_repas_cat['fruits'].values())
        self.nb_fruits_coque_jour=sum(self.nb_repas_cat['fruits_coque'].values())
        self.nb_prod_laitiers_jour=sum(self.nb_repas_cat['prod_laitiers'].values())
        self.nb_mat_grasses_jour=sum(self.nb_repas_cat['mat_grasses'].values())
        self.nb_supplements_jour=sum(self.nb_repas_cat['supplements'].values())
#        self.nb_tot_ingr_jour=self.nb_feculents_jour+self.nb_legumineuses_jour+self.nb_cereales_jour+self.nb_legumes_jour \
#                        +self.nb_vop_jour+self.nb_fruits_jour+self.nb_fruits_coque_jour+self.nb_prod_laitiers_jour \
#                        +self.nb_mat_grasses_jour+self.nb_supplements_jour
        self.nb_tot_ingr_jour=sum(self.nb_cat_repas['pdej'].values())+sum(self.nb_cat_repas['dej'].values())+ \
        sum(self.nb_cat_repas['col'].values())+sum(self.nb_cat_repas['din'].values())
                
        #listes catégories d'ingrédients
        self.list_categories_alim_code={'feculents':list(set(df_ingredients_used.groupby('alim_grp_nom_fr')['alim_code'].apply(list)['Féculents'])),
                              'legumineuses':list(set(df_ingredients_used.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['légumineuses'])),
                              'cereales':list(set(df_ingredients_used.groupby('alim_ssssgrp_nom_fr')['alim_code'].apply(list)['céréales'])) \
                                          +list(set(df_ingredients_used.groupby('alim_ssssgrp_nom_fr')['alim_code'].apply(list)['autres céréales'])),
                              'legumes':list(set(df_ingredients_used.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['légumes'])),
                              'vop':list(set(df_ingredients_used.groupby('alim_grp_nom_fr')['alim_code'].apply(list)['Viandes, œufs, poissons'])),
                              'fruits':list(set(df_ingredients_used.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['fruits'])),
                              'fruits_coque':list(set(df_ingredients_used.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['fruits à coque '])),
                              'prod_laitiers':list(set(df_ingredients_used.groupby('alim_grp_nom_fr')['alim_code'].apply(list)['Laits et produits laitiers'])),
                              'mat_grasses':list(set(df_ingredients_used.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['végétale'])),
                              'ingredients':list(set(df_ingredients_used['alim_code'].to_list()))}
        #liste des index (nouveaux index, pas ceux de Ciqual) 
        self.list_categories_index={'feculents':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['feculents'])].index.to_list(),
                              'legumineuses':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['legumineuses'])].index.to_list(),
                              'cereales':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['cereales'])].index.to_list(),
                              'legumes':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['legumes'])].index.to_list(),
                              'vop':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['vop'])].index.to_list(),
                              'fruits':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['fruits'])].index.to_list(),
                              'fruits_coque':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['fruits_coque'])].index.to_list(),
                              'prod_laitiers':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['prod_laitiers'])].index.to_list(),
                              'mat_grasses':self.df_ingredients_used[self.df_ingredients_used['alim_code'].isin(self.list_categories_alim_code['mat_grasses'])].index.to_list(),
                              'supplements':df_ingredients_used.index.to_list()}

        
        #dimensions des catégories d'ingrédients
        self.dim_feculents=len(self.list_categories_index['feculents'])
#        print('self.dim_feculents=',self.dim_feculents)
        self.dim_legumineuses=len(self.list_categories_index['legumineuses'])
#        print('self.dim_legumineuses=',self.dim_legumineuses)
        self.dim_cereales=len(self.list_categories_index['cereales'])      
#        print('self.dim_cereales=',self.dim_cereales)         
        self.dim_legumes=len(self.list_categories_index['legumes'])
        self.dim_vop=len(self.list_categories_index['vop'])
        self.dim_fruits=len(self.list_categories_index['fruits'])
        self.dim_fruits_coque=len(self.list_categories_index['fruits_coque'])
        self.dim_prod_laitiers=len(self.list_categories_index['prod_laitiers'])
        self.dim_mat_grasses=len(self.list_categories_index['mat_grasses'])
        self.dim_ingredients=len(self.list_categories_index['supplements'])
        
#        self.ingr_jour={'feculents':(self.nb_feculents_jour,self.dim_feculents),'legumineuses':(self.nb_legumineuses_jour,self.dim_legumineuses),
#                        'cereales':(self.nb_cereales_jour,self.dim_cereales),'legumes':(self.nb_legumes_jour,self.dim_legumes),
#                        'vop':(self.nb_vop_jour,self.dim_vop),'fruits':(self.nb_fruits_jour,self.dim_fruits),
#                        'fruits_coque':(self.nb_fruits_coque_jour,self.dim_fruits_coque),'prod_laitiers':(self.nb_prod_laitiers_jour,self.dim_prod_laitiers),
#                        'mat_grasses':(self.nb_mat_grasses_jour,self.dim_mat_grasses),'supplements':(self.nb_supplements_jour,self.dim_ingredients)}
        self.nb_ingr_jour={'feculents':self.nb_feculents_jour,'legumineuses':self.nb_legumineuses_jour,
                        'cereales':self.nb_cereales_jour,'legumes':self.nb_legumes_jour,
                        'vop':self.nb_vop_jour,'fruits':self.nb_fruits_jour,
                        'fruits_coque':self.nb_fruits_coque_jour,'prod_laitiers':self.nb_prod_laitiers_jour,
                        'mat_grasses':self.nb_mat_grasses_jour,'supplements':self.nb_supplements_jour}
        
        self.cat_list=list(self.nb_ingr_jour.keys())
        self.cat_cumsum=[0]+list(np.cumsum(list(self.nb_ingr_jour.values())))[:-1]
        self.cumsum_ingr_jour=dict(zip(self.cat_list,self.cat_cumsum))
        
#        self.comp={'pdej':{'feculents':0,'legumineuses':1,'cereales':1,'legumes':0,'vop':0,'fruits':2,'fruits_coque':1,'prod_laitiers':1, \
#                        'mat_grasses':1,'supplements':2}, \
#                    'dej':{'feculents':0,'legumineuses':1,'cereales':1,'legumes':3,'vop':1,'fruits':0,'fruits_coque':1,'prod_laitiers':0, \
#                        'mat_grasses':1,'supplements':2}, \
#                    'col':{'feculents':1,'legumineuses':0,'cereales':0,'legumes':1,'vop':0,'fruits':1,'fruits_coque':0,'prod_laitiers':1, \
#                        'mat_grasses':1,'supplements':2}, \
#                    'din':{'feculents':0,'legumineuses':1,'cereales':1,'legumes':3,'vop':1,'fruits':0,'fruits_coque':1,'prod_laitiers':0, \
#                        'mat_grasses':1,'supplements':2}}
#                    
#        list_indexes={'pdej':[],'dej':[],'col':[],'din':[]}
#        for rep in list(self.comp.keys()):
#            for key in list(self.comp[rep].keys()):
#                for k in range(self.comp[rep][key]):
#                    list_indexes_pdej.append(self.dim_ingredients*(self.cumsum_ingr_jour+k))
#                      
#        self.index_repas={'p_dej':list_indexes_pdej}
        
        
        #print(len(df_ingredients_used))
        # Observation = Etat des menus   
        #self.observation_space= spaces.Box(low=0, high=np.inf,shape=(self.dim_ingredients,self.nb_ingr_jour,self.nb_jours,self.nb_profils),dtype=np.float16)
        #low = np.array([0]*(self.dim_ingredients*self.nb_ingr_jour*self.nb_jours*self.nb_profils+self.nb_indic_holi+self.nb_nutriments))
        #self.dim_obs_space=(self.dim_ingredients*self.nb_tot_ingr_jour*self.nb_jours*self.nb_profils)#+self.nb_indic_holi*self.nb_profils+self.nb_nutriments*self.nb_profils)
        self.dim_obs_space=sum(self.nb_cat_repas['dej'].values())+self.nb_nutriments
        #sans les menus dans observation_space:
        #self.dim_obs_space=self.nb_nutriments
        
        low_obs=np.array([0]*self.dim_obs_space)
        high_obs=np.array([np.Inf]*self.dim_obs_space)
#        low_act=np.array(self.df_ingredients_used['poids_min'].to_list()*self.nb_tot_ingr_jour*self.nb_jours*self.nb_profils)
#        high_act = np.array(self.df_ingredients_used['poids_max'].to_list()*self.nb_tot_ingr_jour*self.nb_jours*self.nb_profils)
#        print("dim obs space=",self.dim_obs_space)
        self.observation_space = spaces.Box(low_obs, high_obs, dtype=np.float16) #-high instead of low=0 since some indic holi can be negative
        #espace des actions choisis pour réaliser les compositions spécifiques
#        self.action_space = spaces.Dict({'feculents':spaces.Tuple((spaces.MultiDiscrete((self.dim_feculents,1,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'legumineuses':spaces.Tuple((spaces.MultiDiscrete((self.dim_legumineuses,3,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'cereales':spaces.Tuple((spaces.MultiDiscrete((self.dim_cereales,3,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'prod_laitiers':spaces.Tuple((spaces.MultiDiscrete((self.dim_prod_laitiers,2,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'fruits':spaces.Tuple((spaces.MultiDiscrete((self.dim_fruits,3,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'fruits_coque':spaces.Tuple((spaces.MultiDiscrete((self.dim_fruits_coque,3,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'vop':spaces.Tuple((spaces.MultiDiscrete((self.dim_vop,2,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'legumes':spaces.Tuple((spaces.MultiDiscrete((self.dim_legumes,7,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'mat_grasses':spaces.Tuple((spaces.MultiDiscrete((self.dim_mat_grasses,4,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16))),
#                                        'supplements':spaces.Tuple((spaces.MultiDiscrete((self.dim_ingredients,8,self.nb_jours)),spaces.Box(low=0, high=np.inf,shape=(1,),dtype=np.float16)))})
        self.dim_act_space=sum(self.nb_cat_repas['dej'].values()) #self.dim_ingredients*self.nb_tot_ingr_jour*self.nb_jours*self.nb_profils
       #action_space= spaces.Box(low=0, high=np.inf,shape=(self.dim_ingredients,self.nb_ingr_jour,self.nb_jours,self.nb_profils),dtype=np.float16)
        #self.action_space= spaces.Box(low=0, high=np.inf,shape=(self.dim_act_space,),dtype=np.float16)
#        print("dim act space=",self.dim_act_space)
#        self.action_space=spaces.Box(low_act, high_act,dtype=np.float16)
        self.action_space=spaces.MultiDiscrete((self.dim_feculents,self.dim_legumineuses,self.dim_cereales,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_vop,
                                               self.dim_fruits_coque,self.dim_prod_laitiers,self.dim_mat_grasses))
        
#        self.observation_space = spaces.Dict({
#            "action_mask": spaces.Box(0, 1, shape=(self.dim_act_space,)),
#            "avail_actions": spaces.Box(low_act, high_act, dtype=np.float16),
#            "obs_space": spaces.Box(low_obs, high_obs, dtype=np.float16),
#        })
        
        
        
        
        
        
        
        
#        self.actual_action_space = spaces.Dict({'feculents':spaces.Tuple((spaces.MultiDiscrete((self.dim_feculents,self.dim_feculents,self.dim_feculents,self.dim_feculents,self.dim_feculents)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(5,),dtype=np.float16))),
#                                        'legumineuses':spaces.Tuple((spaces.MultiDiscrete((self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,
#                                                                                           self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,
#                                                                                           self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses,self.dim_legumineuses)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(15,),dtype=np.float16))),
#                                        'cereales':spaces.Tuple((spaces.MultiDiscrete((self.dim_cereales,self.dim_cereales,self.dim_cereales,self.dim_cereales,self.dim_cereales,
#                                                                                       self.dim_cereales,self.dim_cereales,self.dim_cereales,self.dim_cereales,self.dim_cereales,
#                                                                                       self.dim_cereales,self.dim_cereales,self.dim_cereales,self.dim_cereales,self.dim_cereales)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(15,),dtype=np.float16))),
#                                        'prod_laitiers':spaces.Tuple((spaces.MultiDiscrete((self.dim_prod_laitiers,self.dim_prod_laitiers,self.dim_prod_laitiers,self.dim_prod_laitiers,self.dim_prod_laitiers,
#                                                                                            self.dim_prod_laitiers,self.dim_prod_laitiers,self.dim_prod_laitiers,self.dim_prod_laitiers,self.dim_prod_laitiers)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(10,),dtype=np.float16))),
#                                        'fruits':spaces.Tuple((spaces.MultiDiscrete((self.dim_fruits,self.dim_fruits,self.dim_fruits,self.dim_fruits,self.dim_fruits,
#                                                                                     self.dim_fruits,self.dim_fruits,self.dim_fruits,self.dim_fruits,self.dim_fruits,
#                                                                                     self.dim_fruits,self.dim_fruits,self.dim_fruits,self.dim_fruits,self.dim_fruits)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(15,),dtype=np.float16))),
#                                        'fruits_coque':spaces.Tuple((spaces.MultiDiscrete((self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,
#                                                                                           self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,
#                                                                                           self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque,self.dim_fruits_coque)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(15,),dtype=np.float16))),
#                                        'vop':spaces.Tuple((spaces.MultiDiscrete((self.dim_vop,self.dim_vop,self.dim_vop,self.dim_vop,self.dim_vop,
#                                                                                  self.dim_vop,self.dim_vop,self.dim_vop,self.dim_vop,self.dim_vop)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(10,),dtype=np.float16))),
#                                        'legumes':spaces.Tuple((spaces.MultiDiscrete((self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,
#                                                                                      self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,
#                                                                                      self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,
#                                                                                      self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,
#                                                                                      self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,
#                                                                                      self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,
#                                                                                      self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_legumes)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(35,),dtype=np.float16))),
#                                        'mat_grasses':spaces.Tuple((spaces.MultiDiscrete((self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,
#                                                                                          self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,
#                                                                                          self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,
#                                                                                          self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses,self.dim_mat_grasses)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(20,),dtype=np.float16))),
#                                        'supplements':spaces.Tuple((spaces.MultiDiscrete((self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,
#                                                                                          self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,
#                                                                                          self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,
#                                                                                          self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,
#                                                                                          self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,
#                                                                                          self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,
#                                                                                          self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,
#                                                                                          self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients,self.dim_ingredients)),
#                                                                    spaces.Box(low=0, high=np.inf,shape=(40,),dtype=np.float16)))})
        
        self.state = self.observation_space.sample()
#        self.cat2dayind={'feculents':[19],'legumineuses':[0,9,26],'cereales':[1,10,27],'legumes':[12,13,14,21,29,30,31],'vop':[11,28],'fruits':[3,4,22]
#                        ,'fruits_coque':[5,15,32],'prod_laitiers':[2,20],'mat_grasses':[6,16,23,33],'supplements':[7,8,17,18,24,25,34,35]}
        self.nb_steps=0
        self.reward_taste=0
        self.reward_nutri=0
        self.reward=0
        self.seed()
#        
#        self.action_mask = np.array([0.] * self.dim_act_space)
#        for key in self.cat_list:
#            for i in range(self.nb_profils):
#                for j in range(self.nb_jours):
#                    for k in range(self.nb_ingr_jour[key]):
#                        #self.state[:,self.cat2dayind[key][j],i]=0
#                        indexes=np.array(self.list_categories_index[key]) \
#                        +(k+self.cumsum_ingr_jour[key])*self.dim_ingredients \
#                        +j*self.nb_tot_ingr_jour*self.dim_ingredients \
#                        +i*self.nb_jours*self.nb_tot_ingr_jour*self.dim_ingredients
#                        self.action_mask[indexes]=1
                        
    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
#    def update_avail_actions(self):
#        #on a juste besoin d'un mask qui determine les actions valides
#        self.action_mask = np.array([0.] * self.dim_act_space)
#        for key in self.cat_list:
#            for i in range(self.nb_profils):
#                for j in range(self.nb_jours):
#                    for k in range(self.nb_ingr_jour[key]):
#                        #self.state[:,self.cat2dayind[key][j],i]=0
#                        indexes=np.array(self.list_categories_index[key]) \
#                        +(k+self.cumsum_ingr_jour[key])*self.dim_ingredients \
#                        +j*self.nb_jours*self.nb_tot_ingr_jour*self.dim_ingredients \
#                        +i*self.nb_profils*self.nb_jours*self.nb_tot_ingr_jour*self.dim_ingredients
#                        self.action_mask[indexes]=1
        
   


    def step(self, action):
        #pass
        #must return observation(=state), reward, done, info       
        #actualise les ingredients
        self.state[:self.dim_act_space]=np.divide(action,np.array([self.dim_feculents,self.dim_legumineuses,self.dim_cereales,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_vop,
                                               self.dim_fruits_coque,self.dim_prod_laitiers,self.dim_mat_grasses])) #self.action_mask*action
        #print("action=",action)
        #actualise la nutrition
        list_indices=[self.list_categories_index[self.cat_dej[i]][action[i]] for i in range(self.dim_act_space)]
#        print("step list_indices=",list_indices)
        self.state[self.dim_act_space:self.dim_act_space+self.nb_nutriments]=self.df_ingredients_used.iloc[list_indices,self.en_col_ind:].sum(axis=0).to_list()
        
        self.nb_steps+=1
#        print(self.nb_steps)
        #taste reward       
        list_indices_ingrs=[self.list_categories_index[self.cat_dej[i]][action[i]] for i in range(self.dim_act_space)]
        list_alim_codes_ingrs=self.df_ingredients_used.loc[list_indices_ingrs,'alim_code']
        self.reward_taste=int(self.pairs_values[list(itertools.combinations(list_alim_codes_ingrs,2))].mean()/0.0050*100)
#        print("reward_taste=",reward_taste)
        #nutri reward
        self.reward_nutri=(self.state[self.dim_act_space:self.dim_act_space+self.nb_nutriments]>1).sum()*2 \
        +(self.state[self.dim_act_space:self.dim_act_space+self.nb_nutriments]<self.df_profil_ratio.iloc[0,:]).sum()
#        print("reward_nutri=",reward_nutri)
        #total reward
        self.reward=self.reward_nutri+self.reward_taste-self.nb_steps
        #determine when we achieved goal   
        done=False
        if (self.nb_steps>100 or ((self.state[self.dim_act_space:self.dim_act_space+self.nb_nutriments]>1).all() 
                                    and (self.state[self.dim_act_space:self.dim_act_space+self.nb_nutriments]<self.df_profil_ratio.iloc[0,:]).all())):
            done=True
            with open('log_recipes6.txt', 'a') as f:
                print("os.getpid()=",os.getpid(),"reward=",self.reward,", done=",done,file=f)#,", state=",self.state,", action=",action) 
                print("os.getpid()=",os.getpid(),"recette:",self.df_ingredients_used.loc[list_indices,['alim_nom_fr','poids']],file=f)
                for i in range(self.dim_act_space,self.dim_act_space+self.nb_nutriments):
                    print("os.getpid()=",os.getpid(),self.df_profil_ratio.columns[i-self.dim_act_space],self.state[i],file=f)
        return self.state, self.reward, done, {}
 
    def reset(self):
        self.nb_steps=0
        act=self.action_space.sample()
#        print("reset act=",act)
        
        self.state[:self.dim_act_space]=np.divide(act,np.array([self.dim_feculents,self.dim_legumineuses,self.dim_cereales,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_vop,
                                               self.dim_fruits_coque,self.dim_prod_laitiers,self.dim_mat_grasses])) #self.action_mask*self.action_space.sample()
        list_indices=[self.list_categories_index[self.cat_dej[i]][act[i]] for i in range(self.dim_act_space)]
        #print("reset list_indices=",list_indices)
        self.state[self.dim_act_space:self.dim_act_space+self.nb_nutriments]=self.df_ingredients_used.iloc[list_indices,self.en_col_ind:].sum(axis=0).to_list()
        return self.state
#        #pass
#        action=self.action_space.sample()
#        for i in range(self.nb_jours):
#            for key in ['cereales','feculents','fruits','fruits_coque','legumes','legumineuses','mat_grasses',
#                        'prod_laitiers','supplements','vop']:
#                nb_ingr_cat=int(len(action[key][0])/self.nb_jours)
#                for j in range(nb_ingr_cat):
#                    print(i,key,j,j+i*nb_ingr_cat)
#                    self.state[:,self.cat2dayind[key][j],i]=0
#                    self.state[self.list_categories_index[key][action[key][0][j+i*nb_ingr_cat]],self.cat2dayind[key][j],i]=action[key][1][j+i*nb_ingr_cat]
#        
#        #print(self.state.shape)
#        return self.state    
 
    def render(self, mode='human', close=False):
        #imprimer les recettes obtenues
        act=np.rint(self.state[:self.dim_act_space]*np.array([self.dim_feculents,self.dim_legumineuses,self.dim_cereales,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_vop,
                                               self.dim_fruits_coque,self.dim_prod_laitiers,self.dim_mat_grasses])).astype(int)
        list_indices=[self.list_categories_index[self.cat_dej[i]][act[i]] for i in range(self.dim_act_space)]
        print(self.df_ingredients_used.loc[list_indices,['alim_nom_fr','poids']])
        print("reward_nutri=",self.reward_nutri,"reward_taste=",self.reward_taste,"reward=",self.reward)
        print("nutrition:")
        #self.df_ingredients_used.iloc[:,self.en_col_ind:].columns
        result_nutri=pd.Series(self.state[self.dim_act_space:self.dim_act_space+self.nb_nutriments],
        index=['Energie, Règlement UE N° 1169/2011', 'Protéines',
               'Lipides', 'AG saturés',
               'Acide laurique + myristique + palmitique ', 'Cholestérol ',
               'AG 18:2 9c,12c (n-6), linoléique ',
               'AG 18:3 c9,c12,c15 (n-3), alpha-linolénique ', 'EPA + DHA ',
               'Glucides ', 'Sucres  totaux hors lactose ',
               'Fibres alimentaires ', 'Vitamine C ',
               'Vitamine B1 ou Thiamine ', 'Vitamine B2 ou Riboflavine ',
               'Vitamine B3 ou PP ou Niacine ',
               'Vitamine B5 ou Acide pantothénique ', 'Vitamine B6 ',
               'Vitamine B9 ou Folates totaux ', 'Vitamine B12 ',
               'Vitamine A  ', 'Vitamine E ', 'Sel chlorure de sodium ',
               'Calcium ', 'Phosphore ', 'Magnésium ', 'Fer ',
               'Iode ', 'Zinc ', 'Cuivre ', 'Sélénium ',
               'Manganèse '])
        print(result_nutri)
#        for icat,cat in enumerate(self.cat_dej):
#            self.list_categories_alim_code[cat][self.action[icat]]
#        
#        self.dim_legumineuses,self.dim_cereales,self.dim_legumes,self.dim_legumes,self.dim_legumes,self.dim_vop,
#                                               self.dim_fruits_coque,self.dim_mat_grasses,self.dim_ingredients,self.dim_ingredients
#        
#        self.nb_cat_repas={'pdej':{'feculents':0,'legumineuses':1,'cereales':1,'legumes':0,
#                        'vop':0,'fruits':2,'fruits_coque':1,'prod_laitiers':1,
#                        'mat_grasses':1,'supplements':2},
        
        
        pass

    def close(self):   
        self.close()
#        pass
#        if self.viewer:
#           self.viewer.close()
#           self.viewer = None
    
    def __del__(self):
        print ("object deleted")
        self.close()
        
    
