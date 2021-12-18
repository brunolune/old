#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:57:49 2019

@author: bruno
"""
import sys
import pandas as pd
import numpy as np
import itertools
import pickle
import gym

from util_bc import flat_list,flat_tuple,intersec
from rec_gram_mat_gen_bc import rec_gram_mat_gen
import re

#df_femmes = pd.read_excel("/home/bruno/wujico/data/femme.xlsm", sheet_name="data retravaillée")
df_ingredients= pd.read_excel("/home/bruno/wujico/data/table_ingredients2b.xls", sheet_name="Ingredients")
df_recettes=pd.read_pickle("/home/bruno/wujico/data/df_recettes_from_Edamam_labels_cuitsfirst.pkl")
df_recettes_fr = df_recettes[df_recettes.apply(lambda x:'French' in x['cuisineType'],axis=1)].reset_index()

#%% 
"""
to select specific dishTypes 
"""
dishtypes=set(flat_list(df_recettes_fr['dishType'].dropna().to_list()))
mains=['main course','meal','supper','lunch','dinner','Main Dish','Main Course','Other Main Dishes','brunch']
salads=['Lettuce Salads','Tuna Salads']
drinks=['Cocktails','digestif','high tea']
soups=['Chicken Soups','Chowders','Cream Soups','Fish Soups','Other Soups','Soups','Vegetable Soup']
breakfasts=['Eggs','Fruit Breakfast','Breakfast Casseroles','Meat Breakfast','Other Breakfast','breakfast','brunch']
appetizers=['Meat Appetizers','Other Appetizers','Other Breads','Seafood Appetizers','Vegetable Appetizers',
            'appetizer','buffet','nibble','starter','teatime','Cheese Appetizers','Appetizers']
desserts=['dessert','Cookies','Desserts','Fruit Desserts','Ice Cream & Ices','Other Desserts','Puddings']
heavys=['Beef','Beef Soups','Chicken','Chicken Soups','Chowders','Cream Soups','Fish','Fish Soups','Pies',
        'Pasta','Pork','Potatoes','Roasts','Savory Pies','Turkey']+mains+soups
lights=['Breads','Cakes','Cookies','Eggs','Other Breads','Other Side Dishes','Rice Sides','Side Dish',
        'Vegetables','buffet','nibble']+appetizers
ind_main_dish=[]
for irec,rec in df_recettes_fr.iterrows():
    if pd.isna(rec['dishType']):
        continue
    if (set(rec['dishType']).issubset(set(mains))):
        ind_main_dish.append(irec)
df_recettes_fr=df_recettes_fr.iloc[ind_main_dish,:].reset_index()
#df_recettes_fr[df_recettes_fr['dishType'].isin()]
#if(set(rec_gen).issubset(set(recipe))):

#%% 
"""
list_rec_ingr_gram: list of (ingredient,weight) for all recipes
df_rec_gram: dataframe of ingredients weights, rows=recipes, columns=ingredients

"""
list_rec_ingr_gram,df_rec_gram=rec_gram_mat_gen(df_recettes_fr, df_ingredients)

#%%
"""
determine most common ingredients
df_top_ingr: sorted list of ingredients by occurence in Edamam's recipes 

"""
df_rec_occ=df_rec_gram.where(df_rec_gram==0,1).astype(int)
top_ingr = list(df_rec_occ.sum(axis = 0))
df_top_ingr = pd.DataFrame({'alim_code': list(df_rec_occ.columns),'count': top_ingr})

df_top_ingr = df_top_ingr.join(df_ingredients[['alim_code','alim_nom_fr']].set_index('alim_code'), on='alim_code')
df_top_ingr = df_top_ingr.sort_values('count', ascending=False)

#%%
# Create correlation matrix
corr_matrix = df_rec_gram.corr()
# to generate a serie with pair of ingredients in indexes and correlation as their associated values
s = corr_matrix.unstack()
pairs_values= s.sort_values(kind="quicksort")
pairs_values= s[~pairs_values.isnull()]
tupleing = list(pairs_values.index)
pairs_values= pairs_values[[item for item in tupleing if item[0] != item[1]]]
#%%
"""
create df_common reduced dataframe of with useful columns of df_ingredients
fill col impt in df_ingredients to keep firstly the most common ingredients

"""
#change name df_common to df_recipes_ingredients
df_recipes_ingredients=df_ingredients.drop(['alim_nom_en_trad','search_keyword2'],axis=1).join(df_top_ingr[['alim_code','count']]
                            .set_index('alim_code'),on='alim_code').replace(np.NaN,0)
#define the importance of an ingredient
df_recipes_ingredients.loc[df_recipes_ingredients['count']>0,'importance']=True
df_recipes_ingredients.loc[df_recipes_ingredients['count']==0,'importance']=False
#change veg_or_an to vegetarian column
df_recipes_ingredients.rename(columns={"veg_or_an": "vegetarian"},inplace=True)
df_recipes_ingredients['vegetarian']=df_recipes_ingredients['vegetarian'].apply(lambda x: True if x==0 else False)
# Calculate the average weight of a specific ingredient in a recipe
df_mean_weight=df_rec_gram.replace(0,np.NaN).mean().reset_index()
df_std_weight=df_rec_gram.replace(0,np.NaN).std().reset_index()
df_mean_weight.columns=['alim_code','mean_weight']
df_std_weight.columns=['alim_code','std_weight']
# Zip the average weight with the most common ingredient
df_recipes_ingredients=df_recipes_ingredients.join(df_mean_weight.set_index('alim_code'),on='alim_code')
df_recipes_ingredients=df_recipes_ingredients.join(df_std_weight.set_index('alim_code'),on='alim_code')
df_recipes_ingredients['mean_weight'].replace(np.NaN,0,inplace=True)
df_recipes_ingredients['ingr_weight'] = list(zip(df_recipes_ingredients.alim_code, df_recipes_ingredients.mean_weight))
#keep only ingredients used in recipes
df_recipes_ingredients=df_recipes_ingredients[df_recipes_ingredients['count']>0]
#convert df_recipes_ingredients to numeric data only
df_recipes_ingredients.replace('traces',0.001,inplace=True) #replace traces by 0
df_recipes_ingredients.replace('-',0,inplace=True)
for label, content in df_recipes_ingredients.iloc[:,10:70].iteritems():
    df_recipes_ingredients[label]=df_recipes_ingredients[label].str.replace(',','.')
    df_recipes_ingredients[label]=df_recipes_ingredients[label].str.replace('<','') #replace <value by value
    df_recipes_ingredients.replace(np.nan,0,inplace=True)
    df_recipes_ingredients[label]=pd.to_numeric(df_recipes_ingredients[label])
#%%
##format with profiles columns
#df_recipes_ingredients_colprof=pd.DataFrame(columns=['alim_grp_nom_fr', 'alim_ssgrp_nom_fr', 'alim_ssssgrp_nom_fr',
#       'alim_code', 'alim_nom_fr', 'alim_nom_eng', 'keyword', 'modifiers',
#       'vegetarian', 'importance','Energie_(kcal)', 'Cholestérol', 'Sucres__totaux_hors_lactose',
#       'Fibres_(g)', 'Vitamine_C_(mg)*', 'Vitamine_B1_(mg)',
#       'Vitamine_B2_(mg)', 'Vitamine_B3_ou_PP_(mg)',
#       'Acide_Pantothénique_-_Vitamine_B5', 'Vitamine_B6_(mg)',
#       'Vitamine_B8_(µg)', 'Vitamine_12__(µg)', 'Vitamine_A__(µg)',
#       'Vitamine_E_(mg)', 'Sel_(g)', 'Calcium_(mg)_≤_24_ans',
#       'Calcium_(mg)_>_24_ans', 'Phosphore_(mg)', 'Fer_(mg)_Règles_faibles',
#       'Fer_(mg)_règles_abondantes', 'Iode_(µg)', 'Cuivre_(mg)',
#       'Sélénium__(µg)', 'Chrome_(µg)', 'Manganèse_(mg)', 'Fluor',
#       'Protéines', 'Lipides', 'AG_saturés',
#       'Acide_laurique_+_myristique_+_palmitique_)', 'Acide_linoléique',
#       'Acide_alpha_linolénique_ALA', 'EPA_+_DHA', 'Glucides__g',
#       'Vitamine_B9_(µg)', 'Magnésium_(mg)', 'Zinc_(mg)','poids_min','poids_max','count',
#       'mean_weight', 'std_weight', 'ingr_weight'])

df_recipes_ingredients_colprof=pd.DataFrame(columns=['alim_grp_nom_fr', 'alim_ssgrp_nom_fr', 'alim_ssssgrp_nom_fr',
       'alim_code', 'alim_nom_fr', 'alim_nom_eng', 'keyword', 'modifiers',
       'vegetarian', 'importance','count','poids_pdej','poids_dej','poids_col','Energie, Règlement UE N° 1169/2011 (kcal)',
       'Protéines (g)', 'Lipides (g)', 'AG saturés (g)',
       'Acide laurique + myristique + palmitique (g)', 'Cholestérol (mg)',
       'AG 18:2 9c,12c (n-6), linoléique (g)',
       'AG 18:3 c9,c12,c15 (n-3), alpha-linolénique (g)', 'EPA + DHA (mg)',
       'Glucides (g)', 'Sucres  totaux hors lactose (g)',
       'Fibres alimentaires (g)', 'Vitamine C (mg)',
       'Vitamine B1 ou Thiamine (mg)',
       'Vitamine B2 ou Riboflavine (mg)',
       'Vitamine B3 ou PP ou Niacine (mg)',
       'Vitamine B5 ou Acide pantothénique (mg)', 'Vitamine B6 (mg)',
       'Vitamine B9 ou Folates totaux (µg)',
       'Vitamine B12 (µg)', 'Vitamine A  (µg)', 'Vitamine E (mg)',
       'Sel chlorure de sodium (g)', 'Calcium (mg) ≤ 24 ans ',
       'Calcium (mg) > 24 ans', 'Phosphore (mg)', 'Magnésium (mg)',
       'Fer (mg) règles faibles ', 'Fer (mg) règles abondantes',
       'Iode (µg)', 'Zinc (mg)', 'Cuivre (mg)',
       'Sélénium (µg)', 'Manganèse (mg)'])
    
#copy first columns from df_recipes_ingredients    
df_recipes_ingredients_colprof['count']=df_recipes_ingredients['count']
df_recipes_ingredients_colprof['poids_pdej']=df_recipes_ingredients['poids_pdej']
df_recipes_ingredients_colprof['poids_dej']=df_recipes_ingredients['poids_dej']
df_recipes_ingredients_colprof['poids_col']=df_recipes_ingredients['poids_col']
df_recipes_ingredients_colprof.iloc[:,:10]=df_recipes_ingredients.iloc[:,:10].copy()
#assign nutriments columns
df_recipes_ingredients_colprof['Energie, Règlement UE N° 1169/2011 (kcal)']=df_recipes_ingredients['Energie, Règlement UE N° 1169/2011 (kcal/100g)']/100
df_recipes_ingredients_colprof['Cholestérol (mg)']=df_recipes_ingredients['Cholestérol (mg/100g)']/100
df_recipes_ingredients_colprof['Sucres  totaux hors lactose (g)']=df_recipes_ingredients['Sucres (g/100g)']/100+df_recipes_ingredients['Amidon (g/100g)']/100
df_recipes_ingredients_colprof['Fibres alimentaires (g)']=df_recipes_ingredients['Fibres alimentaires (g/100g)']/100
df_recipes_ingredients_colprof['Vitamine C (mg)']=df_recipes_ingredients['Vitamine C (mg/100g)']/100
df_recipes_ingredients_colprof['Vitamine B1 ou Thiamine (mg)']=df_recipes_ingredients['Vitamine B1 ou Thiamine (mg/100g)']/100
df_recipes_ingredients_colprof['Vitamine B2 ou Riboflavine (mg)']=df_recipes_ingredients['Vitamine B2 ou Riboflavine (mg/100g)']/100
df_recipes_ingredients_colprof['Vitamine B3 ou PP ou Niacine (mg)']=df_recipes_ingredients['Vitamine B3 ou PP ou Niacine (mg/100g)']/100
df_recipes_ingredients_colprof['Vitamine B5 ou Acide pantothénique (mg)']=df_recipes_ingredients['Vitamine B5 ou Acide pantothénique (mg/100g)']/100
df_recipes_ingredients_colprof['Vitamine B6 (mg)']=df_recipes_ingredients['Vitamine B6 (mg/100g)']/100
df_recipes_ingredients_colprof['Vitamine B12 (µg)']=df_recipes_ingredients['Vitamine B12 (µg/100g)']/100
df_recipes_ingredients_colprof['Vitamine A  (µg)']=df_recipes_ingredients['Rétinol (µg/100g)']/100+df_recipes_ingredients['Beta-Carotène (µg/100g)']/6/100
df_recipes_ingredients_colprof['Vitamine E (mg)']=df_recipes_ingredients['Vitamine E (mg/100g)']/100
df_recipes_ingredients_colprof['Sel chlorure de sodium (g)']=df_recipes_ingredients['Sel chlorure de sodium (g/100g)']/100
df_recipes_ingredients_colprof['Calcium (mg)']=df_recipes_ingredients['Calcium (mg/100g)']/100
df_recipes_ingredients_colprof['Phosphore (mg)']=df_recipes_ingredients['Phosphore (mg/100g)']/100
df_recipes_ingredients_colprof['Fer (mg)']=df_recipes_ingredients['Fer (mg/100g)']/100
df_recipes_ingredients_colprof['Iode (µg)']=df_recipes_ingredients['Iode (µg/100g)']/100
df_recipes_ingredients_colprof['Cuivre (mg)']=df_recipes_ingredients['Cuivre (mg/100g)']/100
df_recipes_ingredients_colprof['Sélénium (µg)']=df_recipes_ingredients['Sélénium (µg/100g)']/100
df_recipes_ingredients_colprof['Manganèse (mg)']=df_recipes_ingredients['Manganèse (mg/100g)']/100
df_recipes_ingredients_colprof['Protéines (g)']=df_recipes_ingredients['Protéines (g/100g)']/100
df_recipes_ingredients_colprof['Lipides (g)']=df_recipes_ingredients['Lipides (g/100g)']/100
df_recipes_ingredients_colprof['AG saturés (g)']=df_recipes_ingredients['AG saturés (g/100g)']/100
df_recipes_ingredients_colprof['Acide laurique + myristique + palmitique (g)']=df_recipes_ingredients['AG 12:0, laurique (g/100g)']/100+df_recipes_ingredients['AG 14:0, myristique (g/100g)']/100+df_recipes_ingredients['AG 16:0, palmitique (g/100g)']/100
df_recipes_ingredients_colprof['AG 18:2 9c,12c (n-6), linoléique (g)']=df_recipes_ingredients['AG 18:2 9c,12c (n-6), linoléique (g/100g)']/100
df_recipes_ingredients_colprof['AG 18:3 c9,c12,c15 (n-3), alpha-linolénique (g)']=df_recipes_ingredients['AG 18:3 c9,c12,c15 (n-3), alpha-linolénique (g/100g)']/100
df_recipes_ingredients_colprof['EPA + DHA (mg)']=(df_recipes_ingredients['AG 20:5 5c,8c,11c,14c,17c (n-3) EPA (g/100g)']/100+df_recipes_ingredients['AG 22:6 4c,7c,10c,13c,16c,19c (n-3) DHA (g/100g)']/100)
df_recipes_ingredients_colprof['Glucides (g)']=df_recipes_ingredients['Glucides (g/100g)']/100
df_recipes_ingredients_colprof['Vitamine B9 ou Folates totaux (µg)']=df_recipes_ingredients['Vitamine B9 ou Folates totaux (µg/100g)']/100
df_recipes_ingredients_colprof['Magnésium (mg)']=df_recipes_ingredients['Magnésium (mg/100g)']/100
df_recipes_ingredients_colprof['Zinc (mg)']=df_recipes_ingredients['Zinc (mg/100g)']/100   

#%% for 1st implementation of RL, need to change df_ingredients to create 3 ingredients per ingredients for weight min,mean,max
df_recipes_ingredients_colprof = df_recipes_ingredients_colprof.reset_index()
df_ingredients_used_3m=pd.DataFrame(0, index=np.arange(len(df_recipes_ingredients_colprof)*3)
                    ,columns=['index','alim_grp_nom_fr', 'alim_ssgrp_nom_fr', 'alim_ssssgrp_nom_fr',
                   'alim_code', 'alim_nom_fr', 'alim_nom_eng', 'keyword', 'modifiers',
                   'vegetarian', 'importance','count','poids',
                   'Energie, Règlement UE N° 1169/2011 (kcal)',
                   'Protéines (g)', 'Lipides (g)', 'AG saturés (g)',
                   'Acide laurique + myristique + palmitique (g)', 'Cholestérol (mg)',
                   'AG 18:2 9c,12c (n-6), linoléique (g)',
                   'AG 18:3 c9,c12,c15 (n-3), alpha-linolénique (g)', 'EPA + DHA (mg)',
                   'Glucides (g)', 'Sucres  totaux hors lactose (g)',
                   'Fibres alimentaires (g)', 'Vitamine C (mg)',
                   'Vitamine B1 ou Thiamine (mg)',
                   'Vitamine B2 ou Riboflavine (mg)',
                   'Vitamine B3 ou PP ou Niacine (mg)',
                   'Vitamine B5 ou Acide pantothénique (mg)', 'Vitamine B6 (mg)',
                   'Vitamine B9 ou Folates totaux (µg)',
                   'Vitamine B12 (µg)', 'Vitamine A  (µg)', 'Vitamine E (mg)',
                   'Sel chlorure de sodium (g)', 'Calcium (mg)',
                   'Phosphore (mg)', 'Magnésium (mg)',
                   'Fer (mg)','Iode (µg)', 'Zinc (mg)', 'Cuivre (mg)',
                   'Sélénium (µg)','Manganèse (mg)'])

#calcule de df_ingredients_used_3m
#to generate dej
for irow in range(df_recipes_ingredients_colprof.shape[0]):
    df_ingredients_used_3m.iloc[irow*3,0:12]=df_recipes_ingredients_colprof.iloc[irow,0:12].copy()
    df_ingredients_used_3m.loc[irow*3,'poids']=df_recipes_ingredients_colprof.loc[irow,'poids_dej']*0.8
    df_ingredients_used_3m.iloc[irow*3,13:]=df_recipes_ingredients_colprof.iloc[irow,15:]*df_recipes_ingredients_colprof.loc[irow,'poids_dej']*0.8
    
    df_ingredients_used_3m.iloc[irow*3+1,0:12]=df_recipes_ingredients_colprof.iloc[irow,0:12].copy()
    df_ingredients_used_3m.loc[irow*3+1,'poids']=df_recipes_ingredients_colprof.loc[irow,'poids_dej']
    df_ingredients_used_3m.iloc[irow*3+1,13:]=df_recipes_ingredients_colprof.iloc[irow,15:]*df_recipes_ingredients_colprof.loc[irow,'poids_dej']
    
    df_ingredients_used_3m.iloc[irow*3+2,0:12]=df_recipes_ingredients_colprof.iloc[irow,0:12].copy()
    df_ingredients_used_3m.loc[irow*3+2,'poids']=df_recipes_ingredients_colprof.loc[irow,'poids_dej']*1.2
    df_ingredients_used_3m.iloc[irow*3+2,13:]=df_recipes_ingredients_colprof.iloc[irow,15:]*df_recipes_ingredients_colprof.loc[irow,'poids_dej']*1.2
    
#%% pour débuguer
#listes catégories d'ingrédients
feculents=df_recipes_ingredients_colprof.groupby('alim_grp_nom_fr')['alim_code'].apply(list)['Féculents']
legumineuses=df_recipes_ingredients_colprof.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['légumineuses']
cereales=df_recipes_ingredients_colprof.groupby('alim_ssssgrp_nom_fr')['alim_code'].apply(list)['céréales'] \
                      +df_recipes_ingredients_colprof.groupby('alim_ssssgrp_nom_fr')['alim_code'].apply(list)['autres céréales']
legumes=df_recipes_ingredients_colprof.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['légumes']
vop=df_recipes_ingredients_colprof.groupby('alim_grp_nom_fr')['alim_code'].apply(list)['Viandes, œufs, poissons']
fruits=df_recipes_ingredients_colprof.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['fruits']
fruits_coque=df_recipes_ingredients_colprof.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['fruits à coque ']
prod_laitiers=df_recipes_ingredients_colprof.groupby('alim_grp_nom_fr')['alim_code'].apply(list)['Laits et produits laitiers']
mat_grasses=df_recipes_ingredients_colprof.groupby('alim_ssgrp_nom_fr')['alim_code'].apply(list)['végétale']
ingredients=df_recipes_ingredients_colprof['alim_code'].to_list()
#%%
##env=gym.make('AI_Recipes:AI_Recipes-v0',df_recipes_ingredients_colprof)
#from stable_baselines.common.policies import MlpPolicy
#from stable_baselines.common.vec_env import DummyVecEnv,SubprocVecEnv
#from stable_baselines import A2C
#from envs.AIRecipes_env import AIRecipesEnv
#
#env = DummyVecEnv([lambda: AIRecipesEnv(df_recipes_ingredients_colprof)])
##env= gym.make('AI_Recipes-v0',df_recipes_ingredients_colprof)
#state=np.squeeze(env.reset()) #comprends pas la dim suppl ajoutee, necessite squeeze
#
#env.close()
#%%
import argparse
import random
import numpy as np
import gym
from gym.spaces import Box, Discrete, Dict

import ray
from ray import tune
from ray.rllib.models import Model, ModelCatalog
from ray.rllib.models.misc import normc_initializer
from ray.tune.registry import register_env
from ray.rllib.utils import try_import_tf
from RLlibAIRecipes_env import RLlibAIRecipesEnv

#from ray.rllib.models.tf.tf_modelv2 import TFModelV2
#from ray.rllib.models.tf.fcnet_v2 import FullyConnectedNetwork

#tf = try_import_tf()
#
#df_femmes_bornes_min = pd.read_pickle("/home/bruno/wujico/data/df_femmes_bornes_min.pkl")
#df_femmes_bornes_max = pd.read_pickle("/home/bruno/wujico/data/df_femmes_bornes_max.pkl")
#df_femmes_bornes_min =df_femmes_bornes_min.drop(['Vitamine_B8_(µg)','Chrome_(µg)','Fluor'],axis=1)
#df_femmes_bornes_max =df_femmes_bornes_max.drop(['Vitamine_B8_(µg)','Chrome_(µg)','Fluor'],axis=1)
df_profilmax=pd.read_excel("/home/bruno/wujico/data/profil_max.xlsx")
df_profilmin=pd.read_excel("/home/bruno/wujico/data/profil_min.xlsx")
#droping profiles syndrome metabolique
df_profilmax.drop([29,30,31],axis=0,inplace=True)
df_profilmin.drop([29,30,31],axis=0,inplace=True)
#droping 'Vitamine_B8_(µg)','Chrome_(µg)','Fluor'
df_profilmax.drop(['nom','Fer (mg) règles faibles ', 'Fer (mg) règles abondantes',
                   'Calcium (mg) ≤ 24 ans ', 'Calcium (mg) > 24 ans','Vitamine B8 (µg)','Chrome (µg)','Fluor (g)'],axis=1,inplace=True)
df_profilmin.drop(['nom','Fer (mg) règles faibles ', 'Fer (mg) règles abondantes',
                   'Calcium (mg) ≤ 24 ans ', 'Calcium (mg) > 24 ans','Vitamine B8 (µg)','Chrome (µg)','Fluor (g)'],axis=1,inplace=True)
#normalisation par rapport à df_profilmin des colonnes de nutriments
#choix du profil
numprofil=0

df_ingredients_used_3m.iloc[:,13:]=df_ingredients_used_3m.iloc[:,13:].div(df_profilmin.iloc[numprofil,:],axis='columns') 
    
    
    #from RLlibAIRecipes_env import RLlibAIRecipesEnv
gym.envs.register(
     id='RLlibAIRecipes-v0',
     entry_point='gym.envs.classic_control:RLlibAIRecipesEnv',
     kwargs={'df_ingredients_used' : df_ingredients_used_3m,'df_profilmin':df_profilmin,
             'df_profilmax':df_profilmax,'pairs_values':pairs_values})
#env = gym.make('RLlibAIRecipes-v0')

#ray.init()
#register_env("RLlibAIRecipes-v0", lambda _: RLlibAIRecipesEnv)


## Instantiate and wrap the env
#env = DummyVecEnv([lambda: RLlibAIRecipesEnv(df_ingredients_used_3m,df_femmes_bornes_min,df_femmes_bornes_max)])
## Define and Train the agent
#model = A2C(MlpPolicy, env).learn(total_timesteps=100)
#%% utiliser RLLib?
#class MyKerasModel(TFModelV2):
#    """Custom model for policy gradient algorithms."""
#
#    def __init__(self, obs_space, action_space, num_outputs, model_config,
#                 name):
#        super(MyKerasModel, self).__init__(obs_space, action_space,
#                                           num_outputs, model_config, name)
#        self.inputs = tf.keras.layers.Input(
#            shape=obs_space.shape, name="observations")
#        layer_1 = tf.keras.layers.Dense(
#            128,
#            name="my_layer1",
#            activation=tf.nn.relu,
#            kernel_initializer=normc_initializer(1.0))(self.inputs)
#        layer_out = tf.keras.layers.Dense(
#            num_outputs,
#            name="my_out",
#            activation=None,
#            kernel_initializer=normc_initializer(0.01))(layer_1)
#        value_out = tf.keras.layers.Dense(
#            1,
#            name="value_out",
#            activation=None,
#            kernel_initializer=normc_initializer(0.01))(layer_1)
#        self.base_model = tf.keras.Model(self.inputs, [layer_out, value_out])
#        self.register_variables(self.base_model.variables)
#
#    def forward(self, input_dict, state, seq_lens):
#        model_out, self._value_out = self.base_model(input_dict["obs"])
#        return model_out, state
#
#    def value_function(self):
#        return tf.reshape(self._value_out, [-1])
#    
#
#ModelCatalog.register_custom_model("keras_model", MyKerasModel)
#tune.run(
#        'A3C',
#        stop={"timesteps_total": 10000},
#        config={
#            "env": "RLlibAIRecipes-v0",r(0.01))(layer_1)r(0.01))(layer_1)
#        self.base_model = tf.keras.Model(self.inputs, [layer_out, value_out])
#        self.register_variables(self.base_model.variables)
#
#    def forward(self, input_dict, state, seq_lens):
#        model_out, self._value_out = self.base_model(input_dict["obs"])
#        self.base_model = tf.keras.Model(self.inputs, [layer_out, value_out])
#        self.register_variables(self.base_model.variables)
#
#    def forward(self, input_dict, state, seq_lens):
#        model_out, self._value_out = self.base_model(input_dict["obs"])
#            "num_gpus": 0,
#            "num_workers":1
##            "model": {
##                "custom_model": "keras_model"
##            },
#        })
    
#%%   piece of code for openai gym multiagent algo

import os
from stable_baselines.common.vec_env import SubprocVecEnv,DummyVecEnv
from stable_baselines.common import set_global_seeds
from stable_baselines import A2C
from stable_baselines.common.policies import FeedForwardPolicy,MlpPolicy
import matplotlib.pyplot as plt
from stable_baselines.bench import Monitor
from stable_baselines.results_plotter import load_results, ts2xy
best_mean_reward, n_steps = 10, 0

def callback(_locals, _globals):
  """
  Callback called at each step (for DQN an others) or after n steps (see ACER or PPO2)
  :param _locals: (dict)
  :param _globals: (dict)

  """
  global n_steps, best_mean_reward
  # Print stats every 1000 calls
  if (n_steps + 1) % 100 == 0:
      # Evaluate policy training performance
      x, y = ts2xy(load_results(log_dir), 'timesteps')
      if len(x) > 0:
          mean_reward = np.mean(y[-100:])
          print(x[-1], 'timesteps')
          print("Best mean reward: {:.2f} - Last mean reward per episode: {:.2f}".format(best_mean_reward, mean_reward))

          # New best model, you could save the agent here
          if mean_reward > best_mean_reward:
              best_mean_reward = mean_reward
              # Example for saving best model
              print("Saving new best model")
              _locals['self'].save(log_dir + 'best_model.pkl')
  n_steps += 1
  return True

# Create log dir
log_dir = "/home/bruno/wujico/generators/RLlib_AI_Recipes/tmp/gym/"
os.makedirs(log_dir, exist_ok=True)



# Custom MLP policy of three layers of size 128 each
class CustomPolicy(FeedForwardPolicy):
    def __init__(self, *args, **kwargs):
        super(CustomPolicy, self).__init__(*args, **kwargs,
                                           net_arch=[128, dict(pi=[128, 128],
                                                          vf=[128, 128, 128])],
                                           feature_extraction="mlp")


def make_env(env_id, rank, seed=0):
    """
    Utility function for multiprocessed env.

    :param env_id: (str) the environment ID
    :param num_env: (int) the number of environments you wish to have in subprocesses
    :param seed: (int) the inital seed for RNG
    :param rank: (int) index of the subprocess
    """
    def _init():
        env = gym.make(env_id)
       # env = Monitor(env, log_dir, allow_early_resets=True) #wrap the environment?
        env.seed(seed + rank)
        return env
    set_global_seeds(seed)
    return _init

env_id = 'RLlibAIRecipes-v0'
num_cpu = 4  # Number of processes to use
# Create the vectorized multiprocessor environment
#env = SubprocVecEnv([lambda: gym.make('RLlibAIRecipes-v0') for i in range(num_cpu)])
#env = SubprocVecEnv([make_env(env_id, i) for i in range(num_cpu)])
#vectorized single cpu 
env = gym.make('RLlibAIRecipes-v0')
env = DummyVecEnv([lambda: env])



model = A2C(CustomPolicy, env,tensorboard_log='/home/bruno/wujico/generators/RLlib_AI_Recipes/tensorboardlog', verbose=1)
#model = A2C(MlpPolicy, env,tensorboard_log='/home/bruno/wujico/generators/RLlib_AI_Recipes/tensorboardlog', verbose=1)
model.learn(total_timesteps=1000)#, callback=callback)
model.save("a2c_RLlibAIRecipes_17_10_19")
#%%
from stable_baselines.common.vec_env import SubprocVecEnv,DummyVecEnv
from stable_baselines.common import set_global_seeds
from stable_baselines import A2C
from stable_baselines.common.policies import FeedForwardPolicy,MlpPolicy
env = gym.make('RLlibAIRecipes-v0')
env = DummyVecEnv([lambda: env])
#env = SubprocVecEnv([lambda: gym.make('RLlibAIRecipes-v0') for i in range(num_cpu)])
model = A2C.load("a2c_RLlibAIRecipes_17_10_19")
##
obs = env.reset()
for _ in range(2):
    action, _states = model.predict(obs)
    obs, rewards, dones, info = env.step(action)
    env.render()