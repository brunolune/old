#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:47:46 2019

@author: bruno
"""

import pandas as pd
import numpy as np
from util import rec_gram_mat_gen

pd.set_option("display.max_columns",7)

path_wujico = "/home/maxence/wujico"
#path_wujico = "/home/bruno/wujico"

df_ingredients = pd.read_excel(path_wujico+"/data/table_ingredients2.xls")
df_recettes= pd.read_pickle(path_wujico+"/data/df_recettes_from_Edamam_labels_cuitsfirst.pkl")
df_recettes_fr = df_recettes[df_recettes.apply(lambda x:'French' in x['cuisineType'],axis=1)].reset_index()

del df_recettes

profil_min = pd.read_excel(path_wujico+"/data/profil_min.xlsx")
profil_max = pd.read_excel(path_wujico+"/data/profil_max.xlsx")
profil_min.set_index('nom',inplace=True)
profil_max.set_index('nom',inplace=True)
nutriscore_tab_N=pd.read_excel(path_wujico+"/data/Nutriscore_tab.xlsx",sheet_name="Negatif")
nutriscore_tab_P=pd.read_excel(path_wujico+"/data/Nutriscore_tab.xlsx",sheet_name="Positif")

#%%

"""
list_rec_ingr_gram: list of (ingredient,weight) for all recipes
df_rec_gram: dataframe of ingredients weights, rows=recipes, columns=ingredients

"""
_,df_rec_gram=rec_gram_mat_gen(df_recettes_fr, df_ingredients)

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
# change veg_or_an to vegetarian column
df_recipes_ingredients.rename(columns={"veg_or_an": "vegetarian"},inplace=True)
df_recipes_ingredients['vegetarian']=df_recipes_ingredients['vegetarian'].apply(lambda x: True if x==0 else False)

df_recipes_ingredients=df_recipes_ingredients[df_recipes_ingredients['count']>0]

"""
Generate correlation matrix

"""

# Create correlation matrix
corr_matrix = df_rec_gram.corr()
# to generate a serie with pair of ingredients in indexes and correlation as their associated values
s = corr_matrix.unstack()
pairs_values= s.sort_values(kind="quicksort")
pairs_values= s[~pairs_values.isnull()]
tupleing = list(pairs_values.index)
pairs_values= pairs_values[[item for item in tupleing if item[0] != item[1]]]

#%%
# convert df_recipes_ingredients to numeric data only
columns_to_clean = df_recipes_ingredients.iloc[:,10:72].columns
df_recipes_ingredients[columns_to_clean] = df_recipes_ingredients[columns_to_clean].replace(['traces', '-'], ['0.0001', '0'])
for column in columns_to_clean:
    df_recipes_ingredients[column] = df_recipes_ingredients[column].astype(str).str.replace(',', '.')
    df_recipes_ingredients[column] = df_recipes_ingredients[column].astype(str).str.replace('<', '').astype(float)
df_recipes_ingredients = df_recipes_ingredients.drop(['poids_min', 'poids_max','count'], axis=1)

# Add new columns for corresponding 
df_recipes_ingredients['Vitamine A  (µg)'] = df_recipes_ingredients['Rétinol (µg/100g)'] + df_recipes_ingredients['Beta-Carotène (µg/100g)']/6
df_recipes_ingredients['Acide laurique + myristique + palmitique (g)'] = df_recipes_ingredients['AG 12:0, laurique (g/100g)'] + df_recipes_ingredients['AG 14:0, myristique (g/100g)'] + df_recipes_ingredients['AG 16:0, palmitique (g/100g)']
df_recipes_ingredients["EPA + DHA (mg)"] = df_recipes_ingredients['AG 20:5 5c,8c,11c,14c,17c (n-3) EPA (g/100g)'] + df_recipes_ingredients['AG 22:6 4c,7c,10c,13c,16c,19c (n-3) DHA (g/100g)']
df_recipes_ingredients["Sucres  totaux hors lactose (g)"] = df_recipes_ingredients['Amidon (g/100g)'] + df_recipes_ingredients['Sucres (g/100g)']

# Change unity from /100g to /1g 
df_recipes_ingredients.iloc[:,10:] = df_recipes_ingredients.iloc[:,10:]/100
df_recipes_ingredients.columns = [i.replace("/100g", "") for i in df_recipes_ingredients.columns]

# Will be used for normalization
columns_nutri = [col for col in profil_min.columns if col in df_recipes_ingredients.columns]
