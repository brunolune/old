#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 17:49:08 2019

@author: bruno
"""
from util import flat_list
import numpy as np
import pandas as pd

def rec_gram_mat_gen(df_recettes, df_ingredients):
    """
    list_rec_ingr_gram: list of (ingredient,weight) for all recipes
    df_rec_gram: dataframe of ingredients weights, rows=recipes, columns=ingredients

    """
    list_rec_ingr_gram=[]
    rec_gram_mat=np.zeros((df_recettes.shape[0],df_ingredients.shape[0]+1))
    for irec,rec in df_recettes.iterrows():
        list_rec_ingr_gram.append([])
        for ingr in flat_list(rec['ingrMatch']):
            if isinstance(ingr[0], float):
                list_rec_ingr_gram[irec].append((np.nan,ingr[2]))
                rec_gram_mat[irec,-1]=ingr[2]
            else:
                list_rec_ingr_gram[irec].append((ingr[0],ingr[2]))
                rec_gram_mat[irec,list(df_ingredients['alim_code']).index(ingr[0])]=ingr[2] 
    df_rec_gram=pd.DataFrame(rec_gram_mat,columns=list(df_ingredients['alim_code'])+[-1])
    df_rec_gram=df_rec_gram.replace(0,np.nan).dropna(how='all',axis=1)
    df_rec_gram=df_rec_gram.replace(np.nan,0)
            
    return list_rec_ingr_gram,df_rec_gram