#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:41:51 2019

@author: bruno
"""
import pandas as pd
import random
import numpy as np
import itertools
import heapq
from genetic_data import df_recipes_ingredients,columns_nutri,profil_min,profil_max,pairs_values,nutriscore_tab_N,nutriscore_tab_P
#%%
ssssgrp = df_recipes_ingredients.groupby('alim_ssssgrp_nom_fr')['alim_code'].apply(list).to_dict()
vp = ['abats','agneau et mouton','autres viandes', 'bœuf et veau', 'dinde, porc, poulet',
       'jambons cuits', 'jambons secs et crus','mollusques et crustacés crus','mollusques et crustacés cuits',
       'omelettes et autres ovoproduits','poissons crus','poissons cuits','saucisses et assimilés']
oeuf = ['œufs crus', 'œufs cuits']
feculent = ['autres céréales','céréales', 'farines', 'pommes de terre ', 'pâtes', 'riz', 'tubercules']
legumineuse = ['légumineuses cuites', 'légumineuses fraîches', 'légumineuses sèches']
fruit = ['fruits crus', 'fruits cuits']
legumes= ['légumes crus', 'légumes cuits']
fruitsec = ['fruits séchés']
graines = ['graines oléagineuses']
fruitcoq = ['fruits oléagineux']
laitage = ['autres laits','fromages blancs et yaourts', 'laits de vache ', 'boissons végétales']
fromage = ['fromages à pâte pressée', 'autres fromages']
matgras = ['beurre','crèmes et autres matières grasses','huiles de poissons','huiles et graisses végétales','margarines']
#%%
#tolerance proportion 20%
tol_prop=0.2

dict_ingredients= {
        "vp": {"alim_code":[], "gramptidej": 0, "gramdej": 100, "gramcol": 0, "percent":1},
        "oeuf": {"alim_code":[], "gramptidej": 50, "gramdej": 50, "gramcol": 0, "percent":1},
        "feculent": {"alim_code":[], "gramptidej": 150 , "gramdej": 200, "gramcol": 100, "percent":0.7},
        "legumineuse": {"alim_code":[], "gramptidej": 150 , "gramdej": 200, "gramcol": 100, "percent":0.3},
        "fruit et legume": {"alim_code":[], "gramptidej": 100, "gramdej": 250, "gramcol": 100, "percent":0.85},
        "fruit séché": {"alim_code":[], "gramptidej": 100, "gramdej": 250, "gramcol": 100, "percent":0.15},
        "graine": {"alim_code":[], "gramptidej": 20, "gramdej": 20, "gramcol": 20, "percent":0.5},
        "fruit à coque": {"alim_code":[], "gramptidej": 20, "gramdej": 20, "gramcol": 20, "percent":0.5},
        "laitage": {"alim_code":[], "gramptidej": 150, "gramdej": 150, "gramcol": 100, "percent":0.8},
        "fromage": {"alim_code":[], "gramptidej": 150, "gramdej": 150, "gramcol": 100, "percent":0.2},
        "matières grasses": {"alim_code":[], "gramptidej": 15, "gramdej": 15, "gramcol": 0, "percent":1},
        }

#fill dict_ingredients
for key,value in ssssgrp.items():
    if key in vp:
        dict_ingredients["vp"]["alim_code"].extend(value)
    if key in oeuf:
        dict_ingredients["oeuf"]["alim_code"].extend(value)
    elif key in feculent:
        dict_ingredients["feculent"]["alim_code"].extend(value)
    elif key in legumineuse:
        dict_ingredients["legumineuse"]["alim_code"].extend(value)
    elif key in fruitsec:
        dict_ingredients["fruit séché"]["alim_code"].extend(value)
    elif key in graines:
        dict_ingredients["graine"]["alim_code"].extend(value)
    elif key in fruitcoq:
        dict_ingredients["fruit à coque"]["alim_code"].extend(value)
    elif key in laitage:
        dict_ingredients["laitage"]["alim_code"].extend(value)
    elif key in fromage:
        dict_ingredients["fromage"]["alim_code"].extend(value)
    elif key in matgras:
        dict_ingredients["matières grasses"]["alim_code"].extend(value)
    if key in fruit or key in legumes:
        dict_ingredients["fruit et legume"]["alim_code"].extend(value)
        
list_family=[["vp"],["oeuf"],["feculent", 'legumineuse'],["fruit et legume", "fruit séché"],
             ["graine","fruit à coque"], ["laitage", 'fromage'], ["matières grasses"]]  
#%%
class ingredient:
    def __init__(self,alim_code=0,weight=0,alim_type=""):
        self.alim_code = alim_code
        self.weight = weight
        self.alim_type = alim_type
        self.keyword = df_recipes_ingredients[df_recipes_ingredients['alim_code']==self.alim_code]['keyword'].to_list()[0]
        
class repas():
    def __init__(self):
        self.ingrs = list()
        self.weightType = ""
        self.regle_compo = {
                "vp": {"nmin": 0, "nmax": 0},
                "oeuf": {"nmin": 0, "nmax": 1},
                "feculent": {"nmin": 1, "nmax": 2},
                "legumineuse": {"nmin": 1, "nmax": 2},
                "fruit et legume": {"nmin": 2, "nmax": 3},
                "fruit séché": {"nmin": 0, "nmax": 1},
                "graine": {"nmin": 0, "nmax": 2},
                "fruit à coque": {"nmin": 0, "nmax": 2},
                "laitage": {"nmin": 1, "nmax": 2},
                "fromage": {"nmin": 0, "nmax": 1},
                "matières grasses": {"nmin": 0, "nmax":1}
                            }
        self.composition = {}
        self.nutri_value = pd.Series()
        self.taste_value = 0
        self.familyStartPos = []        
    
    def setIngrNb(self,cat):      
        if cat == "legumineuse":
                self.composition[cat] = random.randint(
                    self.regle_compo[cat]["nmin"], 
                    self.regle_compo[cat]["nmax"]+1-self.composition["feculent"])
        elif cat == "fromage":
            self.composition[cat] = random.randint(
                    self.regle_compo[cat]["nmin"], 
                    self.regle_compo[cat]["nmax"]+1-self.composition["laitage"])  
        elif cat == "fruit à coque":
            #to force composition to have at least one graine oleagineux
            if self.composition["graine"]==0:
                self.composition[cat]=random.randint(1,2)
            else:
                self.composition[cat] = random.randint(
                        self.regle_compo[cat]["nmin"], 
                        self.regle_compo[cat]["nmax"]-self.composition["graine"])  
        else:
            self.composition[cat] = random.randint(
                    self.regle_compo[cat]["nmin"], 
                    self.regle_compo[cat]["nmax"])
            
    def getFamProp(self,fam):      
        if fam in [["vp"],["oeuf"],["matières grasses"]]: #case of families with single member
            return [[1]]
        elif self.composition[fam[0]]==0 or self.composition[fam[1]]==0: #case that one of the member is not present
            return [[1,1]]
        else:
            #fam[1] is the second member of a family that is always the smallest part in terms of weight
            p=dict_ingredients[fam[1]]["percent"]
            #it is determining the variance of dirichlet sharing by tol_prop % of the smallest component's proportion
            l=p*(1-p)/(tol_prop*p*2)**2 
            
        return np.random.dirichlet((l*(1-p),l*p),1)
    
    def genFamily(self,fam):
        list_ingrs=[]
        #determine proportion of ingredient in family
        p=self.getFamProp(fam)[0]
        # determine the ingredients
        for ifam,cat in enumerate(fam):
            list_prev_ingrs=[]
            for i in range(self.composition[cat]):
                while True:
                    alim_code = random.choice(dict_ingredients[cat]["alim_code"])
                    if not alim_code in list_prev_ingrs:
                        break
                alim_type = cat
                weight = round(dict_ingredients[cat][self.weightType]*p[ifam]/self.composition[cat])
                list_prev_ingrs.append(alim_code)
                list_ingrs.append(ingredient(alim_code,weight,alim_type))
        return list_ingrs
        
    
    def initRepas(self):        
        # determine the nb of ingredients in each families/categories
        for cat in self.regle_compo:
            self.setIngrNb(cat)
        #determine ingredients in each families/categories
        for fam in list_family:
            self.ingrs+=self.genFamily(fam)
        # set nutritional value of the recipe
        self.setNutriValue()
    
    def setNutriValue(self):     
#        np.sum([df_recipes_ingredients.loc[df_recipes_ingredients['alim_code']==i.alim_code][:,10:]*i.weight for i in self.ingrs])
        
        df_ingredients_nutri = df_recipes_ingredients.loc[df_recipes_ingredients['alim_code'].isin([i.alim_code for i in self.ingrs])]
        list_alim_code_df=df_ingredients_nutri['alim_code'].to_list()
        df_ingredients_nutri = df_ingredients_nutri.iloc[:,10:]
        list_alim_code_ingrs=[i.alim_code for i in self.ingrs]
        list_weight=[self.ingrs[list_alim_code_ingrs.index(x)].weight for x in list_alim_code_df]
#        print("reordered list ingrs=",[self.ingrs[list_alim_code_ingrs.index(x)].alim_code for x in list_alim_code_df])
        
#        self.nutri_value=df_ingredients_nutri.mul(list_weight, axis=0).sum()
#        print("list ingrs : ",[i.alim_code for i in self.ingrs])
        #print(self.familyStartPos)
        #print("df_ingredients_nutri : ",df_ingredients_nutri)
        self.nutri_value = np.sum(df_ingredients_nutri.apply(lambda w: np.asarray(w) * np.asarray(list_weight)))
#        return df_ingredients_nutri
         
    def getNutriValueNorm(self,group_number,profil_min_grouped):
        return self.nutri_value[columns_nutri].div(profil_min_grouped[columns_nutri].iloc[group_number,:]*self.apport_journalier)
    
    def getNutriGrade(self,profil_number):
        a=self.nutri_value[columns_nutri]>profil_min[columns_nutri].iloc[profil_number,:]*self.apport_journalier
        b=self.nutri_value[columns_nutri]<profil_max[columns_nutri].iloc[profil_number,:]*self.apport_journalier
        return (a & b).sum()/len(columns_nutri)
    
    def setTasteValue(self):
        self.taste_value = pairs_values[list(itertools.combinations([i.alim_code for i in self.ingrs],2))].mean()
    
    def getProfilScore(self,group_number,ratio_min_max,profil_min_grouped):
        nutri_value_norm = self.getNutriValueNorm(group_number,profil_min_grouped)
        #ratio_min_max = profil_max_grouped[columns_nutri].iloc[group_number,:].div(profil_min_grouped[columns_nutri].iloc[group_number,:])
        nutri_value_norm_mindist=nutri_value_norm - 1
        nutri_value_norm_maxdist=nutri_value_norm-ratio_min_max.iloc[group_number,:]
        return -(np.sum(-nutri_value_norm_mindist.where(nutri_value_norm_mindist <0,0))*2+np.sum(nutri_value_norm_maxdist.where(nutri_value_norm_maxdist>0,0)))
    
    def getTotalWeight(self):
        return np.sum(i.weight for i in self.ingrs)
    
    def getSAIN(self):
        totWeight=self.getTotalWeight()
        sain=(self.nutri_value['Protéines (g)']/0.65+self.nutri_value['Fibres alimentaires (g)']/0.3+self.nutri_value['Vitamine C (mg)']/1.1 \
            +self.nutri_value['Vitamine E (mg)']/0.12+self.nutri_value['Vitamine B1 ou Thiamine (mg)']/0.012+self.nutri_value['Vitamine B2 ou Riboflavine (mg)']/0.016 \
            +self.nutri_value['Vitamine B6 (mg)']/0.017+self.nutri_value['Vitamine B9 ou Folates totaux (µg)']/3.15+self.nutri_value['Calcium (mg)']/9 \
            +self.nutri_value['Fer (mg)']/0.125+self.nutri_value['Magnésium (mg)']/3.9+self.nutri_value['Zinc (mg)']/0.11 \
            +self.nutri_value['Potassium (mg)']/31+self.nutri_value['AG 18:3 c9,c12,c15 (n-3), alpha-linolénique (g)']/0.018+self.nutri_value['AG 22:6 4c,7c,10c,13c,16c,19c (n-3) DHA (g)']/0.0011) \
        /totWeight/0.15/self.nutri_value['Energie, Règlement UE N° 1169/2011 (kcal)']*100     
        return sain 
    
    def getLIM(self):
        list_sucres_ajoutes=[31087,31077,31067,31008,31089,31034,31016,31017]
        sucres_ajoutes=0
        for ingredient in self.ingrs:
            if ingredient.alim_code in list_sucres_ajoutes:
                sucres_ajoutes+= ingredient.weight*df_recipes_ingredients[df_recipes_ingredients['alim_code']==ingredient.alim_code]['Glucides (g)'].to_list()[0]    
        lim=(self.nutri_value['Sel chlorure de sodium (g)']/0.08+self.nutri_value['AG saturés (g)']/0.22+sucres_ajoutes/0.5)/self.getTotalWeight()/3*100
        return lim
    
    def getFF(self):
        calc_FF=41.7/(self.nutri_value['Energie, Règlement UE N° 1169/2011 (kcal)']/self.getTotalWeight()*100)**0.7 \
            +0.05*self.nutri_value['Protéines (g)']/self.getTotalWeight()*100 \
            +6.17E-4*(self.nutri_value['Fibres alimentaires (g)']/self.getTotalWeight()*100)**3 \
            -7.25E-6*(self.nutri_value['Lipides (g)']/self.getTotalWeight()*100)**3+0.617
        if calc_FF > 5:
            calc_FF=5
        if calc_FF < 0.5:
            calc_FF=0.5
        return calc_FF
    
    def getPRAL(self):
        pral=(0.49*self.nutri_value['Protéines (g)']+0.037*self.nutri_value['Phosphore (mg)'] \
                -0.021*self.nutri_value['Potassium (mg)']-0.026*self.nutri_value['Magnésium (mg)']-0.013*self.nutri_value['Calcium (mg)'])*100/self.getTotalWeight()
        return pral
    
    def getNutriScore(self):
        # 'A'=5,'B'=4,'C'=3,'D'=2,'E'=1
        dfN=self.nutri_value[['Energie, N x facteur Jones, avec fibres  (kJ)','AG saturés (g)','Sucres (g)','Sodium (mg)']]/self.getTotalWeight()*100
        dfN.index=['Energie, N x facteur Jones, avec fibres  (kJ/100g)','AG saturés (g/100g)','Sucres (g/100g)','Sodium (mg/100g)']
        dfP=self.nutri_value[['Fibres alimentaires (g)', 'Protéines (g)']]/self.getTotalWeight()*100
        dfP.index=['Fibres alimentaires (g/100g)', 'Protéines (g/100g)']
        dfP['Fruits et légumes, légumineuses et fruits à coque (g/100g) (%)'] \
        =np.sum([i.weight for i in self.ingrs if i.alim_type in ['legumineuse', 'fruit et legume','fruit à coque']])/self.getTotalWeight()*100     
        TabN=nutriscore_tab_N.iloc[:,[1,2,4,5]].lt(dfN,axis=1).sum()
        N=TabN.sum()
        TabP=nutriscore_tab_P.iloc[:,[1,2,3]].lt(dfP,axis=1).sum()
        P=TabP.sum()
        if (N>=11 and TabP['Fruits et légumes, légumineuses et fruits à coque (g/100g) (%)']<5):
            nutriscore = N - (TabP['Fibres alimentaires (g/100g)']+TabP['Fruits et légumes, légumineuses et fruits à coque (g/100g) (%)'])        
        else:
            nutriscore = N-P
        if nutriscore<0:
            return 5
        elif nutriscore <= 2:
            return 4
        elif nutriscore <= 10:
            return 3
        elif nutriscore <= 18:
            return 2
        else:
            return 1
         
    def setFamilyStartPos(self):
        self.familyStartPos = [0]
        for family in list_family[:-1]:
            self.familyStartPos.append(np.sum([self.composition[cat] for cat in family])+self.familyStartPos[-1])
                
    def updateCompo(self):
        for cat in self.composition.keys():
            self.composition[cat] = [ing.alim_type for ing in self.ingrs].count(cat)
            
    def updateValue(self):
        self.setNutriValue()
        self.setTasteValue()
                    
class ptitdej(repas):
    def __init__(self):
        super().__init__()
        self.weightType = "gramptidej"
        self.apport_journalier = 0.25
        
        self.initRepas()
        self.updateValue()
        self.setFamilyStartPos()
        
class dejdin(repas):
    def __init__(self):
        super().__init__()
        self.weightType = "gramdej"
        self.regle_compo["vp"]["nmin"] = 0
        self.regle_compo["vp"]["nmax"] = 1
        self.apport_journalier = 0.325
        
        self.initRepas()
        self.updateValue()
        self.setFamilyStartPos()

class col(repas):
    def __init__(self):
        super().__init__()
        self.weightType = "gramcol"
        self.regle_compo["oeuf"]["nmin"] = 0
        self.regle_compo["oeuf"]["nmax"] = 0
        self.regle_compo["matières grasses"]["nmin"] = 0
        self.regle_compo["matières grasses"]["nmax"] = 0
        self.apport_journalier = 0.1
        
        self.initRepas()
        self.updateValue()
        self.setFamilyStartPos()
        
class jour():
        
    def __init__(self):
        self.repas = [ptitdej(), dejdin(), col(), dejdin()]
        self.optim_infos ={}            
            
    def getScoreNutri(self,group_number,ratio_min_max,profil_min_grouped):
        return np.sum([i.getProfilScore(group_number,ratio_min_max,profil_min_grouped) for i in self.repas])
    
    def getScoreTaste(self):
        return np.sum([i.taste_value for i in self.repas])
    
    def getNutriNormGlobal(self,group_number,profil_min_grouped):
        df = self.repas[0].nutri_value+self.repas[1].nutri_value+self.repas[2].nutri_value+self.repas[3].nutri_value
        return df[columns_nutri].div(profil_min_grouped[columns_nutri].iloc[group_number,:])
     
    def getNoveltyScore(self,archive,notes_archive,group_number,profil_grouped_names,k):
        score_jour=[] #calculates all the scores
        keywords=[set([ingr.keyword for ingr in self.repas[i].ingrs]) for i in range(4)]
        for iday,day in enumerate(archive):
            score_repas=[]
            for i in range(4):   
                note_repas=notes_archive[group_number][iday][i]
                keywords_day_repas=set([ingr.keyword for ingr in day.repas[i].ingrs])
                inter=keywords[i].intersection(keywords_day_repas)
                if len(inter)==0:
                    score_repas.append(note_repas*(len(keywords)+len(keywords_day_repas)-0.5)/0.5)
                else:
                    score_repas.append(note_repas*(len(keywords)+len(keywords_day_repas)-len(inter))/len(inter))
#                if len(inter)>5:
#                    print("pid=",os.getpid(),"len(inter)=",len(inter))
#                    print("pid=",os.getpid(),"score_repas=",(len(keywords)+len(keywords_day_repas)-len(inter))/len(inter))
#                    print("pid=",os.getpid(),"keywords=",keywords[i])
#                    print("pid=",os.getpid(),"keywords_repas=",keywords_day_repas)
#                    print("pid=",os.getpid(),"note=",note_repas)
#                    print("pid=",os.getpid(),"score_repas/note=",note_repas*(len(keywords)+len(keywords_day_repas)-len(inter))/len(inter))
            
            score_jour.append(np.mean(score_repas))
        
        k_nearest=np.partition(score_jour,min(k,len(archive)))[:min(k,len(archive))]
        if k_nearest==[]:
            return 0
        else :
            return np.mean(k_nearest)
