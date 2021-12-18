# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:10:16 2019

@author: ericb & bruno & maxence
"""
import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import util
import copy
import time
import pickle
import os
import shutil
from multiprocessing import Process,Manager

from genetic_data import path_wujico,df_recipes_ingredients
from genetic_classes import list_family,jour
from genetic_profils import profil_min_grouped,profil_grouped_names,ratio_min_max, poids_demo

#%%
'''
All functions for genetic algo
'''
def calcul_probas(population,group_number,coeff_nutri,coeff_taste,coeff_mini,profils_common_ingredients,coeff_novelty,k):
    score_taste = np.array([day.getScoreTaste() for day in population])+1
    score_taste_norm=(score_taste-np.mean(score_taste))/np.std(score_taste)
#    print("score_taste:",score_taste)
    score_nutri = np.array([day.getScoreNutri(group_number,ratio_min_max,profil_min_grouped) for day in population])
    score_nutri_norm=(score_nutri-np.mean(score_nutri))/np.std(score_nutri)
#    print("score_nutri:",score_nutri)
    jour_all_ingrs=[set(util.flat_list([[i.alim_code for i in day.repas[j].ingrs] for j in range(4)])) for day in population]
#    print("jour_all_ingrs:", jour_all_ingrs)
    jour_all_ingrs_communs=[list(jour_all_ingrs[i].intersection(profils_common_ingredients)) for i in range(len(population))]
#    print("jour_all_ingrs_communs:", jour_all_ingrs_communs)
    score_mini= np.array([len(jour_all_ingrs_communs[i])/len(jour_all_ingrs[i]) for i in range(len(population))])
#    print("score_mini:",score_mini)
#    print("len:",[len(jour_all_ingrs[i]) for i in range(len(population))])
    score_mini_norm=(score_mini-np.mean(score_mini))/(np.std(score_mini)+1e-10)
#    print("score_mini:",score_mini)
    score_novelty=np.array([day.getNoveltyScore(archive_recettes,notes_archive,group_number,profil_grouped_names,k) for day in population])
#    print("score novelty=",score_novelty)
    score_novelty_norm=(score_novelty-np.mean(score_novelty))/(np.std(score_novelty)+1e-10)
#    print("score novelty norm=",score_novelty)
    exp_taste=np.exp(coeff_taste*score_taste_norm)
    exp_nutri=np.exp(coeff_nutri*score_nutri_norm)
    exp_mini=np.exp(coeff_mini*score_mini_norm)
    exp_novelty=np.exp(coeff_novelty*score_novelty_norm)
    exp_combined = exp_taste*exp_nutri*exp_mini*exp_novelty
    probabilities = exp_combined/np.sum(exp_combined)
    return probabilities,exp_combined,exp_taste,exp_nutri,exp_mini,score_mini,exp_novelty,score_novelty
    
def selection(population, probabilities):
    return np.random.choice(population,size=2,replace=False,p=probabilities)
  
def crossOver(parents,proba_crossover,proba_double_pivot):
    
    if proba_crossover>np.random.random():
        #to avoid affecting parents from population that may be selected again
        child1=copy.deepcopy(parents[0])
        child2=copy.deepcopy(parents[1])
        
        if 1-proba_double_pivot>np.random.random(): #simple pivot
            
            cut_jour_1 = random.randint(0,3)
            cut_repas_1 = random.randint(0,6)
        
            tmp=copy.deepcopy(child1.repas[cut_jour_1:])
        
            child1.repas[cut_jour_1:] = child2.repas[cut_jour_1:]
            child2.repas[cut_jour_1:]=copy.deepcopy(tmp)
            
            child2.repas[cut_jour_1].ingrs[0:child2.repas[cut_jour_1].familyStartPos[cut_repas_1]]=child1.repas[cut_jour_1].ingrs[0:child1.repas[cut_jour_1].familyStartPos[cut_repas_1]]
            child1.repas[cut_jour_1].ingrs[0:child1.repas[cut_jour_1].familyStartPos[cut_repas_1]]=tmp[0].ingrs[0:tmp[0].familyStartPos[cut_repas_1]]
        
        # Update repas after crossover
            child1.repas[cut_jour_1].updateCompo()
            child2.repas[cut_jour_1].updateCompo()
        
            child1.repas[cut_jour_1].setFamilyStartPos()
            child2.repas[cut_jour_1].setFamilyStartPos()
        
            child1.repas[cut_jour_1].updateValue()
            child2.repas[cut_jour_1].updateValue()   
        
        else: #double pivot
            cut_jour_1 = random.randint(0,3)
            if cut_jour_1==3:
                cut_repas_1 = random.randint(0,5)
            else:
                cut_repas_1 = random.randint(0,6)
            if (cut_repas_1==6):
                cut_jour_2=random.randint(cut_jour_1+1,3)
            else:
                cut_jour_2=random.randint(cut_jour_1,3)
            if cut_jour_1==cut_jour_2:
                cut_repas_2 = random.randint(cut_repas_1+1,6)
            else: 
                cut_repas_2 = random.randint(0,6)
                       
            tmp=copy.deepcopy(child1.repas[cut_jour_1:cut_jour_2+1])
            
            child1.repas[cut_jour_1:cut_jour_2+1]=child2.repas[cut_jour_1:cut_jour_2+1]
            child2.repas[cut_jour_1:cut_jour_2+1]=copy.deepcopy(tmp)
            
            
            child2.repas[cut_jour_1].ingrs[0:child2.repas[cut_jour_1].familyStartPos[cut_repas_1]]=child1.repas[cut_jour_1].ingrs[0:child1.repas[cut_jour_1].familyStartPos[cut_repas_1]]
            child1.repas[cut_jour_1].ingrs[0:child1.repas[cut_jour_1].familyStartPos[cut_repas_1]]=tmp[0].ingrs[0:tmp[0].familyStartPos[cut_repas_1]]
            
            child1.repas[cut_jour_1].updateCompo()
            child2.repas[cut_jour_1].updateCompo()
        
            child1.repas[cut_jour_1].setFamilyStartPos()
            child2.repas[cut_jour_1].setFamilyStartPos()
        
            child1.repas[cut_jour_1].updateValue()
            child2.repas[cut_jour_1].updateValue()  
        
        
            child2.repas[cut_jour_2].ingrs[child2.repas[cut_jour_2].familyStartPos[cut_repas_2]:]=child1.repas[cut_jour_2].ingrs[child1.repas[cut_jour_2].familyStartPos[cut_repas_2]:]
            child1.repas[cut_jour_2].ingrs[child1.repas[cut_jour_2].familyStartPos[cut_repas_2]:]=tmp[-1].ingrs[tmp[-1].familyStartPos[cut_repas_2]:]
            
            child1.repas[cut_jour_2].updateCompo()
            child2.repas[cut_jour_2].updateCompo()
                
            child1.repas[cut_jour_2].setFamilyStartPos()
            child2.repas[cut_jour_2].setFamilyStartPos()
                    
            child1.repas[cut_jour_2].updateValue()
            child2.repas[cut_jour_2].updateValue()   
                       
        return child1,child2
    else:
        return parents

def mutate(individu, max_nb_mutations,proba_mutate):
    
    if proba_mutate>np.random.random(): #proba_mutate*np.exp(-2*iteration/nb_iteration)
        nb_mutations = random.randint(1,max_nb_mutations)
    
        for i in range(nb_mutations):
#            print("new mutation:")
            repas_modif = random.randint(0,3)
            family_modif = random.randint(0,6)
#            print("new mutation: repas, famille =",repas_modif,family_modif)
            #remove previous ingredients of this family
            ind1=individu.repas[repas_modif].familyStartPos[family_modif]
#            ind2=ind1+individu.repas[repas_modif].composition[list_family[family_modif][0]]
            if family_modif!=6:
                ind2=individu.repas[repas_modif].familyStartPos[family_modif+1]
            else:
                ind2=len(individu.repas[repas_modif].ingrs)
#            print("start position : ",  individu.repas[repas_modif].familyStartPos)
#            print("ind1,ind2=",ind1,ind2)
            ingrs_to_remove=individu.repas[repas_modif].ingrs[ind1:ind2]
#            print("ingrs_to_remove",[i.alim_code for i in ingrs_to_remove])
#            print("ingrs_to_remove",[i.alim_type for i in ingrs_to_remove])
#            print("ingrs avant remove=",[i.alim_code for i in individu.repas[repas_modif].ingrs])
#            print("ingrs avant remove=",[i.alim_type for i in individu.repas[repas_modif].ingrs])
            for ing in ingrs_to_remove:
                individu.repas[repas_modif].ingrs.remove(ing)  
#            print("ingrs apres remove=",[i.alim_code for i in individu.repas[repas_modif].ingrs])
#            print("ingrs apres remove=",[i.alim_type for i in individu.repas[repas_modif].ingrs])
            #determine new ingredients numbers in family_modif categories
#            print("compo avant ", individu.repas[repas_modif].composition)
            for cat in list_family[family_modif]:
                individu.repas[repas_modif].setIngrNb(cat)
#            print("compo après ", individu.repas[repas_modif].composition)
            #generate new ingredients for the family_modif
            #individu.repas[repas_modif].ingrs.insert(ind1,individu.repas[repas_modif].genFamily(list_family[family_modif]))
            individu.repas[repas_modif].ingrs[ind1:ind1]=individu.repas[repas_modif].genFamily(list_family[family_modif])
#            print("ingrs apres insertion=",[i.alim_code for i in individu.repas[repas_modif].ingrs])
#            print("ingrs apres insertion=",[i.alim_type for i in individu.repas[repas_modif].ingrs])
            #recalculate crossOverPosition
            individu.repas[repas_modif].setFamilyStartPos()
            individu.repas[repas_modif].updateValue()
            
            return individu
    else: 
        return individu
#%%
'''
Visualisation functions for genetic algo
'''
def plot_convergence(titre,scores,iteration,nb_iteration,record_it,group_number):
    
    fig, ax = plt.subplots(figsize=(9.2, 5))
    if (int(iteration/record_it)-1)==0:
        x=np.arange(record_it+1)
#        print("iteration",iteration,"x",x)
#        print("len(x)",len(x))
#        print("len(scores)",len(scores))
    else:
        x=np.arange((int(iteration/record_it)-1)*record_it+1,int(iteration/record_it)*record_it+1)
#        print("iteration",iteration,"x",x)
#        print("len(x)",len(x))
#        print("len(scores)",len(scores))
    ax = sns.boxplot(x=x, y=scores,color="yellow")
#    ax = sns.swarmplot(x=x, y=scores, color=".25")
    ax.set_xlabel('Generation')
    ax.set_ylabel('Score')
    for ind, label in enumerate(ax.get_xticklabels()):
        if ind % int(iteration/10) == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)
#    ax.xaxis.set_major_locator(plt.MaxNLocator(10))
    if titre=='nutri':       
        ax.set_title('score nutritionnel (1=borne min reached)')
    if titre=='taste':
        ax.set_title('score de goût (somme des paires de correlations en ref. aux recettes du web)')
    plt.show()
    fig.tight_layout()
#    if iteration+1==nb_iteration:
#        fig.savefig('log_GA/'+d+'/plots/plot_scores_'+titre+'/plot_scores_'+titre+'_'+dt+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pm_'+
#                  str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+"_"+tagfile+".svg", format='svg')
#    else:
    fig.savefig('log_GA/'+d+'/plots/plot_scores_'+titre+'/plot_scores_'+titre+'_'+dt+'_gn_'+str(group_number)+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
              str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+"_"+tagfile+".png")
    plt.close(fig)
    
def plot_nutrition(df_nutri,recip,group_number):
        
    fig, ax = plt.subplots(figsize=(13, 10))
    
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(10, 150, n=10, as_cmap=True)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(df_nutri,vmin=0,vmax=1.,square=True,cmap=cmap, annot=True,
                fmt=".0f", linewidths=.5, cbar_kws={"shrink": .5})   
    plt.title("Intake with respect to recommended daily intake",fontsize=16)
    plt.tight_layout()
    fig.savefig("log_GA/"+d+"/plots/df_GA_best_nutri_intake_ratio/df_GA_"+recip+"_nutri_intake_ratio_"+dt+'_gn_'+str(group_number)+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pm_'+
              str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+"_"+tagfile+".png",dpi=150)
    plt.close(fig)
    
def plot_distribs(exp_combined,exp_nutri,exp_taste,exp_mini,exp_novelty,indsbest,group_number):
    fig, axes = plt.subplots(1, 5, figsize=(7, 7))
    
    plt.axes(axes[0])
    sns.distplot(exp_nutri, color="r",hist=False,rug=True,hist_kws={"histtype": "step"})
    for i, txt in enumerate(range(1,4)):
        plt.annotate(txt, (exp_nutri[indsbest[i]], 0))
#    plt.scatter(x=exp_nutri[indsbest[0]], y=0, color='k',marker="o")
#    plt.scatter(x=exp_nutri[indsbest[1]], y=0, color='k',marker="P")
#    plt.scatter(x=exp_nutri[indsbest[2]], y=0, color='k',marker="X")
        
    plt.axes(axes[1])
    sns.distplot(exp_taste, color="g",hist=False,rug=True,hist_kws={"histtype": "step"})
    for i, txt in enumerate(range(1,4)):
        plt.annotate(txt, (exp_taste[indsbest[i]], 0))
        
    plt.axes(axes[2])
    sns.distplot(exp_mini, color="m",hist=False,rug=True,hist_kws={"histtype": "step"})
    for i, txt in enumerate(range(1,4)):
        plt.annotate(txt, (exp_mini[indsbest[i]], 0))
        
    plt.axes(axes[3])
    sns.distplot(exp_novelty, color="k",hist=False,rug=True,hist_kws={"histtype": "step"})
    for i, txt in enumerate(range(1,4)):
        plt.annotate(txt, (exp_novelty[indsbest[i]], 0))
        
    plt.axes(axes[4])
    sns.distplot(exp_combined, color="b",hist=False,rug=True,hist_kws={"histtype": "step"},bins=10)
    for i, txt in enumerate(range(1,4)):
        plt.annotate(txt, (exp_combined[indsbest[i]], 0))

    axes[0].title.set_text('exp_nutri')
    axes[1].title.set_text('exp_taste')
    axes[2].title.set_text('exp_mini')
    axes[3].title.set_text('exp_novelty')
    axes[4].title.set_text('exp_combined')
    #plt.tight_layout()
    fig.savefig("log_GA/"+d+"/plots/plot_distribs/plot_distribs_"+dt+'_gn_'+str(group_number)+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
              str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+"_"+tagfile+".png",dpi=150)
    plt.close(fig)
#%%
def worker(iteration,group_number,pop_ingrs,profils_common_ingredients,
           all_scores_nutri,all_scores_taste):
#   fix seed for reproducibility
#   np.random.seed(123)
      
    with open('log_GA/'+d+'/tmp_populations/population_'+dt+'_gn_'+str(group_number)+'.pkl','rb') as f:
        population=pickle.load(f)
        
    new_population = []
    couples = []
    print("Début de la generation no ", iteration)
    start_it = time.time()
    probabilities,exp_combined,exp_taste,exp_nutri,exp_mini,_,_,_=calcul_probas(population,group_number,coeff_nutri,coeff_taste,coeff_mini,profils_common_ingredients,coeff_novelty,k)
    for new_individu in range(int(population_size/2)):
        parents = selection(population, probabilities)
        couples.append(parents)
    new_population = [list(crossOver(couple,proba_crossover,proba_double_pivot)) for couple in couples]
    new_population = [mutate(child,max_nb_mutations,proba_mutate) for child in util.flat_list(new_population)]
    end_it = time.time()                    
#   tps_calc.append(end_it - start_it)
    
    scores_nutri=[individu.getScoreNutri(group_number,ratio_min_max,profil_min_grouped) for individu in new_population]
    all_scores_nutri[group_number].append(scores_nutri)
    scores_taste=[individu.getScoreTaste() for individu in new_population]
    all_scores_taste[group_number].append(scores_taste)   
    
    #determine ingrs of population
    list_ingrs_repas_pop=[[ [i.alim_code for i in individu.repas[j].ingrs] for j in range(4)] for individu in new_population]
    set_all_ingrs=set(util.flat_list(util.flat_list(list_ingrs_repas_pop)))
    #save population ingrs for each profiles
    pop_ingrs[group_number]=set_all_ingrs
                
    #save new population
    with open('log_GA/'+d+'/tmp_populations/population_'+dt+'_gn_'+str(group_number)+'.pkl','wb') as f:
        pickle.dump(new_population,f)
        
    if iteration%record_it==0 and iteration!=0:   
        probabilities,exp_combined,exp_taste,exp_nutri,exp_mini,score_mini,exp_novelty,score_novelty=calcul_probas(new_population,group_number,coeff_nutri,coeff_taste,coeff_mini,profils_common_ingredients,coeff_novelty,k)
        PRALs=[np.mean([individu.repas[i].getPRAL() for i in range(4)]) for individu in new_population]
        SAINs=[np.mean([individu.repas[i].getSAIN() for i in range(4)]) for individu in new_population]
        LIMs=[np.mean([individu.repas[i].getLIM() for i in range(4)]) for individu in new_population]
        FFs=[np.mean([individu.repas[i].getFF() for i in range(4)]) for individu in new_population]
        NutriScores=[np.mean([individu.repas[i].getNutriScore() for i in range(4)]) for individu in new_population]
        meanPRAL=np.mean(PRALs)
        meanSAIN=np.mean(SAINs)
        meanLIM=np.mean(LIMs)
        meanFF=np.mean(FFs)
        meanNutriScore=np.mean(NutriScores)
        probassortedindices=np.argsort(exp_combined)[::-1]
        nutrisortedindices=np.argsort(scores_nutri)[::-1]
        tastesortedindices=np.argsort(scores_taste)[::-1]
        
        with open('log_GA/'+d+'/log_GA_'+dt+'_profil_'+str(group_number)+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
                  str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_'+tagfile+'.txt', 'a') as f:   
            print("------------------POPULATION----------------------------------------------",file=f)
            print("Generation no",iteration,":",file=f)
            print("Execution time: ",end_it - start_it, " secondes",file=f) 
            print("Mean score_nutri: ", np.mean(scores_nutri),file=f) 
            print("Mean score_taste: ", np.mean(scores_taste),file=f)
            print("Mean exp_combined =", np.mean(exp_combined),file=f)
            print("Mean exp_nutri =", np.mean(exp_nutri),file=f)
            print("Mean exp_taste =", np.mean(exp_taste),file=f)
            print("Mean score_mini =", np.mean(score_mini),file=f)
            print("Mean exp_novelty =", np.mean(exp_novelty),file=f)
            print("Mean score_novelty =", np.mean(score_novelty),file=f)
            print("Mean PRAL (~0?) : ", meanPRAL,file=f)
            print("Mean SAIN (>5?) : ", meanSAIN,file=f)
            print("Mean LIM (<7?) : ",meanLIM,file=f)
            print("Mean FF (>2.5?) : ",meanFF,file=f)
            print("Mean NutriScore (A=5,E=1) : ",meanNutriScore,file=f)
            print("Top 3 worst scores_nutri", np.array(scores_nutri)[nutrisortedindices[0:3]],file=f)
            print("Top 3 best scores_nutri", np.array(scores_nutri)[nutrisortedindices[-3:]][::-1],file=f)
            print("Top 3 best taste scores", np.array(scores_taste)[tastesortedindices[0:4]],file=f)
            print("Top 3 best exp_combined (combined scores)", exp_combined[probassortedindices[0:4]],file=f)
            print("Mean total weight", np.mean([np.sum([individu.repas[i].getTotalWeight() for i in range(4)]) for individu in new_population]),file=f)
            print("nb d'ingredients en commun entre profils (gen-1)=", len(profils_common_ingredients),file=f)
            #print concerned profils
            print("------------------PROFILS-------------------------------------------------",file=f)
            print ("Group's profil no",group_number,file=f)
            print ("profils:",file=f)
            print (profil_grouped_names.loc[group_number,'names'],file=f)
            #analyze diversity of population
            indfirstbest=probassortedindices[0]
            print("Nb of different ingredients in population = ",len(set_all_ingrs),file=f)
            alim_codes_best=[[i.alim_code for i in new_population[indfirstbest].repas[j].ingrs] for j in range(4)]
            ingrs_communs_first= np.array([[len([value for value in individu[i] if value in alim_codes_best[i]]) for i in range(4)] for individu in list_ingrs_repas_pop])  
            rangebests=max(10, int(len(new_population)/10)) #choose the best recipes in the 10% upper scores or at least 10 best scores for small populations
            print("Nb of best day 's ingredients in common with best scores days:",ingrs_communs_first[probassortedindices[0:rangebests]],file=f)         
            #Most different 2 others recipes among the 10 best scores:
            indsecondbest=probassortedindices[np.argsort(ingrs_communs_first[probassortedindices[0:rangebests]].sum(axis=1))[0]]
            alim_codes_second=[[i.alim_code for i in new_population[indsecondbest].repas[j].ingrs] for j in range(4)]
            ingrs_communs_second=np.array([[len([value for value in individu[i] if value in alim_codes_second[i]]) for i in range(4)] for individu in list_ingrs_repas_pop])
            ingrs_communs_first_second=ingrs_communs_first+ingrs_communs_second
            indthirdbest=probassortedindices[np.argsort(ingrs_communs_first_second[probassortedindices[0:rangebests]].sum(axis=1))[0]]
            indsbest=[indfirstbest,indsecondbest,indthirdbest]
            with open('log_GA/'+d+'/tmp_populations/indsbest_'+dt+'_gn_'+str(group_number)+'.pkl','wb') as fff:
                pickle.dump(indsbest,fff)
            #print recipes
            print("------------------BEST MOST DIFFERENT RECIPES--------------------------------------------",file=f)
            recip=["FIRST BEST","SECOND BEST","THIRD BEST"]
            for iind,ind in enumerate(indsbest):
                print("-------------------------"+recip[iind]+" RECIPES---------------------------------------------------",file=f)
                #print (recip[i]+" recipe (combined nutrition and taste):",file=f) 
                print ("position in combined score=",np.where(probassortedindices==ind)[0][0],file=f)
                print ("score_nutri =", new_population[ind].getScoreNutri(group_number,ratio_min_max,profil_min_grouped),file=f)
                print ("score_taste =", new_population[ind].getScoreTaste(),file=f)
                print ("exp_combined =", exp_combined[ind],file=f)
                print ("exp_nutri =", exp_nutri[ind],file=f)
                print ("exp_taste =", exp_taste[ind],file=f)
                print ("score_mini =", score_mini[ind],file=f)
                print ("score_novelty =", score_novelty[ind],file=f)
                print ("exp_novelty =", exp_novelty[ind],file=f)
                print ("PRAL = ",np.mean([new_population[ind].repas[i].getPRAL() for i in range(4)]),file=f)
                print ("SAIN = ",np.mean([new_population[ind].repas[i].getSAIN() for i in range(4)]),file=f)
                print ("LIM = ",np.mean([new_population[ind].repas[i].getLIM() for i in range(4)]),file=f)
                print ("FF = ",np.mean([new_population[ind].repas[i].getFF() for i in range(4)]),file=f)
                print ("NutriScore = ",np.mean([new_population[ind].repas[i].getNutriScore() for i in range(4)]),file=f)
                # print best recipes
                print("+++++++++++++++++++++++++++++++++++++++++++++",file=f)
                rep=["Petit Dejeuner:","Dejeuner:","Collation:","Diner:"]
                for j in range(4):
                    print(rep[j],file=f)
                    for ing in new_population[ind].repas[j].ingrs:    
                        print(ing.alim_type,":",df_recipes_ingredients[df_recipes_ingredients['alim_code']==ing.alim_code]['alim_nom_fr'].to_list()[0],", ",str(ing.weight),"g",file=f)
                    print("+++++++++++++++++++++++++++++++++++++++++++++",file=f)
                #print nutrition 
                df_nutri=pd.DataFrame([new_population[ind].repas[k].getNutriValueNorm(group_number,profil_min_grouped) for k in range(4)])
                df_nutri=df_nutri.transpose()
                df_nutri.columns=rep
                df_nutri['global']=new_population[ind].getNutriNormGlobal(group_number,profil_min_grouped)
                print(df_nutri,file=f)
                #save nutrition to excel file
                df_nutri.to_excel('log_GA/'+d+'/df_best_nutri/df_'+recip[iind]+'_nutri_'+dt+'_profil_'+str(group_number)+'_rec_'+str(iind)+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
                  str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+'_gn_'+str(group_number)+'_'+tagfile+'.xlsx')
#                print("df_nutri",df_nutri)
        if iteration==nb_iteration:
            with open('log_GA/'+d+'/tmp_populations/jour_'+dt+'_gn_'+str(group_number)+'.pkl','wb') as ff:
                pickle.dump(population[indfirstbest],ff)
            for ind in range(population_size):
                new_population[ind].optim_infos={'dt':dt,'pop':population_size,'pc':proba_crossover,'pdp':proba_double_pivot,'pm':proba_mutate,'cn':coeff_nutri,
                               'ct':coeff_taste,'cm':coeff_mini,'cnov':coeff_novelty,'ni':nb_iteration,'gn':group_number,'pn':profil_grouped_names.loc[group_number,'profils']}
            with open('log_GA/'+d+'/saved_populations/final_population_'+dt+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
                  str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_gn_'+str(group_number)+'_'+tagfile+'.pkl','wb') as ff:
                pickle.dump(new_population,ff)

if __name__ == "__main__":
    
    population = []   
    
    with open("archive_recettes.pkl","rb") as f:
        archive_recettes=pickle.load(f)
    notes_archive=[]
    for group_number in range(len(profil_grouped_names)):
        notes_archive.append([[np.mean([day.repas[i].getNutriGrade(j) for j in profil_grouped_names.loc[group_number,'profils']])
        for i in range(4)] for day in archive_recettes])

#    Parameters algo gen
    population_size = 50
    proba_crossover = 0.5
    proba_double_pivot = 0.75
    proba_mutate = 0.5   
    max_nb_mutations = 2
    nb_iteration = 50
    record_it=10
    #coefficients which set the importance of objective between nutrition, taste minimization and novelty
    coeff_nutri=2
    coeff_taste=1.5
    coeff_mini=1
    coeff_novelty=1
    k=4
    
    #file name
    tagfile='testplotconv'
    d=time.strftime("%Y_%m_%d", time.localtime())
    dt=time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
    
    #create a dir date if not exist, and subdirectories
    if not os.path.exists('log_GA/'+d):
        os.makedirs('log_GA/'+d)
        os.makedirs('log_GA/'+d+'/tmp_populations/')
        os.makedirs('log_GA/'+d+'/saved_progs')
        os.makedirs('log_GA/'+d+'/saved_populations')
        os.makedirs('log_GA/'+d+'/df_best_nutri')
        os.makedirs('log_GA/'+d+'/plots')
        os.makedirs('log_GA/'+d+'/plots/df_GA_best_nutri_intake_ratio')
        os.makedirs('log_GA/'+d+'/plots/plot_scores_nutri')
        os.makedirs('log_GA/'+d+'/plots/plot_scores_taste')
        os.makedirs('log_GA/'+d+'/plots/visu_recettes')
        os.makedirs('log_GA/'+d+'/plots/plot_distribs')
        
    #saving prog genetic_multiprocessing.py to be able to not loose track of how we obtained results
    src=path_wujico+'/generators/genrator_genetic/genetic_multiprocessing.py'
    dst=path_wujico+'/generators/genrator_genetic/log_GA/'+d+'/saved_progs/genetic_multiprocessing_'+dt+'_'+tagfile+'.py'
    shutil.copy2(src,dst)
    
    # Genesis: we use the average weight only for the genese
    start = time.time()
    print("Generation des individus à la gen 0")
    for i in range(population_size):
        print(i+1)
        population.append(jour())
        
    end = time.time()
    print("La genese a pris : ",end - start, " secondes")
#   start from given final population
#   with open('log_GA/2019_09_30/saved_populations/random_population_2019_09_30_13_25_07_pop_600_pc_0.95_pm_0.3_cn_0.6_ni_301_iter_-1_gn_0_testsavefiles.pkl','rb') as f:
#       population = pickle.load(f)
        
    #save random population
    iteration=-1
    with open('log_GA/'+d+'/saved_populations/random_population_'+dt+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
              str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+'_gn_-1_'+tagfile+'.pkl','wb') as f:
        pickle.dump(population,f)
    
    for group_number in range(len(profil_grouped_names)):
        with open('log_GA/'+d+'/tmp_populations/population_'+dt+'_gn_'+str(group_number)+'.pkl','wb') as f:
            pickle.dump(population,f)
    
    #determine ingrs of initial population
    list_ingrs_repas_pop=[[ [i.alim_code for i in individu.repas[j].ingrs] for j in range(4)] for individu in population]
    set_all_ingrs=set(util.flat_list(util.flat_list(list_ingrs_repas_pop)))
    
    #multiprocessing arrays (thread-safe)
    manager = Manager()
    pop_ingrs= manager.list([manager.list() for i in range(len(profil_grouped_names))])
    all_scores_nutri=manager.list([manager.list() for i in range(len(profil_grouped_names))])
    all_scores_taste=manager.list([manager.list() for i in range(len(profil_grouped_names))])
    #    tps_calc=manager.list([manager.list()])
    
    profils_common_ingredients=set_all_ingrs
    
    # Main loop
    for iteration in range(nb_iteration+1):
                
        #processes definition
        processes = [Process(name='profil '+str(j),target=worker,
                         args=(iteration,j,pop_ingrs,profils_common_ingredients,
                               all_scores_nutri,all_scores_taste)) for j in range(len(profil_grouped_names))]
        # Run processes
        for p in processes:
            p.start()
        # Exit the completed processes
        for p in processes:
            p.join()
        
        #determine ingredients common to all profiles
        profils_common_ingredients=pop_ingrs[0].intersection(*pop_ingrs[1:])

        
        #we plot graphs here because multiprocessing gives error if plotting in worker
        if iteration%record_it==0 and iteration!=0:    
            for group_number in range(len(profil_grouped_names)):
                #plot scores
                #print("all_scores_nutri[group_number]=",all_scores_nutri[group_number])
                scores_n=copy.deepcopy(all_scores_nutri[group_number])
                scores_t=copy.deepcopy(all_scores_taste[group_number])
                plot_convergence('nutri',scores_n,iteration,nb_iteration,record_it,group_number)
                plot_convergence('taste',scores_t,iteration,nb_iteration,record_it,group_number)
                #reset all_scores_nutri,all_scores_taste
#                all_scores_nutri[group_number]=[]
#                all_scores_taste[group_number]=[]
                #plot df_nutri
                recip=["FIRST BEST","SECOND BEST","THIRD BEST"]
                for irec in range(3):
                    df_nutri=pd.read_excel('log_GA/'+d+'/df_best_nutri/df_'+recip[irec]+'_nutri_'+dt+'_profil_'+str(group_number)+'_rec_'+str(irec)+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
                                           str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+'_gn_'+str(group_number)+'_'+tagfile+'.xlsx',index_col=0)
                    plot_nutrition(df_nutri,recip[irec],group_number) 
                #plot distribs
                with open('log_GA/'+d+'/tmp_populations/indsbest_'+dt+'_gn_'+str(group_number)+'.pkl','rb') as f:
                    indsbest=pickle.load(f)
                with open('log_GA/'+d+'/tmp_populations/population_'+dt+'_gn_'+str(group_number)+'.pkl','rb') as f:
                    population=pickle.load(f)
                probabilities,exp_combined,exp_taste,exp_nutri,exp_mini,score_mini,exp_novelty,score_novelty=calcul_probas(population,group_number,coeff_nutri,coeff_taste,coeff_mini,profils_common_ingredients,coeff_novelty,k)
                plot_distribs(exp_combined,exp_nutri,exp_taste,exp_mini,exp_novelty,indsbest,group_number) 
            all_scores_nutri=manager.list([manager.list() for i in range(len(profil_grouped_names))])
            all_scores_taste=manager.list([manager.list() for i in range(len(profil_grouped_names))])
              
    
with open('log_GA/'+d+'/saved_populations/profil_grouped_names_'+dt+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
          str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+'_'+tagfile+'.pkl','wb') as f:
    pickle.dump(profil_grouped_names,f)
with open('log_GA/'+d+'/saved_populations/poids_demo_'+dt+'_pop_'+str(population_size)+'_pc_'+str(proba_crossover)+'_pdp_'+str(proba_double_pivot)+'_pm_'+
          str(proba_mutate)+'_cn_'+str(coeff_nutri)+'_ct_'+str(coeff_taste)+'_cm_'+str(coeff_mini)+'_cnov_'+str(coeff_novelty)+'_ni_'+str(nb_iteration)+'_iter_'+str(iteration)+'_'+tagfile+'.pkl','wb') as f:
    pickle.dump(poids_demo,f)    
for group_number in range(len(profil_grouped_names)):
    with open('log_GA/'+d+'/tmp_populations/jour_'+dt+'_gn_'+str(group_number)+'.pkl','rb') as f:
         bestjour=pickle.load(f)
         bestjour.optim_infos={'dt':dt,'pop':population_size,'pc':proba_crossover,'pdp':proba_double_pivot,'pm':proba_mutate,'cn':coeff_nutri,
                               'ct':coeff_taste,'cm':coeff_mini,'cnov':coeff_novelty,'ni':nb_iteration,'gn':group_number,'pn':profil_grouped_names.loc[group_number,'profils']}
         archive_recettes.append(bestjour)
with open("archive_recettes.pkl","wb") as f:
    pickle.dump(archive_recettes,f)
with open("log_GA/"+d+"/archive_recettes_"+dt+".pkl","wb") as f:
    pickle.dump(archive_recettes,f)
#print("moyenne tps de calcul",np.mean(tps_calc))

