# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import sys
sys.path.append('../src')
import Nuesslein1
import Nuesslein2
import CJ1
import CJ2
import utils
from pathlib import Path
import numpy as np
import dimod
from dimod import BinaryQuadraticModel
import dwave.embedding
import dwave_networkx as dnx
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
from minorminer import find_embedding
import matplotlib.pyplot as plt
import random
import pandas as pd

# %% [markdown]
# ### Comparison of all results

# %%
files = ['EmbeddingComposite', 'QUBO_scaled', 'bqm_emb_scale']

gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [5, 12, 50] #[5,10,12,20,50]
num_instances = 1
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in range(num_instances):
            for file_names in files:
                #Instance for getting later the clauses
                file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
                instance = getattr(g_module, g)(file_path)
                if file_names == 'bqm_emb_scale':
                    instance.fillQ()
                    Q = instance.Q
                    couplings=[]
                    for k in Q.keys():
                        if k[0]!=k[1]:
                            couplings.append(k)
                    embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            
                sample_set_str=''
                with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_{file_names}.txt', 'r') as f:
                    lines = f.readlines()
                    if file_names == 'EmbeddingComposite':
                        k=3
                    if file_names == 'QUBO_scaled' or file_names == 'bqm_emb_scale':
                        k=4
                    for r in range(k,len(lines)-1):
                        if r==k:
                            sample_set_str+=lines[r][55:]
                        elif r==len(lines)-2:
                            line=lines[r].replace(" ", "")
                            sample_set_str+=line[:-2]
                        else:
                            sample_set_str+=lines[r]
                    sample_set=eval(sample_set_str)
                    variables = [v for v in eval(lines[k-1][2:]).keys()]
                    scale_factor = eval(lines[0].split()[-1])

                o_distr = {}
                for sample in sample_set:
                    if file_names == 'bqm_emb_scale':
                        #This is done for also SPIN sample_sets
                        majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==1:
                                        majority_vote[var][0]+=1
                                    if value!=1:
                                        majority_vote[var][1]+=1
                                    if majority_vote[var][0]>majority_vote[var][2]/2:
                                        assignment[var] = 1
                                        break
                                    if majority_vote[var][1]>majority_vote[var][2]/2:
                                        assignment[var] = 0
                                        break
                                    if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                                        assignment[var] = random.choice([0,1])
                                        break

                    else:
                        assignment={variables[j]:sample[0][j] for j in range(len(sample[0]))}

                    
                    o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
                    if o not in o_distr.keys():
                        o_distr[o] = {sample[-3]: sample[-2]}
                    else:
                        if sample[-3] not in o_distr[o].keys():
                            o_distr[o][sample[-3]] = sample[-2]
                        else:
                            o_distr[o][sample[-3]] += sample[-2]
                o_distr = dict(sorted(o_distr.items()))
                print('----------------------------------')
                print('Gadget: ', g, ' Vars: ', vars, ' Scaling_factor: ', scale_factor, ' ',file_names)
                print('Energy distribution for each optimum:\n', o_distr)
                o_distr_pond_mean = {}
                for o, distr in o_distr.items():
                    sols = sum(value for value in distr.values())
                    o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
                o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
                print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)
        print('==================================')
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

# %% [markdown]
# ### x: Suboptims ; y: Energies ; color: ocurrence ; Legend: different aproaches for solving

# %%
marker=['o', 'v', 'P']

files = ['bqm_emb_scale'] #['EmbeddingComposite', 'QUBO_scaled', 'bqm_emb_scale']

gadgets = {
    #'Nuesslein1': Nuesslein1,
    #'Nuesslein2': Nuesslein2,
    #'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [12] #[5, 12, 50] #[5,10,12,20,50]
num_instances = 1
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in range(num_instances):
            
            plt.figure(figsize=(6, 4))
            m=0
            
            for file_names in files:
                #Instance for getting later the clauses
                file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
                instance = getattr(g_module, g)(file_path)
                if file_names == 'bqm_emb_scale':
                    instance.fillQ()
                    Q = instance.Q
                    couplings=[]
                    for k in Q.keys():
                        if k[0]!=k[1]:
                            couplings.append(k)
                    embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            
                sample_set_str=''
                with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_{file_names}.txt', 'r') as f:
                    lines = f.readlines()
                    if file_names == 'EmbeddingComposite':
                        k=3
                    if file_names == 'QUBO_scaled' or file_names == 'bqm_emb_scale':
                        k=4
                    for r in range(k,len(lines)-1):
                        if r==k:
                            sample_set_str+=lines[r][55:]
                        elif r==len(lines)-2:
                            line=lines[r].replace(" ", "")
                            sample_set_str+=line[:-2]
                        else:
                            sample_set_str+=lines[r]
                    sample_set=eval(sample_set_str)
                    variables = [v for v in eval(lines[k-1][2:]).keys()]
                    scale_factor = eval(lines[0].split()[-1])

                o_distr = {}
                for sample in sample_set:
                    if file_names == 'bqm_emb_scale':
                        #This is done for also SPIN sample_sets
                        majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==1:
                                        majority_vote[var][0]+=1
                                    if value!=1:
                                        majority_vote[var][1]+=1
                                    if majority_vote[var][0]>majority_vote[var][2]/2:
                                        assignment[var] = 1
                                        break
                                    if majority_vote[var][1]>majority_vote[var][2]/2:
                                        assignment[var] = 0
                                        break
                                    if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                                        assignment[var] = random.choice([0,1])
                                        break

                    else:
                        assignment={variables[j]:sample[0][j] for j in range(len(sample[0]))}

                    
                    o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
                    if o not in o_distr.keys():
                        o_distr[o] = {sample[-3]: sample[-2]}
                    else:
                        if sample[-3] not in o_distr[o].keys():
                            o_distr[o][sample[-3]] = sample[-2]
                        else:
                            o_distr[o][sample[-3]] += sample[-2]
                o_distr = dict(sorted(o_distr.items()))

                print('----------------------------------')
                print('Gadget: ', g, ' Vars: ', vars, ' Scaling_factor: ', scale_factor, ' ',file_names) 
                o_distr_pond_mean = {}
                for o, distr in o_distr.items():
                    sols = sum(value for value in distr.values())
                    o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
                o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
                print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)

                #Graficar
                optimos = list(o_distr.keys())
                if file_names != 'EmbeddingComposite':
                    energias = [[o_d*scale_factor for o_d in o_distr[opt].keys()] for opt in optimos]
                else:
                    energias = [[o_d for o_d in o_distr[opt].keys()] for opt in optimos]
                
                ocurrencias = [[o_d for o_d in o_distr[opt].values()] for opt in optimos]                
                
                print(optimos)
                # Dibujar los puntos
                datax = []
                datay = []
                dataoc = []
                for x, opt in enumerate(optimos):
                    for y, energia in enumerate(energias[x]):
                        datax.append(opt)
                        datay.append(energia)
                        dataoc.append(ocurrencias[x][y])
                        #plt.scatter(opt, energia, s=ocurrencias[x][y]*10, c=ocurrencias[x][y], cmap='viridis', alpha=0.6)
                plt.scatter(datax, datay, c=dataoc, cmap='viridis', alpha=0.6, label=file_names, marker=marker[m], vmin=0, vmax=100)
                m+=1
            # Configurar los ejes
            plt.xticks(range(5))
            #plt.yticks(range(len(energias)))
            plt.xlabel('Optimos')
            plt.ylabel('Energias')
            plt.title(f'Solutions for the Max3SAT problem\n{vars} vars with {g}')
            
            # Añadir una barra de colores
            cbar = plt.colorbar()
            cbar.set_label('Ocurrencia')
                
            plt.tight_layout()
            plt.legend()
            plt.show()
        print('==================================')
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

# %% [markdown]
# ### x: Suboptims ; y: Energies ; color: ocurrences ; legend: Gadgets

# %%
marker=['o', 'v', 'P']

files = ['QUBO_scaled'] #['EmbeddingComposite', 'QUBO_scaled', 'bqm_emb_scale'] #JUST 1!!!!

gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    #'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [50] #[5, 12, 50] #[5,10,12,20,50]
num_instances = 1

plt.figure(figsize=(6, 4))
m=0
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in range(num_instances):            
            for file_names in files:
                #Instance for getting later the clauses
                file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
                instance = getattr(g_module, g)(file_path)
                if file_names == 'bqm_emb_scale':
                    instance.fillQ()
                    Q = instance.Q
                    couplings=[]
                    for k in Q.keys():
                        if k[0]!=k[1]:
                            couplings.append(k)
                    embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            
                sample_set_str=''
                with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_{file_names}.txt', 'r') as f:
                    lines = f.readlines()
                    if file_names == 'EmbeddingComposite':
                        k=3
                    if file_names == 'QUBO_scaled' or file_names == 'bqm_emb_scale':
                        k=4
                    for r in range(k,len(lines)-1):
                        if r==k:
                            sample_set_str+=lines[r][55:]
                        elif r==len(lines)-2:
                            line=lines[r].replace(" ", "")
                            sample_set_str+=line[:-2]
                        else:
                            sample_set_str+=lines[r]
                    sample_set=eval(sample_set_str)
                    variables = [v for v in eval(lines[k-1][2:]).keys()]
                    scale_factor = eval(lines[0].split()[-1])

                o_distr = {}
                for sample in sample_set:
                    if file_names == 'bqm_emb_scale':
                        #This is done for also SPIN sample_sets
                        majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==1:
                                        majority_vote[var][0]+=1
                                    if value!=1:
                                        majority_vote[var][1]+=1
                                    if majority_vote[var][0]>majority_vote[var][2]/2:
                                        assignment[var] = 1
                                        break
                                    if majority_vote[var][1]>majority_vote[var][2]/2:
                                        assignment[var] = 0
                                        break
                                    if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                                        assignment[var] = random.choice([0,1])
                                        break

                    else:
                        assignment={variables[j]:sample[0][j] for j in range(len(sample[0]))}

                    
                    o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
                    if o not in o_distr.keys():
                        o_distr[o] = {sample[-3]: sample[-2]}
                    else:
                        if sample[-3] not in o_distr[o].keys():
                            o_distr[o][sample[-3]] = sample[-2]
                        else:
                            o_distr[o][sample[-3]] += sample[-2]
                o_distr = dict(sorted(o_distr.items()))

                print('----------------------------------')
                print('Gadget: ', g, ' Vars: ', vars, ' Scaling_factor: ', scale_factor, ' ',file_names) 
                o_distr_pond_mean = {}
                for o, distr in o_distr.items():
                    sols = sum(value for value in distr.values())
                    o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
                o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
                print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)

                #Graficar
                optimos = list(o_distr.keys())
                if file_names != 'EmbeddingComposite':
                    energias = [[o_d*scale_factor for o_d in o_distr[opt].keys()] for opt in optimos]
                else:
                    energias = [[o_d for o_d in o_distr[opt].keys()] for opt in optimos]
                
                ocurrencias = [[o_d for o_d in o_distr[opt].values()] for opt in optimos]                
                
                print(optimos)
                # Dibujar los puntos
                datax = []
                datay = []
                dataoc = []
                for x, opt in enumerate(optimos):
                    for y, energia in enumerate(energias[x]):
                        datax.append(opt)
                        datay.append(energia)
                        dataoc.append(ocurrencias[x][y])
                        #plt.scatter(opt, energia, s=ocurrencias[x][y]*10, c=ocurrencias[x][y], cmap='viridis', alpha=0.6)
                plt.scatter(datax, datay, c=dataoc, cmap='viridis', alpha=0.6, label=g, marker=marker[m])#, vmin=0, vmax=100)
                m+=1
                print('==================================')
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
# Configurar los ejes
#plt.xticks(range(5))
#plt.yticks(range(len(energias)))
plt.xlabel('Optimos')
plt.ylabel('Energias')
plt.title(f'Optimum distribution (bqm_emb approach)\n(For {vars} vars Max3SAT)')

# Añadir una barra de colores
cbar = plt.colorbar()
cbar.set_label('Ocurrencia')
    
plt.tight_layout()
plt.legend()
plt.show()

# %% [markdown]
# ### x:optimos ; y:E_0-min(E_subo) ; legend:gadgets

# %%
marker=['o', 'v', 'P', '*']

files = ['bqm_emb_scale'] #['EmbeddingComposite', 'QUBO_scaled', 'bqm_emb_scale'] #JUST 1!!!!

gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [12] #[5, 12, 50] #[5,10,12,20,50]
num_instances = 1

plt.figure(figsize=(6, 4))
m=0
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in range(num_instances):            
            for file_names in files:
                #Instance for getting later the clauses
                file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
                instance = getattr(g_module, g)(file_path)
                if file_names == 'bqm_emb_scale':
                    instance.fillQ()
                    Q = instance.Q
                    couplings=[]
                    for k in Q.keys():
                        if k[0]!=k[1]:
                            couplings.append(k)
                    embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            
                sample_set_str=''
                with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_{file_names}.txt', 'r') as f:
                    lines = f.readlines()
                    if file_names == 'EmbeddingComposite':
                        k=3
                    if file_names == 'QUBO_scaled' or file_names == 'bqm_emb_scale':
                        k=4
                    for r in range(k,len(lines)-1):
                        if r==k:
                            sample_set_str+=lines[r][55:]
                        elif r==len(lines)-2:
                            line=lines[r].replace(" ", "")
                            sample_set_str+=line[:-2]
                        else:
                            sample_set_str+=lines[r]
                    sample_set=eval(sample_set_str)
                    variables = [v for v in eval(lines[k-1][2:]).keys()]
                    scale_factor = eval(lines[0].split()[-1])

                o_distr = {}
                for sample in sample_set:
                    if file_names == 'bqm_emb_scale':
                        majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
                        #This is done for also SPIN sample_sets
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==1:
                                        majority_vote[var][0]+=1
                                    if value!=1:
                                        majority_vote[var][1]+=1
                                    if majority_vote[var][0]>majority_vote[var][2]/2:
                                        assignment[var] = 1
                                        break
                                    if majority_vote[var][1]>majority_vote[var][2]/2:
                                        assignment[var] = 0
                                        break
                                    if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                                        assignment[var] = random.choice([0,1])
                                        break
                                    '''if value==-1:
                                        assignment[var] = 0
                                    else:
                                        assignment[var] = value
                                    break'''

                    else:
                        assignment={variables[j]:sample[0][j] for j in range(len(sample[0]))}

                    o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
                    if o not in o_distr.keys():
                        o_distr[o] = {sample[-3]: sample[-2]}
                    else:
                        if sample[-3] not in o_distr[o].keys():
                            o_distr[o][sample[-3]] = sample[-2]
                        else:
                            o_distr[o][sample[-3]] += sample[-2]
                o_distr = dict(sorted(o_distr.items()))
                
                print('----------------------------------')
                print('Gadget: ', g, ' Vars: ', vars, ' Scaling_factor: ', scale_factor, ' ',file_names) 
                o_distr_pond_mean = {}
                for o, distr in o_distr.items():
                    sols = sum(value for value in distr.values())
                    o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
                o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
                print('Distribution\n', o_distr)
                print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)

                dif_e_opt = []
                optimos = list(o_distr.keys())
                for opt in optimos[1:]:
                    #dif_e_opt.append(min(o_distr[optimos[0]].keys())-min(o_distr[opt].keys()))
                    dif_e_opt.append(o_distr_pond_mean[optimos[0]]-o_distr_pond_mean[opt])

                plt.plot(optimos[1:], dif_e_opt, label=g, marker=marker[m], ls='')
                m+=1
                print('==================================')
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
#plt.xticks(range(1,13))
plt.xlabel('Suboptimos')
plt.ylabel('E_min (smallest opt) - E_min (opt)')
#plt.title(f'Difference between the smallest energy of the best suboptim found and \nthe minimum energy of each other suboptim \n(For {vars} vars Max3SAT)')
plt.title(f'Difference between the ponderate mean energy in the first suboptim found and\nthe ponderate mean energy in the others suboptims\n(For {vars} vars Max3SAT)')
plt.legend()
plt.show()

# %%

# %% [markdown]
# ### Mean difference between ponderate means (or minimum) of energies for each suboptim and instance

# %%
marker=['o', 'v', 'P', '*']

files = ['bqm_emb_scale'] #['EmbeddingComposite', 'QUBO_scaled', 'bqm_emb_scale'] #JUST 1!!!!

gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [12] #[5, 12, 50] #[5,10,12,20,50]
num_instances = 20

figure, axis = plt.subplots(2, 1, figsize=(6,10), constrained_layout = True) 
m=0
for g, g_module in gadgets.items():
    for vars in n_vars:
        dif_e_opt_min_all=[]
        dif_e_opt_mean_all=[]
        for i in range(num_instances):            
            for file_names in files:
                #Instance for getting later the clauses
                file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
                instance = getattr(g_module, g)(file_path)
                if file_names == 'bqm_emb_scale':
                    instance.fillQ()
                    Q = instance.Q
                    couplings=[]
                    for k in Q.keys():
                        if k[0]!=k[1]:
                            couplings.append(k)
                    embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            
                sample_set_str=''
                with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/bqm_emb_scale/p{vars}_{i}_{file_names}.txt', 'r') as f:
                    lines = f.readlines()
                    if file_names == 'EmbeddingComposite':
                        k=3
                    if file_names == 'QUBO_scaled' or file_names == 'bqm_emb_scale':
                        k=4
                    for r in range(k,len(lines)-1):
                        if r==k:
                            sample_set_str+=lines[r][55:]
                        elif r==len(lines)-2:
                            line=lines[r].replace(" ", "")
                            sample_set_str+=line[:-2]
                        else:
                            sample_set_str+=lines[r]
                    sample_set=eval(sample_set_str)
                    variables = [v for v in eval(lines[k-1][2:]).keys()]
                    scale_factor = eval(lines[0].split()[-1])

                o_distr = {}
                for sample in sample_set:
                    if file_names == 'bqm_emb_scale':
                        majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
                        #This is done for also SPIN sample_sets
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==1:
                                        majority_vote[var][0]+=1
                                    if value!=1:
                                        majority_vote[var][1]+=1
                                    if majority_vote[var][0]>majority_vote[var][2]/2:
                                        assignment[var] = 1
                                        break
                                    if majority_vote[var][1]>majority_vote[var][2]/2:
                                        assignment[var] = 0
                                        break
                                    if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                                        assignment[var] = random.choice([0,1])
                                        break
                                    '''if value==-1:
                                        assignment[var] = 0
                                    else:
                                        assignment[var] = value
                                    break'''

                    else:
                        assignment={variables[j]:sample[0][j] for j in range(len(sample[0]))}

                    o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
                    if o not in o_distr.keys():
                        o_distr[o] = {sample[-3]: sample[-2]}
                    else:
                        if sample[-3] not in o_distr[o].keys():
                            o_distr[o][sample[-3]] = sample[-2]
                        else:
                            o_distr[o][sample[-3]] += sample[-2]
                o_distr = dict(sorted(o_distr.items()))
                
                o_distr_pond_mean = {}
                for o, distr in o_distr.items():
                    sols = sum(value for value in distr.values())
                    o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
                o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))

                dif_e_opt_min = {}
                dif_e_opt_mean = {}
                optimos = list(o_distr.keys())
                for opt in optimos[1:]:
                    dif_e_opt_min[opt] = min(o_distr[optimos[0]].keys())-min(o_distr[opt].keys())
                    dif_e_opt_mean[opt] = o_distr_pond_mean[optimos[0]]-o_distr_pond_mean[opt]

            dif_e_opt_min_all.append(dif_e_opt_min)
            dif_e_opt_mean_all.append(dif_e_opt_mean)

        df = pd.DataFrame(dif_e_opt_min_all)
        mean_dif_e_opt_min = dict(df.mean())
        df = pd.DataFrame(dif_e_opt_mean_all)
        mean_dif_e_opt_mean = dict(df.mean())
        
        axis[0].plot(mean_dif_e_opt_min.keys(), mean_dif_e_opt_min.values(), label=g, marker=marker[m], ls='')
        axis[1].plot(mean_dif_e_opt_mean.keys(), mean_dif_e_opt_mean.values(), label=g, marker=marker[m], ls='')
        m+=1
        print('==================================')
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
#plt.xticks(range(1,13))
axis[0].set_xlabel('Suboptimos')
axis[1].set_xlabel('Suboptimos')
axis[0].set_ylabel('E_min (smallest opt) - E_min (opt)')
axis[1].set_ylabel('E_pm0 - E_pmi')
axis[0].set_title(f'Mean difference between the smallest energy of the best suboptim found and \nthe minimum energy of each other suboptim \n(For {vars} vars Max3SAT)')
axis[1].set_title(f'Mean difference between the ponderate mean energy in the first suboptim found and\nthe ponderate mean energy in the others suboptims\n(For {vars} vars Max3SAT)')
plt.legend()
plt.show()

# %%

# %% [markdown]
# ### x: Energies ; y: Ocurrences ; legend: gadgets (distribution of energies per suboptim)

# %%
marker=['o', 'v', 'P']

files = ['bqm_emb_scale'] #['EmbeddingComposite', 'QUBO_scaled', 'bqm_emb_scale']
suboptims = [0, 1, 2]

gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [12] #[5, 12, 50] #GO 1 AT A TIME!!!
num_instances = 1
for subo in suboptims:
    m=0
    for g, g_module in gadgets.items():
        for vars in n_vars:
            for i in range(num_instances):                
                for file_names in files:
                    #Instance for getting later the clauses
                    file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
                    instance = getattr(g_module, g)(file_path)
                    if file_names == 'bqm_emb_scale':
                        instance.fillQ()
                        Q = instance.Q
                        couplings=[]
                        for k in Q.keys():
                            if k[0]!=k[1]:
                                couplings.append(k)
                        embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
                
                    sample_set_str=''
                    with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_{file_names}.txt', 'r') as f:
                        lines = f.readlines()
                        if file_names == 'EmbeddingComposite':
                            k=3
                        if file_names == 'QUBO_scaled' or file_names == 'bqm_emb_scale':
                            k=4
                        for r in range(k,len(lines)-1):
                            if r==k:
                                sample_set_str+=lines[r][55:]
                            elif r==len(lines)-2:
                                line=lines[r].replace(" ", "")
                                sample_set_str+=line[:-2]
                            else:
                                sample_set_str+=lines[r]
                        sample_set=eval(sample_set_str)
                        variables = [v for v in eval(lines[k-1][2:]).keys()]
                        scale_factor = eval(lines[0].split()[-1])
    
                    o_distr = {}
                    for sample in sample_set:
                        if file_names == 'bqm_emb_scale':
                            #This is done for also SPIN sample_sets
                            majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
                            assignment = {}
                            for j, value in enumerate(sample[0]):
                                for var, q in embedding.items():
                                    if variables[j] in q:
                                        if value==1:
                                            majority_vote[var][0]+=1
                                        if value!=1:
                                            majority_vote[var][1]+=1
                                        if majority_vote[var][0]>majority_vote[var][2]/2:
                                            assignment[var] = 1
                                            break
                                        if majority_vote[var][1]>majority_vote[var][2]/2:
                                            assignment[var] = 0
                                            break
                                        if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                                            assignment[var] = random.choice([0,1])
                                            break
    
                        else:
                            assignment={variables[j]:sample[0][j] for j in range(len(sample[0]))}
    
                        
                        o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
                        if o not in o_distr.keys():
                            o_distr[o] = {sample[-3]: sample[-2]}
                        else:
                            if sample[-3] not in o_distr[o].keys():
                                o_distr[o][sample[-3]] = sample[-2]
                            else:
                                o_distr[o][sample[-3]] += sample[-2]
                    o_distr = dict(sorted(o_distr.items()))
    
                    print('----------------------------------')
                    print('Gadget: ', g, ' Vars: ', vars, ' Scaling_factor: ', scale_factor, ' ',file_names) 
                    print('Energy distribution:\n', o_distr)
                    o_distr_pond_mean = {}
                    for o, distr in o_distr.items():
                        sols = sum(value for value in distr.values())
                        o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
                    o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
                    print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)
    
                    #Graficar
                    optimos = list(o_distr.keys())
                    if file_names != 'EmbeddingComposite':
                        energias = [[o_d*scale_factor for o_d in o_distr[opt].keys()] for opt in optimos]
                    else:
                        energias = [[o_d for o_d in o_distr[opt].keys()] for opt in optimos]
                    
                    ocurrencias = [[o_d for o_d in o_distr[opt].values()] for opt in optimos]                
                    
                    plt.plot(o_distr[subo].keys(), o_distr[subo].values(), label=g, marker=marker[m], ls='')
                    m+=1
                    print('==================================')
    plt.xlabel('E')
    plt.ylabel('# ocurrences')
    plt.title('Suboptim: '+str(subo))
    plt.legend()
    plt.show()


# %%

# %% [markdown]
# ### Optimums for the Clauses + Auxiliar_variables

# %%
def find_suboptimums_clauses_aux(assignment, clauses, V):
    count = 0
    for i, c in enumerate(clauses):
        if assignment[V+i]==0:
            negated_aux = 1
        if assignment[V+i]==1:
            negated_aux = 0
        
        if ((c[0] > 0 and assignment[np.abs(c[0])-1] == assignment[V+i]) or (c[0] < 0 and assignment[np.abs(c[0])-1] == negated_aux)) or ((c[1] > 0 and assignment[np.abs(c[1])-1] == assignment[V+i]) or (c[1] < 0 and assignment[np.abs(c[1])-1] == negated_aux)):
            count += 1
        if (assignment[V+i] == 1) or ((c[2] > 0 and assignment[np.abs(c[2])-1] == 1) or (c[2] < 0 and assignment[np.abs(c[2])-1] == 0)):
            count += 1
                
    return str(2*len(clauses)-count)


# %%
files = ['EmbeddingComposite', 'QUBO_scaled', 'bqm_emb_scale']

gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [5, 12, 50] #[5,10,12,20,50]
num_instances = 1
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in range(num_instances):
            for file_names in files:
                #Instance for getting later the clauses
                file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
                instance = getattr(g_module, g)(file_path)
                if file_names == 'bqm_emb_scale':
                    instance.fillQ()
                    Q = instance.Q
                    couplings=[]
                    for k in Q.keys():
                        if k[0]!=k[1]:
                            couplings.append(k)
                    embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            
                sample_set_str=''
                with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_{file_names}.txt', 'r') as f:
                    lines = f.readlines()
                    if file_names == 'EmbeddingComposite':
                        k=3
                    if file_names == 'QUBO_scaled' or file_names == 'bqm_emb_scale':
                        k=4
                    for r in range(k,len(lines)-1):
                        if r==k:
                            sample_set_str+=lines[r][55:]
                        elif r==len(lines)-2:
                            line=lines[r].replace(" ", "")
                            sample_set_str+=line[:-2]
                        else:
                            sample_set_str+=lines[r]
                    sample_set=eval(sample_set_str)
                    variables = [v for v in eval(lines[k-1][2:]).keys()]
                    scale_factor = eval(lines[0].split()[-1])

                o_distr = {}
                for sample in sample_set:
                    if file_names == 'bqm_emb_scale':
                        #This is done for also SPIN sample_sets
                        majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==1:
                                        majority_vote[var][0]+=1
                                    if value!=1:
                                        majority_vote[var][1]+=1
                                    if majority_vote[var][0]>majority_vote[var][2]/2:
                                        assignment[var] = 1
                                        break
                                    if majority_vote[var][1]>majority_vote[var][2]/2:
                                        assignment[var] = 0
                                        break
                                    if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                                        assignment[var] = random.choice([0,1])
                                        break

                    else:
                        assignment={variables[j]:sample[0][j] for j in range(len(sample[0]))}

                    
                    o = int(find_suboptimums_clauses_aux(assignment, instance.clauses, instance.V))
                    if o not in o_distr.keys():
                        o_distr[o] = {sample[-3]: sample[-2]}
                    else:
                        if sample[-3] not in o_distr[o].keys():
                            o_distr[o][sample[-3]] = sample[-2]
                        else:
                            o_distr[o][sample[-3]] += sample[-2]
                o_distr = dict(sorted(o_distr.items()))
                print('----------------------------------')
                print('Gadget: ', g, ' Vars: ', vars, ' Scaling_factor: ', scale_factor, ' ',file_names)
                print('Energy distribution for each optimum:\n', o_distr)
                o_distr_pond_mean = {}
                for o, distr in o_distr.items():
                    sols = sum(value for value in distr.values())
                    o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
                o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
                print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)
        print('==================================')
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

# %%
