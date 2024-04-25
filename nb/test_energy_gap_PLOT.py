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
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==-1:
                                        assignment[var] = 0
                                    else:
                                        assignment[var] = value
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

# %%
marker=['o', 'v', 'P']

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
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==-1:
                                        assignment[var] = 0
                                    else:
                                        assignment[var] = value
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
                
                print(ocurrencias)
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
            plt.title('Gráfico de Oportunidades')
            
            # Añadir una barra de colores
            cbar = plt.colorbar()
            cbar.set_label('Ocurrencia')
                
            plt.tight_layout()
            plt.legend()
            plt.show()
        print('==================================')
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')


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
                        assignment = {}
                        for j, value in enumerate(sample[0]):
                            for var, q in embedding.items():
                                if variables[j] in q:
                                    if value==-1:
                                        assignment[var] = 0
                                    else:
                                        assignment[var] = value
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
