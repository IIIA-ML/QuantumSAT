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
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import time
sys.path.append('../src')

import utils
import CJ2

from dwave.system import DWaveSampler, EmbeddingComposite
from minorminer import find_embedding

# %%

# %% [markdown]
# ### Generate problem

# %%
random.seed(901)
num_vars = 50
p = utils.generate_3sat(num_vars, ratio=4.2)

# %%

# %% [markdown]
# #### Save problem in exp/embedding/problems/

# %%
dir = "../exp/embedding/problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)

with open(p_dir / ("p"+str(num_vars)+".cnf"),"w") as f:
    f.write(p)

# %%

# %%

# %% [markdown]
# ### Perform embedding with different seeds
# ##### When solving with DWave, do the solutions depend strongly on the embedding?

# %%
dwave_token = "DEV-291d80af600d6eb433a8019c579070ba37436e9a"

# %%
instance = CJ2.CJ2("../exp/embedding/problems/p"+str(num_vars)+".cnf")
instance.fillQ()

# %%
for i in range(100,110):
    print(f"Finding Embedding #{i}...")
    sampler = EmbeddingComposite(DWaveSampler(token=dwave_token), embedding_parameters={"random_seed":i})
    print(f"Embedding #{i} found.")

    try:
    
        response = sampler.sample_qubo(instance.Q, num_reads=100, annealing_time=100, \
                                   reduce_intersample_correlation=True, return_embedding=True)
    
        print(f"Solved with Embedding #{i}. Appending in txt...")    
        with open(p_dir / ("p"+str(num_vars)+"_solutions.txt"),"a") as f:
            f.write("o "+str(utils.count_unsatisfied_clauses(response.first.sample, instance.clauses))+"\n")
            f.write("e "+str(response.first.energy)+"\n")
            f.write("v "+str(response.first.sample)+"\n")
            f.write(str(response.info['embedding_context'])+"\n")
            f.write(str(response.samples)+"\n")

    except:
        print("en verda no lmao")
        continue

# %%

# %%

# %% [markdown]
# ### Perform Embedding from SubEmbedding

# %% [markdown]
# #### Most seen variables in SAT isntance

# %% [markdown]
# ##### First, we will only embed the clauses with the variable that appears the most

# %%
u, count = np.unique(np.absolute(instance.clauses), return_counts=True)
count_sort_ind = np.argsort(-count)
u[count_sort_ind]

# %%
sub_clauses = []
for c in instance.clauses:
    if u[count_sort_ind][0] in c or -u[count_sort_ind][0] in c:
        sub_clauses.append(np.array(c))
    #elif u[count_sort_ind][1] in c or -u[count_sort_ind][1] in c:
    #    sub_clauses.append(np.array(c))

# %%
############### OHO! The aux variables need to be the same between the Q and the sub_Q!
############### Since they are indexed in order, we will sort the final Q having sub_Q at the beginning:

# %%
SC_tuple = [tuple(array) for array in sub_clauses]
C_tuple = [tuple(array) for array in instance.clauses]

incomplete_instance_tuple = [array for array in C_tuple if array not in SC_tuple]
incomplete_instance = [np.array(array) for array in incomplete_instance_tuple]

# %%
instance_sorted = CJ2.CJ2(clauses=np.concatenate((sub_clauses, incomplete_instance)), V=num_vars)

# %%
instance_sorted.fillQ()

# %%
print(len(SC_tuple))
print(len(C_tuple))
print(len(incomplete_instance))

# %%

# %%
sub_instance = CJ2.CJ2(clauses=sub_clauses, V=num_vars)

# %%
sub_instance.fillQ()

# %%
sub_embedding = utils.get_embedding(sub_instance.Q, dwave_token, random_seed=0)
sub_embedding

# %%
print(max(len(v) for v in sub_embedding.values()))

# %%
for i in range(100,110):
    print(f"Finding Embedding #{i}...")
    sampler = EmbeddingComposite(DWaveSampler(token=dwave_token), \
                             embedding_parameters={"random_seed":i, "initial_chains":sub_embedding})
    print(f"Embedding #{i} found.")

    try:    
        response = sampler.sample_qubo(instance_sorted.Q, num_reads=100, annealing_time=100, \
                               reduce_intersample_correlation=True, return_embedding=True)
    
        print(f"Solved with Embedding #{i}. Appending in txt...")    
        with open(p_dir / ("p"+str(num_vars)+"_subembedding_sorted"+"_solutions.txt"),"a") as f:
            f.write("o "+str(utils.count_unsatisfied_clauses(response.first.sample, instance.clauses))+"\n")
            f.write("e "+str(response.first.energy)+"\n")
            f.write("v "+str(response.first.sample)+"\n")
            f.write(str(response.info['embedding_context'])+"\n")
            f.write(str(response.samples)+"\n")
    except:
        print("en verdad no lmao")
        continue

# %%
with open('../exp/embedding/problems/p50_subembedding_sorted_solutions.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('o '):
            print(line)

# %%

# %% [markdown]
# #### Most seen variables in Q

# %%
non_zero_couplings = {}
for key, value in instance.Q.items():
        if key[0] != key[1]:
            if value != 0:
                non_zero_couplings[key] = value

u2, count2 = np.unique(np.absolute([key for key in non_zero_couplings]), return_counts=True)
count_sort_ind2 = np.argsort(-count2)
u2[count_sort_ind2]

# %%
sub_couplings = {}
for key, value in non_zero_couplings.items():
    if key[0] == u2[count_sort_ind2][0] or key[1] == u2[count_sort_ind2][0]:
        sub_couplings[key] = value
    elif key[0] == u2[count_sort_ind2][1] or key[1] == u2[count_sort_ind2][1]:
        sub_couplings[key] = value

# %%
sub_embedding2 = utils.get_embedding(sub_couplings, dwave_token, random_seed=0)
sub_embedding2

# %%
print(max(len(v) for v in sub_embedding2.values()))

# %%
for i in range(100,110):
    print(f"Finding Embedding #{i}...")
    sampler = EmbeddingComposite(DWaveSampler(token=dwave_token), \
                             embedding_parameters={"random_seed":i, "initial_chains":sub_embedding2})
    print(f"Embedding #{i} found.")

    try:
        response = sampler.sample_qubo(instance.Q, num_reads=100, annealing_time=100, \
                               reduce_intersample_correlation=True, return_embedding=True)
    
        print(f"Solved with Embedding #{i}. Appending in txt...")    
        with open(p_dir / ("p"+str(num_vars)+"_subembedding2"+"_solutions.txt"),"a") as f:
            f.write("o "+str(utils.count_unsatisfied_clauses(response.first.sample, instance.clauses))+"\n")
            f.write("e "+str(response.first.energy)+"\n")
            f.write("v "+str(response.first.sample)+"\n")
            f.write(str(response.info['embedding_context'])+"\n")
            f.write(str(response.samples)+"\n")
            
    except:
        print("En verdad no lmao")
        continue

# %%
with open('../exp/embedding/problems/p50_subembedding2_solutions_v2.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('o '):
            print(line)

# %%
o={}
most_seen_v = [2,4,10]
for i in most_seen_v:
    o_emb=[]
    with open(f'../exp/embedding/problems/p50_subembedding2_solutions_v{i}.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('o '):
                o_emb.append(int(line.split()[1]))
    o_emb.append(np.mean(o_emb))
    o[f'Emb2_v{i}'] = o_emb
iterations = [i for i in range(len(o_emb)-1)]
iterations.append('Mean')
pd.DataFrame(o, iterations)

# %% [markdown]
# ### Most seen (v_i,v_j) in SAT clauses

# %%
clauses, V = utils.parse_cnf_file("../exp/embedding/problems/p"+str(num_vars)+".cnf")
k = [[0,1],[0,2],[1,2]]
relevant_keys={}
for c in clauses:
    for key in k:
        if tuple(sorted([np.abs(c[key[0]])-1, np.abs(c[key[1]])-1])) not in relevant_keys.keys():
            relevant_keys[tuple(sorted([np.abs(c[key[0]])-1, np.abs(c[key[1]])-1]))] = 1
        else:
            relevant_keys[tuple(sorted([np.abs(c[key[0]])-1, np.abs(c[key[1]])-1]))] += 1
relevant_keys = dict(sorted(relevant_keys.items(), key=lambda x:x[1], reverse=True))
relevant_keys = dict(list(relevant_keys.items())[:len(relevant_keys)//8])

relevant_clauses = []
secundary_clauses = []
for c in clauses:
    if (tuple(sorted([np.abs(c[0])-1,np.abs(c[1])-1])) or tuple(sorted([np.abs(c[0])-1,np.abs(c[2])-1])) or tuple(sorted([np.abs(c[1])-1,np.abs(c[2])-1]))) in relevant_keys:
        relevant_clauses.append(c)
    else:
        secundary_clauses.append(c)

relevant_keys_Q = list(relevant_keys.keys())

for key in relevant_keys:
    for b in range(num_vars,num_vars+len(relevant_clauses)):
        relevant_keys_Q.append(tuple([key[0],b]))
        relevant_keys_Q.append(tuple([key[1],b]))

# %%
print(len(relevant_clauses))

# %%
instance = CJ2.CJ2(clauses=relevant_clauses+secundary_clauses, V=num_vars)
instance.fillQ()

# %%
for i in range(100,110):
    print(f"Finding Embedding #{i}...")
    couplings = []
    for key, value in instance.Q.items():
        if key in relevant_keys_Q:
            #if value != 0:
                couplings.append(key)
    initial_chains = find_embedding(couplings, DWaveSampler().edgelist, random_seed=0)
    sampler = EmbeddingComposite(DWaveSampler(token=dwave_token), embedding_parameters={"random_seed":i, "initial_chains":initial_chains})
    print(f"Embedding #{i} found.")

    try:
        #Q = {k: v for k, v in instance.Q.items() if (k[0] != k[1] and v != 0) or k[0] == k[1]}
        response = sampler.sample_qubo(instance.Q, num_reads=100, annealing_time=100, reduce_intersample_correlation=True, return_embedding=True)
        
        print(f"Solved with Embedding #{i}. Appending in txt...")    
        with open(p_dir / ("p"+str(num_vars)+"_solutions_approach_emb.txt"),"a") as f:
            f.write("o "+str(utils.count_unsatisfied_clauses(response.first.sample, instance.clauses))+"\n")
            f.write("e "+str(response.first.energy)+"\n")
            f.write("v "+str(response.first.sample)+"\n")
            f.write(str(response.info['embedding_context'])+"\n")
            f.write(str(response.samples)+"\n")
    except Exception as e:
        continue

# %%
count=0
for v in initial_chains.values():
    if len(v) > 40:
        count+=1
print('Num of chains that chain_length>40: ', count)
print('Maximum chain_length: ', max(len(v) for v in initial_chains.values()))

# %%
with open('../exp/embedding/problems/p50_solutions_approach_emb.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('o '):
            print(line)

# %% [markdown]
# ### Optimum distribution for each approach

# %%
o_v = {}
file_path = f'../exp/embedding/problems/p50.cnf'
clauses, v = utils.parse_cnf_file(file_path)
for i in most_seen_v:
    o_list = {}
    num_sols=0
    sample_set_str=""
    with open(f'../exp/embedding/problems/p50_subembedding2_solutions_v{i}.txt', 'r') as f:
        lines=f.readlines()
        r=0
        while r<len(lines):
            sample_set_str=''
            if lines[r].startswith('<bound'):
                for k in range(r,len(lines)):
                    if k==r:
                        sample_set_str+=lines[r][55:]
                    if lines[k].split()[0].startswith('(['):
                        if lines[k+1].split()[0].startswith('(['):
                            sample_set_str+=lines[k]
                        if not lines[k+1].split()[0].startswith('(['):
                            line=lines[k].replace(" ", "")
                            sample_set_str+=line[:-2]
                    if k!=r and not lines[k].split()[0].startswith('(['):
                        break
                r=k
                num_sols+=1 #To know how many solutions we have generated
                sample_set=eval(sample_set_str)
                for sample in sample_set:
                    assignment={i:sample[0][i] for i in range(len(sample[0]))}
                    o_found = int(utils.count_unsatisfied_clauses(assignment, clauses))
                    if o_found not in o_list.keys():
                        o_list[o_found]=sample[-2]
                    else:
                        o_list[o_found]+=sample[-2]
            r+=1
    o_v[f'V most seen: {i}']={key:value/num_sols for key, value in o_list.items()} #value/(num_sols*num_reads)*100%

# %%
print(o_v)

# %%
fig, ax = plt.subplots(figsize=(5, 5))

j = -1
k=0
for o_v_name, o_optims in o_v.items():
    x_pos=[x+j*0.2 for x in list(o_optims.keys())]
    ax.bar(x_pos, o_optims.values(), label=o_v_name, width=0.2)
    j += 1
ax.set_title('Optimum distribution for approach v_i in QUBO\n num_vars=50')
ax.set_xlabel('Optimum difference')
ax.set_ylabel('Percentage (%)')
#ax.set_ylim(0,20)
#ax.set_xlim([1.4, 10.6])
ax.set_xticks(np.arange(17))
ax.legend()

plt.show()

# %%
o_a = {}
z=0
file_path = f'../exp/embedding/problems/p50.cnf'
clauses, v = utils.parse_cnf_file(file_path)

approach_file = ['_solutions', '_subembedding_sorted_solutions_v1', '_subembedding2_solutions_v2', '_solutions_approach_emb_div8']
approach = ['MinorMiner', 'v_i in SAT', 'v_i in QUBO', '(v_i,v_j) in SAT']

for approach_f in approach_file:
    o_list = {}
    num_sols=0
    sample_set_str=""
    with open(f'../exp/embedding/problems/p50{approach_f}.txt', 'r') as f:
        lines=f.readlines()
        r=0
        while r<len(lines):
            sample_set_str=''
            if lines[r].startswith('<bound'):
                for k in range(r,len(lines)):
                    if k==r:
                        sample_set_str+=lines[r][55:]
                    if lines[k].split()[0].startswith('(['):
                        if lines[k+1].split()[0].startswith('(['):
                            sample_set_str+=lines[k]
                        if not lines[k+1].split()[0].startswith('(['):
                            line=lines[k].replace(" ", "")
                            sample_set_str+=line[:-2]
                    if k!=r and not lines[k].split()[0].startswith('(['):
                        break
                r=k
                num_sols+=1 #To know how many solutions we have generated
                sample_set=eval(sample_set_str)
                for sample in sample_set:
                    assignment={i:sample[0][i] for i in range(len(sample[0]))}
                    o_found = int(utils.count_unsatisfied_clauses(assignment, clauses))
                    if o_found not in o_list.keys():
                        o_list[o_found]=sample[-2]
                    else:
                        o_list[o_found]+=sample[-2]
            r+=1
    o_a[approach[z]]={key:value/num_sols for key, value in o_list.items()} #value/(num_sols*num_reads)*100%
    z+=1

# %%
fig, ax = plt.subplots(figsize=(5, 5))

j = -1
k=0
for o_a_name, o_optims in o_a.items():
    x_pos=[x+j*0.2 for x in list(o_optims.keys())]
    ax.bar(x_pos, o_optims.values(), label=approach[k], width=0.2)
    j += 1
    k+=1
ax.set_title('Optimum distribution for approach v_i in QUBO\n num_vars=50')
ax.set_xlabel('Optimum difference')
ax.set_ylabel('Percentage (%)')
#ax.set_ylim(0,20)
#ax.set_xlim([1.4, 10.6])
ax.set_xticks(np.arange(17))
ax.legend()

plt.show()

# %% [markdown]
# ### Scaling bqm and test Enery gap

# %%
instance = CJ2.CJ2("../exp/embedding/problems/p"+str(num_vars)+".cnf")
instance.fillQ()
h, J, e = dimod.qubo_to_ising(instance.Q)
bqm = dimod.BinaryQuadraticModel.from_ising(h,J)

embedding = []
with open("../exp/embedding/problems/p"+str(num_vars)+"_solutions.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith("{'embedding':"):
            line=eval(line)
            embedding.append(line['embedding'])

for i in range(len(embedding)):
    bqm_emb = dwave.embedding.embed_bqm(source_bqm=bqm, embedding=embedding[i], target_adjacency=DWaveSampler().adjacency)
    
    J_per_qubit = {}
    for key, value in bqm_emb.quadratic.items():
        if key[0] in J_per_qubit.keys():
            J_per_qubit[key[0]]+=value
        else:
            J_per_qubit[key[0]]=value
        if key[1] in J_per_qubit.keys():
            J_per_qubit[key[1]]+=value
        else:
            J_per_qubit[key[1]]=value            
    coupling_limit = max(max(max(J_per_qubit.values())/15,0),max(min(J_per_qubit)/(-18),0))
    auto_scale = max(max(max(bqm_emb.linear.values())/4,0),max(min(bqm_emb.linear.values())/(-4),0),max(max(bqm_emb.quadratic.values())/1,0),max(min(bqm_emb.quadratic.values())/(-2),0),coupling_limit)
    
    bqm_emb.scale(1/(auto_scale))
    
    response = DWaveSampler(token='DEV-291d80af600d6eb433a8019c579070ba37436e9a').sample(bqm_emb, num_reads=100, annealing_time=100, reduce_intersample_correlation=True, auto_scale=False)
    with open(p_dir / ("p"+str(num_vars)+"_solutions_scaled.txt"),"w") as f:
        f.write("o "+str(utils.count_unsatisfied_clauses(response.first.sample, instance.clauses))+"\n")
        f.write("e "+str(response.first.energy)+"\n")
        f.write("v "+str(response.first.sample)+"\n")
        f.write(str(response.info['embedding_context'])+"\n")
        f.write(str(response.samples)+"\n")

# %%

# %%

# %%
