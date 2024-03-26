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
num_vars = 150
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
dwave_token = ""

# %%
instance = CJ2.CJ2("../exp/embedding/problems/p"+str(num_vars)+".cnf")
instance.fillQ()

# %%
for i in range(43,50):
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

# %%
##### First, we will only embed the clauses with the variable that appears the most

# %%
variables_appearence = np.unique(np.absolute(instance.clauses), return_counts=True)
most_appearing_variable = variables_appearence[0][np.argmax(variables_appearence[1])]

print(f"{most_appearing_variable} is the variable the appears the most")

# %%
variables_appearence

# %%
u, count = np.unique(np.absolute(instance.clauses), return_counts=True)
count_sort_ind = np.argsort(-count)
u[count_sort_ind]

# %%
sub_clauses = []
for c in instance.clauses:
    if most_appearing_variable in c or -most_appearing_variable in c:
        sub_clauses.append(np.array(c))
    elif 23 in c or -23 in c:
        sub_clauses.append(np.array(c))
    elif 129 in c or -129 in c:
        sub_clauses.append(np.array(c))
    elif 125 in c or -125 in c:
        sub_clauses.append(np.array(c))

sub_clauses


# %%
############### OHO! The aux variables need to be the same between the Q and the sub_Q!
############### Since they are indexed in order, we will sort the final Q having sub_Q at the beginning:

# %%
class sub_CJ2:

    def __init__(self, sub_clauses):
        self.clauses, self.V = sub_clauses, 150
        self.Q = {}

    def add(self, x, y, value):
        x = np.abs(x) - 1
        y = np.abs(y) - 1
        if x > y:
            x,y = y,x
        if (x,y) in self.Q.keys():
            self.Q[(x,y)] += value
        else:
            self.Q[(x,y)] = value

    def fillQ(self):
        for i, c in enumerate(self.clauses):
            s = [1 if l>0 else -1 for l in c]
            var = [abs(l) for l in c]
            self.add(var[0], var[0], s[0]-s[0]*s[2])
            self.add(var[1], var[1], -2*s[1])
            self.add(var[2], var[2], s[2]-s[0]*s[2])
            self.add(self.V+i+1, self.V+i+1, -1+s[0]-s[1]+s[2])
            self.add(var[0], var[2], 2*s[0]*s[2])
            self.add(var[0], self.V+i+1, -2*s[0])
            self.add(var[1], self.V+i+1, 2*s[1])
            self.add(var[2], self.V+i+1, -2*s[2])


# %%
SC_tuple = [tuple(array) for array in sub_clauses]
C_tuple = [tuple(array) for array in instance.clauses]

incomplete_instance_tuple = [array for array in C_tuple if array not in SC_tuple]
incomplete_instance = [np.array(array) for array in incomplete_instance_tuple]

# %%
instance_sorted = sub_CJ2(np.concatenate((sub_clauses, incomplete_instance)))
len(instance.clauses)

# %%
instance_sorted.fillQ()

# %%

# %%
sub_instance = sub_CJ2(sub_clauses)
sub_instance.Q

# %%
sub_instance.fillQ()

# %%
sub_embedding = utils.get_embedding(sub_instance.Q, dwave_token, random_seed=0)
sub_embedding

# %%

# %%
for i in range(0,10):
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
        if key[0] == 134 or key[1] == 134:
            sub_couplings[key] = value
        if key[0] == 148 or key[1] == 148:
            sub_couplings[key] = value
        if key[0] == 22 or key[1] == 22:
            sub_couplings[key] = value
        if key[0] == 124 or key[1] == 124:
            sub_couplings[key] = value

# %%
sub_embedding2 = utils.get_embedding(sub_couplings, dwave_token, random_seed=0)
sub_embedding2

# %%
with open(p_dir / ("p"+str(num_vars)+"_subembedding2"+"_solutions.txt"),"a") as f:
    f.write("Sub_embedding: "+str(sub_embedding2)+"\n")

# %%
for i in range(0,15):
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

# %%
