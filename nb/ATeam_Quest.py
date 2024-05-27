# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
from IPython.display import Image
from src.Instance import Instance
import utils
import pickle

import numpy as np
import matplotlib.pyplot as plt
import random
import time
from pathlib import Path
import dimod

import itertools
import multiprocessing

from gurobipy import GRB, Model, LinExpr, sys

from minorminer import find_embedding
from dwave.system import DWaveSampler, FixedEmbeddingComposite
import dwave.embedding

# %% [markdown]
# # Chapter 1
# # From SAT to scaled BQM embedded 

# %%
Image(filename='../figs/scaled_BQM_embedded_quest.png') 

# %% [markdown]
# ## Generate SAT Problem

# %%
random.seed(901)
num_vars = 4

dir = "../exp/eBeyond/"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)

p = utils.generate_3sat(num_vars, ratio=4.2)
file_name = "p"+str(num_vars)+".cnf"
with open(p_dir / file_name,"w") as f:
    f.write(p)

# %%
file_path = dir + file_name

# %% [markdown]
# ## Initialize instance

# %%
instance = Instance(file_path)

# %% [markdown]
# ## Select Gadget

# %%
gadgets = [
    "Nuesslein1",
    "Nuesslein2",
    "CJ1",
    "CJ2",
    "CJ1_bian",
    "CJ2_bian"
]

# %% [markdown]
# ## Compute BQM (N + aux)

# %%
for gadget in gadgets:
    getattr(instance, gadget).compute_bqm()
    print(f"BQM computed using {gadget}")

# %%
# instance.CJ1.bqm

# %%
###############################################################################################################################################

# %%
###############################################################################################################################################

# %% [markdown]
# ## Brute force study of all solution space with energy landscape

# %%
Image(filename='../figs/Solution_space_BQM_N_aux.png')


# %%
def energy_landscape_aux(possible_solutions, return_dict, bqm):
    all_solutions = []
    all_solutions_aux = []
    for solution in possible_solutions:
        energy = 0
        for (u, v), bias in bqm.quadratic.items():
            energy += bias * solution[u] * solution[v]
        for v, bias in bqm.linear.items():
            energy += bias * solution[v]
        all_solutions_aux.append((solution, energy))
    
    all_solutions.append(all_solutions_aux)
    return_dict[possible_solutions] = all_solutions


# %%
def compute_energy_landscape_in_parallel(gadget, parallel_batches, bqm):
    print(f"Using {gadget}...")
    possible_solutions = itertools.product([1, -1], repeat=len(bqm.variables))
    
    batches = {}
    batch_size = 2**len(bqm.variables) / parallel_batches
    assert batch_size.is_integer(), (f"The total number of solutions 2**{len(bqm.variables)} "
                                     f"is not divisible by the number of parallel batches ({parallel_batches}).")
    print(f"Number of solutions per batch: {int(batch_size)}")
    assert int(parallel_batches*batch_size) - int(2**len(bqm.variables)) == 0, "Batches not adding up"
    
    start_time = time.time()
    for i in range(parallel_batches):
        batches['batch'+str(i)] = itertools.islice(possible_solutions, int(i*batch_size), int((i+1)*batch_size))
        
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    for key, value in batches.items():
        p = multiprocessing.Process(target=energy_landscape_aux, args=(value, return_dict, bqm))
        jobs.append(p)
        p.start()
    
    for proc in jobs:
        proc.join()
    end_time = time.time()
    print("TIME:", round(end_time - start_time, 1), "seconds")
    
    all_possible_solutions = []
    for lst in return_dict.values():
        all_possible_solutions.extend(lst[0])
    # print(len(all_possible_solutions))
    assert 2**len(bqm.variables) - len(all_possible_solutions) == 0, "Solutions missing"
    return all_possible_solutions


# %%
def plot_energy_landscape(gadget, all_possible_solutions, bqm, instance, save_png=False):
    logical_variables = sorted(bqm.variables)
    o_distr = {}
    for solution, energy in all_possible_solutions:
        assignment={logical_variables[j]:solution[j] if solution[j]!=-1 else 0 for j in range(len(solution))}
        o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
        if o not in o_distr.keys():
            o_distr[o] = {energy: 1}
        else:
            if energy not in o_distr[o].keys():
                o_distr[o][energy] = 1
            else:
                o_distr[o][energy] += 1
                
    hist_distr = dict(sorted(o_distr.items()))    
    energy_values = [inner_dict_value for outer_dict_value in o_distr.values() for inner_dict_value in outer_dict_value.keys()]
    x_min = int(np.min(energy_values))
    x_max = int(np.max(energy_values))
    plt.figure(figsize=(20,4))
    plt.xticks(range(x_min,x_max+1,1), rotation=45)

    for i in range(len(list(hist_distr.keys())[:2])):
        plt.bar(list(hist_distr.values())[i].keys(), list(hist_distr.values())[i].values(), alpha=1-0.5*i, label=f'o{list(hist_distr.keys())[i]}')

    plt.title(f'{gadget}')
    plt.xlabel('Energy')
    plt.legend()
    if save_png is True:
        plt.savefig(f'../figs/{gadget}_p{instance.N}.png')
    plt.show()
 


# %%
brute_force_gadgets = [
    "Nuesslein1",
    "Nuesslein2",
    "CJ1",
    "CJ2"
]
parallel_batches = 8

# %%
print("Total number of solutions:")
print("--------------------------")
for g in brute_force_gadgets:
    print(f"2**{len(getattr(instance, g).bqm.variables)}={2**len(getattr(instance, g).bqm.variables)} for {g}")

# %%
all_possible_solutions = {}
for gadget in brute_force_gadgets:
    all_possible_solutions[gadget] = compute_energy_landscape_in_parallel(gadget, parallel_batches, getattr(instance, gadget).bqm)

# %%
for gadget in brute_force_gadgets:
    plot_energy_landscape(gadget, all_possible_solutions[gadget], getattr(instance, gadget).bqm, instance, save_png=False)

# %%
for gadget in brute_force_gadgets:
    getattr(instance, gadget).all_solution_space['bqm'] = all_possible_solutions[gadget]

# %%
##############################################################################################################################################

# %%
##############################################################################################################################################

# %% [markdown]
# ## Compute BQM embedded (N + aux + chains)

# %%
token = "DEV-2302e31be58c968e12b87cdb35d8b396d5c16ecd"

# %%
for gadget in gadgets:
    getattr(instance, gadget).compute_bqm_embedded(token)

# %%
# instance.CJ1.bqm_embedded

# %% [markdown]
# ## Compute scaling factor and scaled BQM embedded

# %%
for gadget in gadgets:
    getattr(instance, gadget).compute_scaled_bqm_embedded(token)
    print(gadget, getattr(instance, gadget).scale_factor)

# %% [markdown]
# # Chapter 2
# # From scaled BQM embedded to Solutions

# %%
Image(filename='../figs/From_scaled_BQM_embedded_to_solutions.png') 

# %% [markdown]
# ## Get D-Wave solutions (no auto-scale & no auto-embedding)

# %%
secure_question = input("Do you really want to use D-Wave? (y/n)")
if secure_question in ["Yes", "Y", "y", "yes"]:
    for gadget in gadgets:
        getattr(instance, gadget).solve_dwave_from_scaled_bqm_embedded(token, num_reads=100, annealing_time=100)

# %%
# instance.CJ1.response_dwave_from_scaled_bqm_embedded.samples

# %%

# %%

# %% [markdown]
# ## Save pickle file

# %%
with open(file_path+'.pkl', 'wb') as f:
    pickle.dump(instance, f)

# %%
