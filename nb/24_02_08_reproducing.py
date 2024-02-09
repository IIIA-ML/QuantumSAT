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

# %% [markdown]
# # Reproducing Nuesslein's experiments

# %%
from generation import generate_3sat
import random
import numpy as np
from pathlib import Path
import os
import Nuesslein
from dwave.system import DWaveSampler, EmbeddingComposite
import time
from Gadgets_Carlos_Jordi import tseitin, three_sat_max2xor, regular_like, tree_like, clique
import dimod

# %% [markdown]
# ### Generate Random 3-SAT

# %%
random.seed(178)
n_vars = [5, 10, 12]

for vars in n_vars:
    dir = "../exp/iccs24/nuesslein_problems/p"+str(vars)
    p_dir = Path(dir)
    p_dir.mkdir(parents=True, exist_ok=True)
    for iter in range(20):
        p = generate_3sat(num_vars=vars, num_clauses=int(vars*4.2))
        with open(p_dir / ("p"+str(vars)+"_"+str(iter)+".cnf"),"w") as f:
            f.write(p)


# %% [markdown]
# ### Parse cnf

# %%
def parsear_cnf_file(file_name):
    clauses=[]
    with open(file_name, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith('c') and not line.startswith('p') and line:
                clauses.append(np.array(list(map(int, line.split()[:-1]))))
            if line.startswith('p'):
                b = int(line.split()[2]) +1
    return clauses, b


# %% [markdown]
# ### Utils

# %%
def count_unsatisfied_clauses(assignment, clauses):
    satisfied_count = 0
    for clause in clauses:
        for literal in clause:
            if (literal > 0 and assignment[literal - 1] == 1) or (literal < 0 and assignment[-literal - 1] == 0):
                # At least one literal is satisfied
                satisfied_count += 1
                break
    return len(clauses) - satisfied_count


# %%

# %% [markdown]
# ### Solve with MAXSATZ

# %%
def maxsatz(file_path):
    solution = ! /home/rocco/Maxsatz/maxsatz {file_path}
    return solution


# %%
# Solving all problmes with maxsatz
for folder in ["p5", "p10", "p12"]:
    directory = "../exp/iccs24/nuesslein_problems/"+folder
    for file in os.listdir(directory):
        if not file.startswith('.'):
            solution = maxsatz("../exp/iccs24/nuesslein_problems/"+folder+"/"+file)
            dir = "../exp/iccs24/nuesslein_maxsatz/"+folder
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            with open(p_dir / (file.replace(".cnf","_maxsatz.txt")),"w") as f:
                f.write(str(solution))

# %%
# Nuesslein's Tabe3 Maxsatz
unsolved_clauses_folder = []
for folder in ["p5", "p10", "p12"]:
    directory = "../exp/iccs24/nuesslein_maxsatz/"+folder
    unsolved_clauses_file = []
    for file in os.listdir(directory):
        if not file.startswith('.'):
            with open(directory+"/"+file, "r") as file:
                unsolved_clauses = int(file.readlines()[0].split("Optimal Solution (minimum number of unsatisfied clauses) = ")[1].split("'")[0])
            unsolved_clauses_file.append(unsolved_clauses)
    unsolved_clauses_folder.append(np.mean(unsolved_clauses_file))

print(unsolved_clauses_folder)


# %%

# %% [markdown]
# ### Get mean satisfied clauses with Simulated Annealing

# %%
def get_mean_results_sa(gadget: str):
    results = []
    if gadget in ["nuesslein1", "nuesslein2", "chancellor", "choi"]:
        module = __import__("Nuesslein")
        for folder in ["p5", "p10", "p12"]:
            directory = "../exp/iccs24/nuesslein_problems/"+folder
            solved_clauses = []
            for file in os.listdir(directory):
                if not file.startswith('.'):
                    clauses, V = parsear_cnf_file("../exp/iccs24/nuesslein_problems/"+folder+"/"+file)
                    instance = getattr(module, gadget)(clauses, V-1)
                    solved_clauses.append(instance.solve()[1])
            results.append(np.mean(solved_clauses))
    
    else:
        for folder in ["p5", "p10", "p12"]:
            directory = "../exp/iccs24/nuesslein_problems/"+folder
            solved_clauses = []
            for file in os.listdir(directory):
                if not file.startswith('.'):
                    clauses, V = parsear_cnf_file("../exp/iccs24/nuesslein_problems/"+folder+"/"+file)
                    h, J, aux = globals()[gadget](clauses, V)
                    sampler = dimod.SimulatedAnnealingSampler()
                    sample_set = sampler.sample_ising(h, J, num_reads=10)
                    solution = sample_set.first.sample
                    for key, val in solution.items():
                        if val == -1:
                            solution[key] = 0
                    solved_clauses.append(len(clauses)-count_unsatisfied_clauses(list(solution.values()), clauses))
            results.append(np.mean(solved_clauses))
                    
    return results


# %%
gadgets = ["nuesslein1", "nuesslein2", "chancellor", "choi", \
           "tseitin", "three_sat_max2xor", "regular_like", "tree_like", "clique"]
table3_sa = {}

for gadget in gadgets:
    print(f"Computing {gadget}...") 
    mean_results = get_mean_results_sa(gadget)
    table3_sa[gadget] = mean_results

print(table3_sa)


# %%

# %%

# %%

# %% [markdown]
# ### Get mean satisfied clauses with DWave

# %%
def solve_Q_dwave(Q, num_reads=100, anneal_time=1):
    token = 'Your token'
    sampler = EmbeddingComposite(DWaveSampler(token=token))
    sample_set = sampler.sample_qubo(Q, num_reads=num_reads, annealing_time=anneal_time, \
                                     return_embedding=True, reduce_intersample_correlation=True)

    return sample_set


# %%
def get_mean_results_dwave(gadget: str):
    module = __import__("Nuesslein")
    results = []
    for folder in ["p5", "p10", "p12"]:
        print(f"Computing folder {folder}")
        directory = "../exp/iccs24/nuesslein_problems/"+folder
        solved_clauses = []
        for file in os.listdir(directory):
            if not file.startswith('.'):
                clauses, V = parsear_cnf_file("../exp/iccs24/nuesslein_problems/"+folder+"/"+file)
                instance = getattr(module, gadget)(clauses, V-1)
                instance.fillQ()
                sample_set = solve_Q_dwave(instance.Q)
                unsatisfied_clauses = int(count_unsatisfied_clauses(sample_set.first.sample, clauses))
                solved_clauses.append(len(clauses)-unsatisfied_clauses)

                dir = "../exp/iccs24/nuesslein_dwave/"+gadget+"/"+folder
                save_dir = Path(dir)
                save_dir.mkdir(parents=True, exist_ok=True)
                with open((dir+"/"+file).replace('.cnf', '_dwave.txt'), 'a') as archivo:
                    archivo.write("o " + str(unsatisfied_clauses) + "\n")
                    archivo.write("e " + str(sample_set.first.energy) + "\n")
                    archivo.write("v " + str(sample_set.first.sample) + "\n")
                    archivo.write(str(sample_set.samples))
        results.append(np.mean(solved_clauses))

    return results


# %%
gadgets = ["nuesslein1", "nuesslein2", "chancellor", "choi"]
# gadgets = ["nuesslein2"]
table3_dwave = {}

for gadget in gadgets:
    print(f"Computing {gadget}...") 
    mean_results = get_mean_results_dwave(gadget)
    print(f"{gadget} successfully computed!")
    table3_dwave[gadget] = mean_results

# %%
table3_dwave

# %%

# %%
#ToDo: implement tseitin with dwave. Are the results even good??

# %%

# %%

# %%

# %%

# %%
