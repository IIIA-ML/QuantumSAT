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
import pandas as pd
import matplotlib.pyplot as plt

# %% [markdown]
# ### Generate Random 3-SAT

# %%
random.seed(178)
n_vars = [5, 10, 12, 20]

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
for folder in ["p5", "p10", "p12", "p20"]:
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
for folder in ["p5", "p10", "p12", "p20"]:
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
        for folder in ["p5", "p10", "p12", "p20"]:
        #for folder in ["p20"]:
            directory = "../exp/iccs24/nuesslein_problems/"+folder
            solved_clauses = []
            for file in os.listdir(directory):
                if not file.startswith('.'):
                    clauses, V = parsear_cnf_file("../exp/iccs24/nuesslein_problems/"+folder+"/"+file)
                    instance = getattr(module, gadget)(clauses, V-1)
                    solution = instance.solve()

                    dir = "../exp/iccs24/nuesslein_sa/"+gadget+"/"+folder
                    save_dir = Path(dir)
                    save_dir.mkdir(parents=True, exist_ok=True)
                    with open((dir+"/"+file).replace('.cnf', '_sa.txt'), 'w') as archivo:
                        archivo.write(str(solution))
                    
                    solved_clauses.append(solution[1])
            results.append(np.mean(solved_clauses))
    
    else:
        for folder in ["p5", "p10", "p12", "p20"]:
        #for folder in ["p20"]:
            directory = "../exp/iccs24/nuesslein_problems/"+folder
            solved_clauses = []
            for file in os.listdir(directory):
                if not file.startswith('.'):
                    clauses, V = parsear_cnf_file("../exp/iccs24/nuesslein_problems/"+folder+"/"+file)
                    h, J, aux = globals()[gadget](clauses, V)
                    sampler = dimod.SimulatedAnnealingSampler()
                    sample_set = sampler.sample_ising(h, J, num_reads=100)
                    solution = sample_set.first.sample
                    for key, val in solution.items():
                        if val == -1:
                            solution[key] = 0
                            
                    dir = "../exp/iccs24/nuesslein_sa/"+gadget+"/"+folder
                    save_dir = Path(dir)
                    save_dir.mkdir(parents=True, exist_ok=True)
                    with open((dir+"/"+file).replace('.cnf', '_sa.txt'), 'w') as archivo:
                        archivo.write(str(sample_set.samples))
                        
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
#hauríem de fer un loop que ens llegís els fitxers de SA i contés les clàusules per fer una taula tipo Taula3

# %%
# gadgets = ["nuesslein1", "nuesslein2", "chancellor", "choi", \
#             "tseitin", "three_sat_max2xor", "regular_like", "tree_like", "clique"]
# rows_df_sa = []
# for gadget in gadgets:
#     # print("Gadget: ", gadget)
#     satisfied_clauses_list = []
#     #for problem in ["p5", "p10", "p12"]:
#     for problem in ["p5", "p10", "p12", "p20"]:
#         # print("Problem size: ", problem)
#         num_satisfied_clausules = []
#         for i in range(0,20):
#             file = open(f'../exp/iccs24/nuesslein_sa/{gadget}/{problem}/{problem}_{i}_sa.txt', 'r')
#             for line in file:
#                 if line.startswith('o'):
#                     num_satisfied_clausules.append(int(int(problem.split('p')[1])*4.2) - float(line.split(' ')[1]))
#         # print("\t Satisfied clausules: ", np.mean(num_satisfied_clausules))
#         satisfied_clauses_list.append(np.mean(num_satisfied_clausules))
#     rows_df.append(np.concatenate(([gadget], satisfied_clauses_list)))

# %%

# %%

# %%

# %%

# %% [markdown]
# ### Get mean satisfied clauses with DWave

# %%
def solve_Q_dwave(Q, num_reads=100, anneal_time=100):
    token = 'DEV-eaed625d70ca6269d29ddf201caa6400e9eede52'
    sampler = EmbeddingComposite(DWaveSampler(token=token))
    sample_set = sampler.sample_qubo(Q, num_reads=num_reads, annealing_time=anneal_time, \
                                     return_embedding=True, reduce_intersample_correlation=True)

    return sample_set


# %%
def solve_h_J_dwave(h, J, num_reads=100, anneal_time=100):
    token = 'DEV-eaed625d70ca6269d29ddf201caa6400e9eede52'
    sampler = EmbeddingComposite(DWaveSampler(token=token))
    sample_set = sampler.sample_ising(h, J, num_reads=num_reads, annealing_time=anneal_time, \
                                     return_embedding=True, reduce_intersample_correlation=True)

    return sample_set


# %%
def get_mean_results_dwave(gadget: str):
    results = []
    if gadget in ["nuesslein1", "nuesslein2", "chancellor", "choi"]:
        module = __import__("Nuesslein")
        for folder in ["p5", "p10"]:
        #for folder in ["p20"]:
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
    
                    dir = "../exp/iccs24/nuesslein_dwave/100_100/"+gadget+"/"+folder
                    save_dir = Path(dir)
                    save_dir.mkdir(parents=True, exist_ok=True)
                    with open((dir+"/"+file).replace('.cnf', '_dwave.txt'), 'a') as archivo:
                        archivo.write("o " + str(unsatisfied_clauses) + "\n")
                        archivo.write("e " + str(sample_set.first.energy) + "\n")
                        archivo.write("v " + str(sample_set.first.sample) + "\n")
                        archivo.write(str(sample_set.samples))
            results.append(np.mean(solved_clauses))

    else:
        for folder in ["p5", "p10"]:
        #for folder in ["p20"]:
            print(f"Computing folder {folder}")
            directory = "../exp/iccs24/nuesslein_problems/"+folder
            solved_clauses = []
            for file in os.listdir(directory):
                if not file.startswith('.'):
                    clauses, V = parsear_cnf_file("../exp/iccs24/nuesslein_problems/"+folder+"/"+file)
                    h, J, aux = globals()[gadget](clauses, V)
                    sample_set = solve_h_J_dwave(h,J)
                    solution = sample_set.first.sample
                    for key, val in solution.items():
                        if val == -1:
                            solution[key] = 0
                    unsatisfied_clauses = count_unsatisfied_clauses(list(solution.values()), clauses)
                    solved_clauses.append(len(clauses)-unsatisfied_clauses)
   
                    dir = "../exp/iccs24/nuesslein_dwave/100_100/"+gadget+"/"+folder
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
# gadgets = ["nuesslein1", "nuesslein2", "chancellor", "choi", "tseitin"]
# gadgets = ["tseitin", "three_sat_max2xor", "nuesslein1", "nuesslein2"]
gadgets = ["tseitin"]
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
rows_df = []
for gadget in ["nuesslein1", "nuesslein2", "tseitin", "three_sat_max2xor"]:
    # print("Gadget: ", gadget)
    satisfied_clauses_list = []
    #for problem in ["p5", "p10", "p12"]:
    for problem in ["p5", "p10", "p12", "p20"]:
        # print("Problem size: ", problem)
        num_satisfied_clausules = []
        for i in range(0,20):
            file = open(f'../exp/iccs24/nuesslein_dwave/100_100/{gadget}/{problem}/{problem}_{i}_dwave.txt', 'r')
            for line in file:
                if line.startswith('o'):
                    num_satisfied_clausules.append(int(int(problem.split('p')[1])*4.2) - float(line.split(' ')[1]))
        # print("\t Satisfied clausules: ", np.mean(num_satisfied_clausules))
        satisfied_clauses_list.append(np.mean(num_satisfied_clausules))
    rows_df.append(np.concatenate(([gadget], satisfied_clauses_list)))

# %%
table3_df = pd.DataFrame(rows_df, columns = ["Gadget", "(V=5, C=21)", "(V=10, C=42)", "(V=12, C=50)", "(V=20, C=84)"])
table3_df.set_index("Gadget", inplace=True)

# %%
table3_df

# %%
# Incorporating MaxSatZ results to the df...:
solved_clauses_maxsatz = []
for folder in ["p5", "p10", "p12", "p20"]:
    directory = "../exp/iccs24/nuesslein_maxsatz/"+folder
    unsolved_clauses_file = []
    for file in os.listdir(directory):
        if not file.startswith('.'):
            with open(directory+"/"+file, "r") as file:
                unsolved_clauses = int(file.readlines()[0].split("Optimal Solution (minimum number of unsatisfied clauses) = ")[1].split("'")[0])
            unsolved_clauses_file.append(unsolved_clauses)
    solved_clauses_maxsatz.append(int(int(folder.split('p')[1])*4.2) - np.mean(unsolved_clauses_file))

table3_df.loc["maxsatz"] = solved_clauses_maxsatz

# %%
table3_df

# %%

# %% [markdown]
# ## Beyond Nuesslein's paper

# %%
print("p10 problem comparison with 10 reads vs 100 reads")
print("-------------------------------------------------")
for gadget in ["nuesslein1", "nuesslein2", "tseitin"]:
    for combination in ['10_100', '100_100']:
        num_satisfied_clausules = []
        for i in range(0,20):
            file = open(f'../exp/iccs24/nuesslein_dwave/{combination}/{gadget}/p10/p10_{i}_dwave.txt', 'r')
            for line in file:
                if line.startswith('o'):
                    num_satisfied_clausules.append(42 - float(line.split(' ')[1]))
        print(combination, gadget, np.mean(num_satisfied_clausules))
    print("-----------")

# %%

# %%
# To Do dimecres 15:
    # - Fer una taula3 com he fet a dwave pero amb simulated annealing
    # - Passarli el paper a firmar al Jesús¿?

# %%
for folder in ["p5", "p10", "p12", "p20"]:
    count_optimum_found_plot = []
    count_second_optimum_found_plot = []
    count_third_optimum_found_plot = []
    for gadget in ["nuesslein1", "nuesslein2", "tseitin", "three_sat_max2xor"]:
        print("Gadget:", gadget, "PROBLEM:", folder)
        count_optimum_found = 0
        count_second_optimum_found = 0
        count_third_optimum_found = 0
        for exp in range(0,20):
            file_ex = open(f'../exp/iccs24/nuesslein_dwave/100_100/{gadget}/{folder}/{folder}_{exp}_dwave.txt', 'r')
            content = file_ex.readlines() 
            sample_set = content[3:-1]
            solutions = []
            num_occurrences_list = []
            for shot in sample_set:
                shot = shot.replace(' ', '')
                if "array" in shot:
                    shot = shot.split('rec.array([')[1]
                # solutions_str.append(shot.split(']')[0][2:])
                solution = [int(numeric_str) for numeric_str in shot.split(']')[0][2:].split(',')]
                solution = (np.array(solution) + 1)/2
                solutions.append(solution)
                # num_occurrences_str.append(shot.split(',')[-3])
                num_occurrences = int(shot.split(',')[-3])
                num_occurrences_list.append(num_occurrences)
            
            clauses, V = parsear_cnf_file(f'../exp/iccs24/nuesslein_problems/{folder}/{folder}_{exp}.cnf')
            complete_opts = []
            for i in range(len(solutions)):
                opt = count_unsatisfied_clauses(solutions[i], clauses)
                complete_opts = np.concatenate((complete_opts, [opt]*num_occurrences_list[i]))
            
            maxsatz_file = open(f'../exp/iccs24/nuesslein_maxsatz/{folder}/{folder}_{exp}_maxsatz.txt', 'r')
            maxsatz_opt = int(maxsatz_file.readlines()[0].split("Optimal Solution (minimum number of unsatisfied clauses) = ")[1].split("'")[0])
            
            for opt in complete_opts:
                if opt == maxsatz_opt:
                    count_optimum_found += 1
                if opt == maxsatz_opt+1:
                    count_second_optimum_found += 1
                if opt == maxsatz_opt+2:
                    count_third_optimum_found += 1
        
            
        print(f"OPTIMUM FOUND {count_optimum_found} TIMES")

        count_optimum_found_plot.append(count_optimum_found)
        count_second_optimum_found_plot.append(count_second_optimum_found)
        count_third_optimum_found_plot.append(count_third_optimum_found)

    
    print(count_optimum_found_plot, count_second_optimum_found_plot, count_third_optimum_found_plot)
    plt.bar(["nuesslein1", "nuesslein2", "tseitin", "three_sat_max2xor"], count_optimum_found_plot, color="green",alpha=0.9)
    # plt.bar(["nuesslein1", "nuesslein2", "tseitin", "three_sat_max2xor"], count_second_optimum_found_plot, color="blue", alpha=0.7)
    # plt.bar(["nuesslein1", "nuesslein2", "tseitin", "three_sat_max2xor"], count_third_optimum_found_plot, color="red", alpha=0.5)
    plt.show()

