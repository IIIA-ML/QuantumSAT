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
import utils
import pickle

import numpy as np
import matplotlib.pyplot as plt
import random
from pathlib import Path
import dimod
import math
import multiprocessing

#from gurobipy import GRB, Model, LinExpr, sys

from minorminer import find_embedding
from dwave.system import DWaveSampler, FixedEmbeddingComposite
import dwave.embedding

# %% [markdown]
# # Generate pickle

# %% [markdown]
# Generate first the problem that is commented # that has 3 variables and 4 clauses (does not need dwave) and put peculiarity="_Manual" (when definint the instance). Then generate any problem you want with utils.generate_3sat().
# This is because the study in the next section starts with ploting the matrix of the Ising model of that problem with 3 vars.

# %% [markdown]
# ### Problem

# %%
random.seed(901)
num_vars = 50

dir = "../exp/eBeyond/Bian_study/Problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)

p = utils.generate_3sat(num_vars, ratio=4.2)
#p = "c generated problem\np cnf 3 4\n1 2 3 0\n-1 2 3 0\n1 -2 3 0\n 1 -2 -3 0\n"
file_name = "p"+str(num_vars)+".cnf"
file_path = f"{dir}/{file_name}"
with open(file_path,"w") as f:
    f.write(p)

# %% [markdown]
# ### Pickle data

# %%
token = "Your token"
#Other gadgets can be implemented if wanted
gadgets = [
    "CJ1",
    "CJ2",
    "CJ1_bian",
    "CJ2_bian"
]

# %%
# Generate a new Instance_edited.py with the changes needed and import it in the Notebook
# (.py example: chain_strength fixed to 2 ; CJ1_bian with weights in repeated literals of -1)
from Instance import Instance
created_instance = Instance(file_path)

# Write the peculiarity of the bian that will be writen in .pkl, ex: peculiarity='chain_strength2' -> p50_chain_strength2.pkl
peculiarity = ""

# %%
secure_question = input("Do you really want to use D-Wave? (y/n)")
for gadget in gadgets:
    getattr(created_instance, gadget).compute_scaled_bqm_embedded(token)
    if secure_question in ["Yes", "Y", "y", "yes"]:
        getattr(created_instance, gadget).solve_dwave_from_scaled_bqm_embedded(token, num_reads=100, annealing_time=100)

# %% [markdown]
# ### Save pickle

# %%
dir = "../exp/eBeyond/Bian_study/Pickles"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)
file_name = "p"+str(num_vars)+ ("_"+peculiarity+".pkl" if peculiarity != "" else ".pkl")
with open(p_dir / file_name, 'wb') as f:
    pickle.dump(created_instance, f)

# %% [markdown]
# # =======================================================================

# %% [markdown]
# # Study

# %%
token = "Your token"
# Other gadgets can be added if they have been implemented in pkl
gadgets = [
    "CJ1",
    "CJ2",
    "CJ1_bian",
    "CJ2_bian"
]

# %%
file_name = "../exp/eBeyond/Bian_study/Pickles/p50.pkl"
with open(file_name, 'rb') as f:
    instance = pickle.load(f)

# %% [markdown]
# ## No dwave solution

# %% [markdown]
# ### ISING

# %%
gadget = "CJ1"

file_name = "../exp/eBeyond/Bian_study/Pickles/p3_Manual.pkl"
with open(file_name, 'rb') as f:
    problem = pickle.load(f)
print("Num_vars: "+str(problem.N)+" ; Clauses: "+str(len(problem.clauses)))
num_vars = len(getattr(problem,gadget).bqm.linear)
matrix = np.zeros((num_vars, num_vars))
for i, value in getattr(problem,gadget).bqm.linear.items():
    matrix[i, i] = value
for (i, j), value in getattr(problem,gadget).bqm.quadratic.items():
    matrix[j, i] = value
fig, ax = plt.subplots()
ax.axis('tight')
ax.axis('off')

#Define labels for indexes so auxiliars are 'bi'
if gadget=="CJ1_bian" or gadget=="CJ2_bian":
    label = list(range(len(problem.clauses)*3))
    for i in range(len(problem.clauses)):
        label.append("b"+str(i))
else:
    label = list(range(problem.N))
    for i in range(len(problem.clauses)):
        label.append("b"+str(i))
table = ax.table(cellText=matrix, cellLoc='center', loc='center', colLabels=label, rowLabels=label, colWidths=[0.1]*num_vars)
table.scale(1, 1.5)
plt.title('BQM (Ising)')
plt.show()

print(getattr(problem,gadget).clauses)

# %% [markdown]
# ### Number of variables for each format (SAT, QUBO, QPU)

# %%
num_variables = {"SAT": [], "QUBO": [], "EMBEDDED": []}
for gadget in gadgets:
    num_variables["SAT"].append(instance.N)
    num_variables["QUBO"].append(len(set(getattr(instance,gadget).bqm.variables)))
    num_variables["EMBEDDED"].append(len(set(getattr(instance,gadget).bqm_embedded["bqm_embedded"].variables)))
pd.DataFrame(num_variables, index=gadgets)

# %% [markdown]
# ### Chains from bian, embedding and final

# %%
#Define real embedding
real_embedding = {}
for gadget in gadgets:
    if gadget in ["CJ1_bian", "CJ2_bian"]:
        real_embedding[gadget] = {}
        for key, values in getattr(instance,gadget).refactored_bian_embedding.items():
            combined_values = []
            for value in values:
                combined_values.extend(getattr(instance,gadget).bqm_embedded['embedding'].get(value, []))
            real_embedding[gadget][key] = combined_values
    else:
        real_embedding[gadget] = getattr(instance,gadget).bqm_embedded["embedding"]

# %%
chain_len = {("BIAN", "Max."): [], ("BIAN", "Mean"): [], ("EMBEDDING", "Max."): [], ("EMBEDDING", "Mean"): [], ("FINAL", "Max."): [], ("FINAL", "Mean"): []}
for gadget in gadgets:
    chain_len[("BIAN", "Max.")].append(max(len(chain) for chain in getattr(instance,gadget).refactored_bian_embedding.values()) if gadget=="CJ1_bian" or gadget=="CJ2_bian" else '-')
    chain_len[("BIAN", "Mean")].append(round(np.mean(list(len(chain) for chain in getattr(instance,gadget).refactored_bian_embedding.values())),2) if gadget=="CJ1_bian" or gadget=="CJ2_bian" else '-')
    chain_len[("EMBEDDING", "Max.")].append(max(len(chain) for chain in getattr(instance,gadget).bqm_embedded["embedding"].values()))
    chain_len[("EMBEDDING", "Mean")].append(round(np.mean(list(len(chain) for chain in getattr(instance,gadget).bqm_embedded["embedding"].values())),2))
    chain_len[("FINAL", "Max.")].append(max(len(chain) for chain in real_embedding[gadget].values()))
    chain_len[("FINAL", "Mean")].append(round(np.mean(list(len(chain) for chain in real_embedding[gadget].values())),2))
pd.DataFrame(chain_len, index=gadgets)

# %%
fig, ax = plt.subplots(figsize=(20, 8))
num_gadgets = len(gadgets)
bar_width = 1 / (num_gadgets + 1)  # Width of each bar

# Loop through each gadget and plot the bars with an offset
for i, gadget in enumerate(gadgets):
    chain_length = {}
    emb = getattr(instance, gadget).bqm_embedded['embedding']
    for chain in emb.values():
        chain_l = len(chain)
        if chain_l in chain_length:
            chain_length[chain_l] += 1
        else:
            chain_length[chain_l] = 1
    # Calculate the positions of the bars with an offset
    positions = np.array(list(chain_length.keys())) + i * bar_width - (num_gadgets - 1) * bar_width / 2
    
    # Plot the bars with the calculated positions
    ax.bar(positions, list(chain_length.values()), width=bar_width, label=gadget)

# Set the x-axis limits
# Set the x-axis ticks to be integers only
#ax.set_xticks(np.arange(0, max_len_all + 1))
#ax.set_xlim(7.5, 16.5)
#ax.set_ylim(0,20)
ax.set_title('Chain_length distribution')
ax.set_xlabel('Chain_length')
ax.set_ylabel('#')

# Add a legend
ax.legend(fontsize=12)
plt.show()


# %% [markdown]
# ## DWave solution study

# %%
def count_unsatisfied_clauses(assignment, clauses):
    count = 0
    for c in clauses:
        for l in c:
            if (l > 0 and assignment[np.abs(l)-1] == 1) or (l < 0 and assignment[np.abs(l)-1] == -1):
                count += 1
                break
                
    return str(len(clauses)-count)


# %%
def unchain_dwave_solutions(sample_set, embedding):
    majority_voting_dict = {}
    for key, value in embedding.items():
        chain_length = len(value)
        chain_df = dimod.keep_variables(sample_set, value).to_pandas_dataframe().iloc[:,0:chain_length]
        majority_voting_dict[key] = chain_df.mode(axis=1)[0]            
    solutions_bqm_embedded_df = pd.DataFrame(majority_voting_dict)
    solutions_bqm_embedded_df['Real Energy'] = sample_set.to_pandas_dataframe()['energy']
    solutions_bqm_embedded_df['Occurrences'] = sample_set.to_pandas_dataframe()['num_occurrences']
    return solutions_bqm_embedded_df


# %%
def chain_break_count(sample_set, embedding):
    chain_break = {}
    key_break = {}
    for key, value in embedding.items():
        chain_length = len(value)
        chain_df = dimod.keep_variables(sample_set, value).to_pandas_dataframe().iloc[:,0:chain_length]
        #Compute how many chains of len X break and how many chains from the variables Y breaks
        for chain in chain_df.to_numpy():
            if len(set(chain))!=1:
                if chain_length in chain_break.keys():
                    chain_break[chain_length] += 1/100 #Divided by 100 because it is computed for all the solutions from D-Wave
                else:
                    chain_break[chain_length] = 1/100 #Divided by 100 because it is computed for all the solutions from D-Wave
                if key in key_break.keys():
                    key_break[key] += 1
                else:
                    key_break[key] = 1
    return chain_break, key_break


# %% [markdown]
# ### Mean optimum found with D-Wave

# %%
unsatisfied_clauses = {}
for gadget in gadgets:
    unsatisfied_clauses[gadget] = []
    solutions_bqm_embedded_df_aux = unchain_dwave_solutions(getattr(instance,gadget).response_dwave_from_scaled_bqm_embedded, real_embedding[gadget])
    solutions_bqm_embedded_df = solutions_bqm_embedded_df_aux.loc[solutions_bqm_embedded_df_aux.index.repeat(solutions_bqm_embedded_df_aux.Occurrences)].reset_index(drop=True).drop(columns=['Occurrences'])
    for _, row in solutions_bqm_embedded_df.iterrows():
        assigment = dict(row)
        o = count_unsatisfied_clauses(assigment, instance.clauses)
        unsatisfied_clauses[gadget].append(int(o))
    print(f"{gadget} o_mean: {np.mean(unsatisfied_clauses[gadget])}")

# %% [markdown]
# ### DWave embedding: Chain_len, chain_breaks and study of which variables are broken

# %% [markdown]
# This is the same plot as the histogram we saw before, but has also the bars of the chain_breaks, although each time a chain is broken it only adds 1/100 (see function chain_break_count), because we are counting all the chains broken in the 100 runs of the solution.

# %%
chain_break = {}
variable_break = {}
fig, ax = plt.subplots(figsize=(20, 8))
num_gadgets = len(gadgets)
bar_width = 1 / (num_gadgets + 1)  # Width of each bar

for i, gadget in enumerate(gadgets):
    #Data for a histogram of chain_break (x-axis is chain_len) and Data for occurence of chain break for a QUBO variable
    chain_break[gadget], variable_break[gadget] = chain_break_count(getattr(instance,gadget).response_dwave_from_scaled_bqm_embedded, getattr(instance,gadget).bqm_embedded['embedding'])
    #Data for a histogram of chain_len
    chain_length = {}
    emb = getattr(instance, gadget).bqm_embedded['embedding']
    for chain in emb.values():
        chain_l = len(chain)
        if chain_l in chain_length:
            chain_length[chain_l] += 1
        else:
            chain_length[chain_l] = 1
    # Calculate the positions of the bars and plot the bar
    positions_chain_len = np.array(list(chain_length.keys())) + i * bar_width - (num_gadgets - 1) * bar_width / 2   
    ax.bar(positions_chain_len, list(chain_length.values()), width=bar_width, label=gadget)

    if chain_break[gadget]!={}:
        positions_chain_break = np.array(list(chain_break[gadget].keys())) + i * bar_width - (num_gadgets - 1) * bar_width / 2
        ax.bar(positions_chain_break, list(chain_break[gadget].values()), width=bar_width)

# Set the x-axis limits
# Set the x-axis ticks to be integers only
#ax.set_xticks(np.arange(0, max_len_all + 1))
#ax.set_xlim(7.5, 16.5)
#ax.set_ylim(0,20)
ax.set_title('Chain_length distribution, also ploted superposed the chain_break occurence')
ax.set_xlabel('Chain_length')
ax.set_ylabel('#')

# Add a legend
ax.legend(fontsize=12)
plt.show()

# %% [markdown]
# ### Which variable is broken (i.e. has its chain of dwave embedding broken), (occurrences, chain_length_DWave_embed.)

# %%
print("Table of the variable that is broken (rows), with (ocurrence, chain_length of dwave_embedding).\nThe occurrence is the sum of each time is broken for the 100 runs of the solution")
variable_break_df = {}
for gadget in gadgets:
    variable_break_df[gadget] = {}
    for variable, ocur in variable_break[gadget].items():
        variable_break_df[gadget][variable] = (ocur, len(getattr(instance,gadget).bqm_embedded['embedding'][variable]))
display(pd.DataFrame(variable_break_df))

# %% [markdown]
# ### Which literal is broken (l1 OR l2 OR l3)

# %%
literal_break = {}
for gadget in gadgets:
    lit_break = [0,0,0]
    if gadget=="CJ1_bian" or gadget=="CJ2_bian":
        for variable, ocur in variable_break[gadget].items():
            if variable%3==0:
                lit_break[2] += ocur
            elif variable%2==0:
                lit_break[1] += ocur
            else:
                lit_break[0] += ocur            
        literal_break[gadget] = lit_break
print("Literal that has its chain broken ; (l_1 OR l_2 OR l_3)")
pd.DataFrame(literal_break).transpose()

# %% [markdown]
# ### Study on Bian embedding

# %%
chain_break_bian = {}
variable_break_bian = {}
for gadget in ["CJ1_bian", "CJ2_bian"]:
    solutions_bian_embedding_df_aux = unchain_dwave_solutions(getattr(instance,gadget).response_dwave_from_scaled_bqm_embedded, getattr(instance,gadget).bqm_embedded['embedding'])
    chain_break_bian[gadget], variable_break_bian[gadget] = chain_break_count(getattr(instance,gadget).response_dwave_from_scaled_bqm_embedded, getattr(instance,gadget).refactored_bian_embedding)

# %%
#WE NEED TO CHANGE IN INSTANCE THE 'from_serializable(response_aux)' TO DO IT HERE ONCE WE LOAD THE PICKLE!
import dwave.inspector
response = dimod.SampleSet.from_serializable(instance.CJ2.response_dwave_from_scaled_bqm_embedded)
dwave.inspector.show(response)

# %%
