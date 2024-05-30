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
import pickle
import numpy as np
import matplotlib.pyplot as plt
import dimod
import pandas as pd
from IPython.display import Image
from pathlib import Path


# %%
num_vars = 50

# %%
with open(f'../exp/eBeyond/p{num_vars}.cnf.pkl', 'rb') as f:
    instance = pickle.load(f)

# %%
Image(filename='../figs/From_scaled_BQM_embedded_to_solutions.png') 

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
# # Number of physical qubits

# %%
for gadget in gadgets:
    emb = getattr(instance, gadget).bqm_embedded['embedding']
    print(f'{gadget}: {len(set([value for sublist in emb.values() for value in sublist]))}')
    plt.bar(gadget, len([value for sublist in emb.values() for value in sublist]))
    plt.title("Number of physical qubits")

# %% [markdown]
# # Chain_length distribution

# %%
max_len = 0
num_gadgets = len(gadgets)
bar_width = 1 / (num_gadgets + 1)  # Width of each bar

# Create a figure and axis for plotting
fig, ax = plt.subplots(figsize=(20, 8))

# Loop through each gadget and plot the bars with an offset
for i, gadget in enumerate(gadgets):
    chain_length = {}
    emb = getattr(instance, gadget).bqm_embedded['embedding']
    for chain in emb.values():
        chain_len = len(chain)
        if chain_len in chain_length:
            chain_length[chain_len] += 1
        else:
            chain_length[chain_len] = 1
        if chain_len > max_len:
            max_len = chain_len
    
    # Calculate the positions of the bars with an offset
    positions = np.array(list(chain_length.keys())) + i * bar_width - (num_gadgets - 1) * bar_width / 2
    
    # Plot the bars with the calculated positions
    ax.bar(positions, list(chain_length.values()), width=bar_width, label=gadget)

# Set the x-axis limits
# Set the x-axis ticks to be integers only
ax.set_xticks(np.arange(0, max_len + 1))
ax.set_title('Chain_length distribution')
ax.set_xlabel('Chain_length')
ax.set_ylabel('#')

# Add a legend
ax.legend(fontsize=12)
plt.show()


# %% [markdown]
# # D-Wave results

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
def count_unsatisfied_clauses_bian(assignment, clauses):
    count = 0
    i = 0
    for l in clauses:
            if (l[0] > 0 and assignment[i] == 1) or (l[0] < 0 and assignment[i] == -1) or (l[1] > 0 and assignment[i+1] == 1) or (l[1] < 0 and assignment[i+1] == -1) or (l[2] > 0 and assignment[i+2] == 1) or (l[2] < 0 and assignment[i+2] == -1):
                count += 1
            i += 3
                
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
def groupby_subopt(sample_set, embedding, SAT_num_vars, gadget, instance):
    solutions_bqm_embedded_df_aux = unchain_dwave_solutions(sample_set, embedding)
    solutions_bqm_embedded_df = solutions_bqm_embedded_df_aux.loc[solutions_bqm_embedded_df_aux.index.repeat(solutions_bqm_embedded_df_aux.Occurrences)].reset_index(drop=True).drop(columns=['Occurrences'])
    if gadget != "CJ1_bian" and gadget != "CJ2_bian":
        assigment_df = solutions_bqm_embedded_df[np.arange(0,SAT_num_vars)]
    else:
        assigment_df = solutions_bqm_embedded_df
    num_unsatisfied_clauses = []
    for _, row in assigment_df.iterrows():
        assigment = dict(row)
        if gadget != "CJ1_bian" and gadget != "CJ2_bian":
            o = count_unsatisfied_clauses(assigment, instance.clauses)
        else:
            o = count_unsatisfied_clauses_bian(assigment, instance.clauses)
        num_unsatisfied_clauses.append(o)
    solutions_bqm_embedded_df['o'] = num_unsatisfied_clauses
    grouped_dfs = solutions_bqm_embedded_df.groupby('o')
    return grouped_dfs


# %%
for gadget in gadgets:
    grouped_dfs = groupby_subopt(sample_set=getattr(instance, gadget).response_dwave_from_scaled_bqm_embedded, embedding=getattr(instance, gadget).bqm_embedded['embedding'], SAT_num_vars=num_vars, gadget=gadget, instance=instance)
    for group_key, sub_opt_df in grouped_dfs:
        if group_key in ['0', '1']:
            plt.bar(np.unique(np.round(sub_opt_df['Real Energy'].values, 2), return_counts=True)[0],
                    np.unique(np.round(sub_opt_df['Real Energy'].values, 2), return_counts=True)[1], width=0.01,
                    alpha=1 - 0.2 * int(group_key), label=group_key)
            plt.legend()
    plt.title(f"{gadget}")
    plt.show()

# %% [markdown]
# ## GUROBI results

# %%
from gurobipy import Model, LinExpr, GRB, sys


# %%
def solve_gurobi(scaled_bqm_embedded, num_solutions=100):
    qubits = list(scaled_bqm_embedded.variables)
    q_names = list(map(str,scaled_bqm_embedded.variables))

    model = Model(name = 'linear program')
    q=model.addVars(qubits, name=q_names, vtype = GRB.INTEGER, lb=-1, ub=1)
   
    for qub in qubits:
        model.addConstr(q[qub]*q[qub] == 1, name='non_zero')

    obj_fn = LinExpr()
    for (u, v), bias in scaled_bqm_embedded.quadratic.items():
        obj_fn += bias * q[u] * q[v]
    for v, bias in scaled_bqm_embedded.linear.items():
        obj_fn += bias * q[v]
    model.setObjective(obj_fn, GRB.MINIMIZE)

    solutions = []
    model.setParam(GRB.Param.PoolSolutions, num_solutions)
    model.setParam(GRB.Param.PoolSearchMode, 2)
    #model.setParam(GRB.Param.PoolGap, 100)
    model.setParam(GRB.Param.Cutoff, GRB.INFINITY)
    model.optimize()
   
    nSolutions = model.SolCount
    # print(f"Number of solutions found: {nSolutions}")
   
    for e in range(nSolutions):
        model.setParam(GRB.Param.SolutionNumber, e)
   
        # Status checking
        status = model.Status
        if status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
            print("The model cannot be solved because it is infeasible or unbounded")
            sys.exit(1)
        if status != GRB.OPTIMAL:
            print(f"Optimization was stopped with status {status}")
            sys.exit(1)
   
        energy = model.PoolObjVal
   
        assignment = {}
        for v in model.getVars():
            if round(v.Xn)==1:
                assignment[int(v.varName)] = 1
            else:
                assignment[int(v.varName)] = -1
   
        solutions.append([list(assignment.values()), energy])

    return solutions


# %%
solutions = solve_gurobi(instance.Nuesslein2.scaled_bqm_embedded['scaled_bqm_embedded'], num_solutions=100)

# %%
print(len(solutions[0][0]))

# %%
sample_set = instance.Nuesslein2.response_dwave_from_scaled_bqm_embedded
print(len(sample_set.to_pandas_dataframe().axes[1]))
print(sample_set.to_pandas_dataframe().axes[1])

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## Compute Table 3 paper with Bian

# %%
# Aixo es fent la mitja de tooots els shots (no es la Taula 3)
gadgets = ["CJ2_bian"]
for gadget in gadgets:
    mean_satisfied_clauses_per_instance = []
    for i in range(1):    
        with open(f'../exp/eBeyond/Pickles/p5/p5-{i}.pkl', 'rb') as f:
            p5 = pickle.load(f)
        solutions_bqm_embedded_df_aux = unchain_dwave_solutions(sample_set=getattr(p5, gadget).response_dwave_from_scaled_bqm_embedded, embedding=getattr(p5, gadget).bqm_embedded['embedding'])
        solutions_bqm_embedded_df = solutions_bqm_embedded_df_aux.loc[solutions_bqm_embedded_df_aux.index.repeat(solutions_bqm_embedded_df_aux.Occurrences)].reset_index(drop=True).drop(columns=['Occurrences'])
        print(solutions_bqm_embedded_df)
    #     if gadget != "CJ1_bian" and gadget != "CJ2_bian":
    #         assigment_df = solutions_bqm_embedded_df[np.arange(0,50)]
    #     else:
    #         assigment_df = solutions_bqm_embedded_df
    #     num_unsatisfied_clauses = []
    #     for _, row in assigment_df.iterrows():
    #         assigment = dict(row)
    #         # print(assigment)
    #         if gadget != "CJ1_bian" and gadget != "CJ2_bian":
    #             o = count_unsatisfied_clauses(assigment, p50.clauses)
    #         else:
    #             o = count_unsatisfied_clauses_bian(assigment, p50.clauses)
    #         num_unsatisfied_clauses.append(210-int(o))
    #     
    #     # print(np.unique(num_unsatisfied_clauses, return_counts=True))
    #     # print(np.mean(num_unsatisfied_clauses))
    #     mean_satisfied_clauses_per_instance.append(np.mean(num_unsatisfied_clauses))
    # print(np.round(np.mean(mean_satisfied_clauses_per_instance), 2))

# %%
print("BIAN EMBEDDING:",p5.CJ2_bian.bian_embedding) #{variable literal: variable sat}
# Necessitem {variable sat: [variable_literal_1, variable_literal_2, ...]} per poder cridar el unchain_dwave
refactored_bian_embedding = {}
for key, value in p5.CJ2_bian.bian_embedding.items():
    if value not in refactored_bian_embedding.keys():
        refactored_bian_embedding[value] = [key]
    else:
        refactored_bian_embedding[value].append(key)

print("REFACTORED BIAN EMBEDDING:",refactored_bian_embedding)
print("DWAVE EMBEDDING:", p5.CJ2_bian.bqm_embedded["embedding"])
real_embedding = {}
for key, values in refactored_bian_embedding.items():
    # Initialize a list to store the combined values
    combined_values = []
    # Iterate over each value in the list from the bian dictionary
    for value in values:
        # Extend the combined_values list with the list from the emb dictionary
        combined_values.extend(p5.CJ2_bian.bqm_embedded['embedding'].get(value, []))
    # Assign the combined list to the key in the new dictionary
    real_embedding[key] = combined_values
print("FINAL EMBEDDING:",real_embedding)

unchain_dwave_solutions(p5.CJ2_bian.response_dwave_from_scaled_bqm_embedded, real_embedding)

# %%
majority_voting_dict = {}
for key, value in real_embedding.items():
    print(key, value)
    chain_length = len(value)
    print(chain_length)
    chain_df = dimod.keep_variables(p5.CJ2_bian.response_dwave_from_scaled_bqm_embedded, value).to_pandas_dataframe().iloc[:,0:chain_length]
    print(chain_df)
    print(np.unique(chain_df.iloc[1,:].values, return_counts=True))
    break
    majority_voting_dict[key] = chain_df.mode(axis=1)[0]
    
            
    solutions_bqm_embedded_df = pd.DataFrame(majority_voting_dict)

# %%
solutions_bqm_embedded_df_aux = unchain_dwave_solutions(p5.CJ2_bian.response_dwave_from_scaled_bqm_embedded, real_embedding)
solutions_bqm_embedded_df = solutions_bqm_embedded_df_aux.loc[solutions_bqm_embedded_df_aux.index.repeat(solutions_bqm_embedded_df_aux.Occurrences)].reset_index(drop=True).drop(columns=['Occurrences'])


# %%
def count_unsatisfied_clauses_bian(assignment, clauses):
    count = 0
    for c in clauses:
        for l in c:
            if (l > 0 and assignment[np.abs(l)] == 1) or (l < 0 and assignment[np.abs(l)] == -1):
                count += 1
                break
    return str(len(clauses)-count)


# %%
unsatisfied_clauses = []
for _, row in solutions_bqm_embedded_df.iterrows():
    assigment = dict(row)
    # print(assigment)
    # print(assigment[1])
    o = count_unsatisfied_clauses_bian(assigment, p50.clauses)

    unsatisfied_clauses.append(int(o))

np.mean(unsatisfied_clauses)

# %%
assigment = solutions_bqm_embedded_df.iloc[1,]
print(assigment)

count_unsatisfied_clauses_bian(assigment, p50.clauses)

# %%
#S'ha de mirar un probema més petit i seguir què fan bian pas per pas i comparar-ho amb CJ2

# %%

# %%

# %%

# %%
# La Taula 3 es fent la mitja del millor shot per cada instancia
gadgets = ["CJ2"]
for gadget in gadgets:
    mean_satisfied_clauses_per_instance = []
    for i in range(20):    
        with open(f'../exp/eBeyond/Pickles/p50/p50-{i}.pkl', 'rb') as f:
            p50 = pickle.load(f)
        solutions_bqm_embedded_df_aux = unchain_dwave_solutions(sample_set=getattr(p50, gadget).response_dwave_from_scaled_bqm_embedded, embedding=getattr(p50, gadget).bqm_embedded['embedding'])
        solutions_bqm_embedded_df = solutions_bqm_embedded_df_aux.loc[solutions_bqm_embedded_df_aux.index.repeat(solutions_bqm_embedded_df_aux.Occurrences)].reset_index(drop=True).drop(columns=['Occurrences'])
        # if gadget != "CJ1_bian" and gadget != "CJ2_bian":
        #     assigment_df = solutions_bqm_embedded_df[np.arange(0,50)]
        # else:
        #     assigment_df = solutions_bqm_embedded_df
        assigment_df = solutions_bqm_embedded_df
        num_unsatisfied_clauses = []
        row = solutions_bqm_embedded_df.iloc[0,]
        assigment = dict(row)
        # if gadget != "CJ1_bian" and gadget != "CJ2_bian":
        #     o = count_unsatisfied_clauses(assigment, p50.clauses)
        # else:
        #     o = count_unsatisfied_clauses_bian(assigment, p50.clauses)
        o = count_unsatisfied_clauses(assigment, p50.clauses)
        print(o)
        num_unsatisfied_clauses.append(210-int(o))
        
        # print(np.unique(num_unsatisfied_clauses, return_counts=True))
        # print(np.mean(num_unsatisfied_clauses))
        mean_satisfied_clauses_per_instance.append(np.mean(num_unsatisfied_clauses))
    print(np.round(np.mean(mean_satisfied_clauses_per_instance), 2))

# %%
solutions_bqm_embedded_df.iloc[0,]

# %%
#Exact solver (Maxsatz)
exact_num_satisfied_clauses_per_instance = []
for i in range(20):
    dir = '../exp/eBeyond/Problems/p5'
    p_dir = Path(dir)
    file_name = f'{dir}/p{5}-{i}.cnf'
    response = !./../src/maxsatz {file_name}
    o = int((str(response).split("(minimum number of unsatisfied clauses) = ")[1].split("'")[0]))
    exact_num_satisfied_clauses_per_instance.append(21-o)


# %%
np.mean(exact_num_satisfied_clauses_per_instance)

# %%
