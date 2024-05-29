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
import random


# %%
num_vars = 5

# %%
with open(f'../exp/eBeyond/p{num_vars}.cnf.pkl', 'rb') as f:
#with open(f'../exp/eBeyond/Pickles/p50/p{num_vars}-0.pkl', 'rb') as f:
    instance = pickle.load(f)

# %%
Image(filename='../figs/From_scaled_BQM_embedded_to_solutions.png') 

# %% [markdown]
# ## Number physical qubits

# %%
gadgets = [
    "Nuesslein1",
    "Nuesslein2",
    "CJ1",
    "CJ2",
    "CJ1_bian",
    "CJ2_bian"
]

# %%
bian_log_qubits = 4*int(num_vars*4.2)
gadget_log_qubits = num_vars + int(num_vars*4.2)
for gadget in gadgets:
    emb = getattr(instance, gadget).bqm_embedded['embedding']
    if gadget=="CJ1_bian" or gadget=="CJ2_bian":
        extra_qubits=len(set([value for sublist in emb.values() for value in sublist]))-bian_log_qubits
        print(f'{gadget}: {len(set([value for sublist in emb.values() for value in sublist]))} || Extra qubits: {extra_qubits}')
    else:
        extra_qubits=len(set([value for sublist in emb.values() for value in sublist]))-gadget_log_qubits
        print(f'{gadget}: {len(set([value for sublist in emb.values() for value in sublist]))} || Extra qubits: {extra_qubits}')
    plt.bar(gadget, len([value for sublist in emb.values() for value in sublist]), label='Total qubits', alpha=0.7)
    plt.bar(gadget, extra_qubits, label='Extra_qubits', alpha=0.7)

    plt.title("Number of physical qubits and extra used in the embedding")

# %% [markdown]
# ## Chain_length distribution

# %%
'''max_len_dict = {}
for gadget in gadgets:
    max_len = 0
    chain_length = {}
    emb = getattr(instance, gadget).bqm_embedded['embedding']
    for chain in emb.values():
        if len(chain) in chain_length:
            chain_length[len(chain)] += 1
        if not len(chain) in chain_length:
            chain_length[len(chain)] = 1
        if len(chain)>max_len:
            max_len = len(chain)
    max_len_dict[gadget] = max_len
    plt.bar(list(chain_length.keys()), list(chain_length.values()), label=gadget)
    
max_len_df = pd.DataFrame(data=list(max_len_dict.values()), index=list(max_len_dict.keys()), columns=['Max_chain_length'])
display(max_len_df)
plt.xlim(1.5, max_len+0.5)
plt.title('Chain_length distribution')
plt.xlabel('Chain_length')
plt.ylabel('#')
plt.legend()
plt.show()'''

max_len_dict = {}
max_len_all = 0
num_gadgets = len(gadgets)
bar_width = 1 / (num_gadgets + 1)  # Width of each bar

# Create a figure and axis for plotting
fig, ax = plt.subplots(figsize=(20, 8))

# Loop through each gadget and plot the bars with an offset
for i, gadget in enumerate(gadgets):
    max_len = 0
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
        if chain_len > max_len_all:
            max_len_all = chain_len
    # Dict of max_chain_len
    max_len_dict[gadget] = max_len
    # Calculate the positions of the bars with an offset
    positions = np.array(list(chain_length.keys())) + i * bar_width - (num_gadgets - 1) * bar_width / 2
    
    # Plot the bars with the calculated positions
    ax.bar(positions, list(chain_length.values()), width=bar_width, label=gadget)

# Creat DataFrame for max_chain_len
max_len_df = pd.DataFrame(data=list(max_len_dict.values()), index=list(max_len_dict.keys()), columns=['Max_chain_length'])
display(max_len_df)

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
# ## D-Wave results

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
        if ((l[0] > 0 and assignment[i] == 1) or (l[0] < 0 and assignment[i] == -1) or (l[1] > 0 and assignment[i+1] == 1) or (l[1] < 0 and assignment[i+1] == -1) or (l[2] > 0 and assignment[i+2] == 1) or (l[2] < 0 and assignment[i+2] == -1)):
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
def groupby_subopt(sample_set, embedding, SAT_num_vars, gadget):
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
gadgets = [
    "Nuesslein1",
    "Nuesslein2",
    "CJ1",
    "CJ2",
    "CJ1_bian",
    "CJ2_bian"
]

# %%
for gadget in gadgets:
    grouped_dfs = groupby_subopt(sample_set=getattr(instance, gadget).response_dwave_from_scaled_bqm_embedded, embedding=getattr(instance, gadget).bqm_embedded['embedding'], SAT_num_vars=num_vars, gadget=gadget)
    for group_key, sub_opt_df in grouped_dfs:
        if group_key in ['0', '1', '2', '3']:
            plt.bar(np.unique(np.round(sub_opt_df['Real Energy'].values, 2), return_counts=True)[0],
                    np.unique(np.round(sub_opt_df['Real Energy'].values, 2), return_counts=True)[1], width=0.01,
                    alpha=1 - 0.2 * int(group_key), label=group_key)
            plt.legend()
    plt.title(f"{gadget}")
    plt.show()


# %% [markdown]
# ## DWave compare with exact

# %%
def majority_voting(embedding, sample):
    assignment = {}
    majority_vote = {var: [0,0,len(q)] for var, q in embedding.items()}
    for qubit, value in sample.items():
        for var, q in embedding.items():
            if qubit in q:
                if value==1:
                    majority_vote[var][0]+=1
                if value!=1:
                    majority_vote[var][1]+=1
                if majority_vote[var][0]>majority_vote[var][2]/2:
                    assignment[var] = 1
                    break
                if majority_vote[var][1]>majority_vote[var][2]/2:
                    assignment[var] = -1
                    break
                if majority_vote[var][0]==majority_vote[var][2]/2 and majority_vote[var][0]+majority_vote[var][1]==majority_vote[var][2]:
                    assignment[var] = random.choice([-1,1])
                    break

    return assignment


# %%
gadgets = [
    #"Nuesslein1",
    "Nuesslein2",
    "CJ1",
    "CJ2",
    "CJ1_bian",
    "CJ2_bian"
]
n_vars=[50]
iterations=20

o={}
for vars in n_vars:
    o_g={}
    o_ver={}
    for gadget in gadgets:
        o_list={}
        for i in range(iterations):
            #Maxsatz exact solution
            response = !./../src/maxsatz {f'../exp/e3/problems/p{vars}-{i}.cnf'}
            for l in reversed(response):
                if l.split()[0] == 'o':
                    optimum = int(l.split()[1])
                    break
            
            #Dwave pickle optimum distribution
            file_path = f'../exp/eBeyond/Pickles/p{vars}/p{vars}-{i}.pkl'
            with open(file_path, 'rb') as f:
                instance = pickle.load(f)
            sample_set = getattr(instance, gadget).response_dwave_from_scaled_bqm_embedded
            embedding = getattr(instance, gadget).bqm_embedded['embedding']
            for sample in sample_set.samples():
                assignment = majority_voting(embedding, sample)

                #For bian, do another majority voting to trsnform from bian_qubo variables to sat variables
                if gadget=='CJ1_bian' or gadget=='CJ2_bian':
                    #Rewrite the embedding because the one in the instance is not in a good format
                    embedding_bian = {}
                    for sat_v in sorted(set(getattr(instance, gadget).bian_embedding.values())):
                        embedding_bian[sat_v-1] = []
                        for bian_v, sat_var_emb in getattr(instance, gadget).bian_embedding.items():
                            if sat_var_emb==sat_v:
                                embedding_bian[sat_v-1].append(bian_v)
                    assignment_def = majority_voting(embedding_bian, assignment)
                    o_found = int(count_unsatisfied_clauses(assignment_def, getattr(instance, gadget).clauses))
                else:
                    assignment_def = assignment
                    o_found = int(count_unsatisfied_clauses(assignment_def, getattr(instance, gadget).clauses))
                
                if o_found-optimum not in o_list.keys():
                    o_list[o_found-optimum]=1
                else:
                    o_list[o_found-optimum]+=1
        o_g[gadget]={key:value/iterations for key, value in o_list.items()} #value/(iterations*num_reads)*100%
    o[f'(V={vars}, C={int(4.2*vars)})']=o_g


# %%
def sort_nested_dicts(d):
    # Create a new dictionary to hold the sorted structure
    sorted_dict = {}
    for key, value in d.items():
        if isinstance(value, dict):
            # Recursively sort nested dictionaries
            sorted_dict[key] = sort_nested_dicts(value)
        else:
            sorted_dict[key] = value
    # Sort the dictionary itself if it is a nested one
    return dict(sorted(sorted_dict.items()))

# Sort each inner dictionary
sorted_data = {outer_key: {inner_key: sort_nested_dicts(inner_value)
                           for inner_key, inner_value in outer_value.items()}
               for outer_key, outer_value in o.items()}
for k,v in sorted_data.items():
    for kk, vv in v.items():
        print(kk, ':', vv)

# %%
fig, axs = plt.subplots(3, 2, figsize=(8, 12))

i = 0
for o_vars in o.values():
    j = -1
    k=0
    for o_g_name, o_g in o_vars.items():
        x_pos=[x+j*0.18 for x in list(o_g.keys())]
        axs[i // 2, i % 2].bar(x_pos, o_g.values(), label=gadgets[k], width=0.2)
        j += 1
        k += 1
    axs[i // 2, i % 2].set_title(list(o.keys())[i])
    axs[i // 2, i % 2].set_xlabel('Optimum difference')
    axs[i // 2, i % 2].set_ylabel('Percentage (%)')
    if i==0:
        axs[i // 2, i % 2].set_ylim(0,90)
    if i in [1,2,3]:
        axs[i // 2, i % 2].set_ylim(0,50)
    if i!=4:
        axs[i // 2, i % 2].set_xlim([-0.6, 2.6])
        axs[i // 2, i % 2].set_xticks(np.arange(3))
    if i==4:
        axs[i // 2, i % 2].set_xlim([-0.6, 5.6])
        axs[i // 2, i % 2].set_xticks(np.arange(6))
        axs[i // 2, i % 2].set_ylim(0,16)
    axs[i // 2, i % 2].legend()
    i += 1

plt.tight_layout()
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
