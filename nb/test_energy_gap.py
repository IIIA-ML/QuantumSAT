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
# ### Testing if bqm/scaling = bqm_emb/scaling

# %%
file = '../exp/e3/problems/p12-0.cnf'
instance = CJ1.CJ1(file)
instance.fillQ()
bqm = dimod.BinaryQuadraticModel.from_qubo(instance.Q)
couplings = []
for key in bqm.quadratic.keys():
    if key[0] != key[1]:
        couplings.append(key)
embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)

bqm_emb = dwave.embedding.embed_bqm(source_bqm=bqm, embedding=embedding, target_adjacency=DWaveSampler().adjacency)
bqm_emb.scale(1/2)

bqm.scale(1/2)
bqm_sc = dwave.embedding.embed_bqm(source_bqm=bqm, embedding=embedding, target_adjacency=DWaveSampler().adjacency)

#print(bqm_emb)
#print(bqm_sc)
if bqm_emb==bqm_sc:
    print('YES!')


# %%

# %% [markdown]
# ### Getting all possible solutions from a bqm

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
import itertools

def generate_possible_solutions(variables):
    return itertools.product([1, -1], repeat=variables)

def evaluate_energy(bqm, solution):
    energy = 0
    for (u, v), bias in bqm.quadratic.items():
        energy += bias * solution[u] * solution[v]
    for v, bias in bqm.linear.items():
        energy += bias * solution[v]
    return energy

def solve_bqm_all_solutions(bqm):
    all_solutions = []
    
    # Generate all possible combinations of values for the variables
    possible_solutions = generate_possible_solutions(len(bqm.variables))
    
    for solution in possible_solutions:
        energy = evaluate_energy(bqm, solution)
        all_solutions.append((solution, energy))
    
    return all_solutions

with open('../exp/e3/CNF_test.cnf', 'w') as f:
    f.write('c generated problem\np cnf 3 7\n')
    f.write('1 2 3 0\n1 2 -3 0\n1 -2 3 0\n1 -2 -3 0\n-1 2 3 0\n-1 -2 3 0\n -1 2 -3 0')

gadgets = {
    #'Nuesslein1': Nuesslein1,
    #'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    #'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[5] #[5,10,12,20,50]
num_instances = [0] #20
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in num_instances: #range(num_instances):
            file_path = '../exp/e3/CNF_test.cnf'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q
            h, J, e = dimod.qubo_to_ising(Q)
            print(h)
            print(J)
            bqm = dimod.BinaryQuadraticModel.from_ising(h,J)
            print(bqm)
            variables = sorted(bqm.variables)
            #bqm = dimod.BinaryQuadraticModel.from_qubo(Q)

# Solve for all possible solutions
all_solutions = solve_bqm_all_solutions(bqm)

# Print all solutions
o_distr = {}
for solution, energy in all_solutions:
    print(f"Solution: {solution}, Energy: {energy}")
    
    assignment={variables[j]:solution[j] if solution[j]!=-1 else 0 for j in range(len(solution))}
    o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
    if o not in o_distr.keys():
        o_distr[o] = {energy: 1}
    else:
        if energy not in o_distr[o].keys():
            o_distr[o][energy] = 1
        else:
            o_distr[o][energy] += 1
o_distr = dict(sorted(o_distr.items()))
print('Energy distribution for each optimum:\n', o_distr)
o_distr_pond_mean = {}
for o, distr in o_distr.items():
    sols = sum(value for value in distr.values())
    o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)

print('-----------------------------------------------------\nStudying suboptimal distribution for all the clauses (including Auxiliar Variables!)')
o_distr_clausesaux = {}
for solution, energy in all_solutions:
    o = int(find_suboptimums_clauses_aux(assignment, instance.clauses, instance.V))
    if o not in o_distr_clausesaux.keys():
        o_distr_clausesaux[o] = {energy: 1}
    else:
        if energy not in o_distr_clausesaux[o].keys():
            o_distr_clausesaux[o][energy] = 1
        else:
            o_distr_clausesaux[o][energy] += 1
o_distr = dict(sorted(o_distr_clausesaux.items()))
print('Energy distribution for each optimum:\n', o_distr_clausesaux)
o_distr_clausesaux_pond_mean = {}
for o, distr in o_distr_clausesaux.items():
    sols = sum(value for value in distr.values())
    o_distr_clausesaux_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
o_distr_clausesaux_pond_mean = dict(sorted(o_distr_clausesaux_pond_mean.items()))
print('Energy ponderate mean for each optimum:\n', o_distr_clausesaux_pond_mean)


# %%

# %% [markdown]
# ### Solve QUBO Embedding composite with autoscaling

# %%
gadgets = {
    'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars= [5,10,12,20,50]
num_instances = 1
for g, g_module in gadgets.items():
    for vars in n_vars:
        auto_scale=[]
        for i in range(num_instances):
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q
            sampler = EmbeddingComposite(DWaveSampler(token='Your token'))
            response = sampler.sample_qubo(Q, num_reads=100, annealing_time=100, reduce_intersample_correlation=True)
            dir = f'../exp/e3/dwave_scaling/{g}/p{vars}'
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            with open(f'{dir}/p{vars}_{i}_EmbeddingComposite.txt', 'w') as file:
                file.write("o "+utils.count_unsatisfied_clauses(response.first.sample, instance.clauses)+"\n")
                file.write("e "+str(response.first.energy)+"\n")
                file.write("v "+str(response.first.sample)+"\n")
                file.write(str(response.samples))

# %% [markdown]
# #### Study energy distribution for each optimum

# %%
gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    #'CJ1': CJ1,
    #'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[12] #[5,10,12,20,50]
num_instances = [6] #20
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in num_instances: #range(num_instances):
            #Instance for getting later the clauses
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            
            sample_set_str=''
            with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_dwave_autoscale.txt', 'r') as f:
                lines = f.readlines()
                for r in range(4,len(lines)-1):
                    if r==4:
                        sample_set_str+=lines[r][55:]
                    elif r==len(lines)-2:
                        line=lines[r].replace(" ", "")
                        sample_set_str+=line[:-2]
                    else:
                        sample_set_str+=lines[r]
                sample_set=eval(sample_set_str)
                variables = [v for v in eval(lines[3][2:]).keys()]
                scale_factor = eval(lines[0].split()[-1])
                
            o_distr = {}
            for sample in sample_set:
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
            print('Energy distribution for each optimum:\n', o_distr)
            o_distr_pond_mean = {}
            for o, distr in o_distr.items():
                sols = sum(value for value in distr.values())
                o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
            o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
            print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)

# %%

# %% [markdown]
# ### Solve QUBO scaled manually (embedding recicle or get new)

# %%
gadgets = {
    'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[5, 12, 50] #[5,10,12,20,50]
num_instances = 1 #20
for g, g_module in gadgets.items():
    for vars in n_vars:
        auto_scale=[]
        for i in range(num_instances):
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q
            h, J, e = dimod.qubo_to_ising(Q)
            bqm = dimod.BinaryQuadraticModel.from_ising(h,J)
            couplings=[]
            for k in Q.keys():
                if k[0]!=k[1]:
                    couplings.append(k)
            embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            '''with open(f'../exp/e3/DWave/{g}/p{vars}/p{vars}_{i}_dwave.txt', 'r') as f:
                lines=f.readlines()
                line=lines[-1]
                if g=='Nuesslein1' or g=='Nuesslein2':
                    line = eval(line[line.find('{'):-12])
                    embedding = line['embedding_context']['embedding']
                else:
                    line = eval(line[line.find('{'):-10])
                    embedding = {}
                    for k, v in line['embedding_context']['embedding'].items():
                        if k<=vars:
                            embedding[k-1]=v
                        else:
                            embedding[k-2]=v'''
            bqm_emb = dwave.embedding.embed_bqm(source_bqm=bqm, embedding=embedding, target_adjacency=DWaveSampler().adjacency)
            J_max = max(bqm_emb.quadratic.values())
            J_min = min(bqm_emb.quadratic.values())
            h_max = max(bqm_emb.linear.values())
            h_min = min(bqm_emb.linear.values())
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
            aut_sc = max(max(max(bqm_emb.linear.values())/4,0),max(min(bqm_emb.linear.values())/(-4),0),max(max(bqm_emb.quadratic.values())/1,0),max(min(bqm_emb.quadratic.values())/(-2),0),coupling_limit)
            auto_scale.append(1/aut_sc)
            if aut_sc>1.0:
                Q = {k: v/aut_sc for k,v in Q.items()}
                sampler = FixedEmbeddingComposite(DWaveSampler(token='DEV-291d80af600d6eb433a8019c579070ba37436e9a'), embedding=embedding)
                response = sampler.sample_qubo(Q, num_reads=100, annealing_time=100, return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            else:
                sampler = FixedEmbeddingComposite(DWaveSampler(token='Your token'), embedding=embedding)
                response = sampler.sample_qubo(Q, num_reads=100, annealing_time=100, return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            dir = f'../exp/e3/dwave_scaling/{g}/p{vars}'
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            with open(f'{dir}/p{vars}_{i}_QUBO_scaled.txt', 'w') as file:
                file.write("Scaling factor: " + str(aut_sc) + "\n")
                file.write("o "+utils.count_unsatisfied_clauses(response.first.sample, instance.clauses)+"\n")
                file.write("e "+str(response.first.energy)+"\n")
                file.write("v "+str(response.first.sample)+"\n")
                file.write(str(response.samples))

# %% [markdown]
# #### Distribution of all energies for each suboptimum (and ponderate mean)

# %%
gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    #'CJ1': CJ1,
    #'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[12] #[5,10,12,20,50]
num_instances = [6] #20
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in num_instances: #range(num_instances):
            #Instance for getting later the clauses
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            
            sample_set_str=''
            with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}.txt', 'r') as f:
                lines = f.readlines()
                for r in range(4,len(lines)-1):
                    if r==4:
                        sample_set_str+=lines[r][55:]
                    elif r==len(lines)-2:
                        line=lines[r].replace(" ", "")
                        sample_set_str+=line[:-2]
                    else:
                        sample_set_str+=lines[r]
                sample_set=eval(sample_set_str)
                variables = [v for v in eval(lines[3][2:]).keys()]
                scale_factor = eval(lines[0].split()[-1])
                
            o_distr = {}
            for sample in sample_set:
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
            print('Energy distribution for each optimum:\n', o_distr)
            o_distr_pond_mean = {}
            for o, distr in o_distr.items():
                sols = sum(value for value in distr.values())
                o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
            o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
            print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)

# %%

# %% [markdown]
# ### Solve BQM_emb scaled manually (new embedding)

# %%
gadgets = {
    'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[5, 12, 50] #[5,10,12,20,50]
num_instances = 1
for g, g_module in gadgets.items():
    for vars in n_vars:
        auto_scale=[]
        for i in range(num_instances):
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q
            h, J, e = dimod.qubo_to_ising(Q)
            bqm = dimod.BinaryQuadraticModel.from_ising(h,J)
            couplings=[]
            for k in Q.keys():
                if k[0]!=k[1]:
                    couplings.append(k)
            embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            '''with open(f'../exp/e3/DWave/{g}/p{vars}/p{vars}_{i}_dwave.txt', 'r') as f:
                lines=f.readlines()
                line=lines[-1]
                if g=='Nuesslein1' or g=='Nuesslein2':
                    line = eval(line[line.find('{'):-12])
                    embedding = line['embedding_context']['embedding']
                else:
                    line = eval(line[line.find('{'):-10])
                    embedding = {}
                    for k, v in line['embedding_context']['embedding'].items():
                        if k<=vars:
                            embedding[k-1]=v
                        else:
                            embedding[k-2]=v'''
            bqm_emb = dwave.embedding.embed_bqm(source_bqm=bqm, embedding=embedding, target_adjacency=DWaveSampler().adjacency)
            J_max = max(bqm_emb.quadratic.values())
            J_min = min(bqm_emb.quadratic.values())
            h_max = max(bqm_emb.linear.values())
            h_min = min(bqm_emb.linear.values())
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
            aut_sc = max(max(max(bqm_emb.linear.values())/4,0),max(min(bqm_emb.linear.values())/(-4),0),max(max(bqm_emb.quadratic.values())/1,0),max(min(bqm_emb.quadratic.values())/(-2),0),coupling_limit)
            auto_scale.append(aut_sc)

            #Generate embedding such as embed = {qubit: qubit}
            embedding_bqm_emb = {}
            for k in bqm_emb.linear.keys():
                if k not in embedding_bqm_emb.keys():
                    embedding_bqm_emb[k] = [k]
            for k in bqm_emb.quadratic.keys():
                if k[0] not in embedding_bqm_emb.keys():
                    embedding_bqm_emb[k[0]] = [k[0]]
                if k[1] not in embedding_bqm_emb.keys():
                    embedding_bqm_emb[k[1]] = [k[1]]
            if aut_sc>1.0:
                bqm_emb.scale(1/aut_sc)
                sampler = FixedEmbeddingComposite(DWaveSampler(token='Your token'), embedding=embedding_bqm_emb)
                response = sampler.sample(bqm_emb, num_reads=100, annealing_time=100, return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            else:
                sampler = FixedEmbeddingComposite(DWaveSampler(token='Your token'), embedding=embedding_bqm_emb)
                response = sampler.sample(bqm_emb, num_reads=100, annealing_time=100, return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            
            assignment = {}
            for qubit, value in response.first.sample.items():
                for var, q in embedding.items():
                    if qubit in q:
                        if value==-1:
                            assignment[var] = value+1
                        else:
                            assignment[var] = value
                        #if var in assignment.keys() and value!= assignment[var]:
                            #print('ERROR')
                        break
            dir = f'../exp/e3/dwave_scaling/{g}/p{vars}'
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            with open(f'{dir}/p{vars}_{i}_bqm_emb_scale.txt', 'w') as file:
                file.write("Scaling factor: " + str(aut_sc) + "\n")
                file.write("o "+utils.count_unsatisfied_clauses(assignment, instance.clauses)+"\n")
                file.write("e "+str(response.first.energy)+"\n")
                file.write("v "+str(response.first.sample)+"\n")
                file.write(str(response.samples))

# %% [raw]
# ADALT els ERRORs signifiquen que hi ha cadenes trencades, per tant no se quin valor posar a l'hora de generar l'assignment i calcular l'optim

# %% [markdown]
# #### Distribution of all energies for each suboptimum (and ponderate mean) (E FOR SPIN BQM -> NOT QUBO!)

# %%
gadgets = {
    #'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[12] #[5,10,12,20,50]
num_instances = [6] #20
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in num_instances: #range(num_instances):
            #Instance for getting later the clauses
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q
            couplings=[]
            for k in Q.keys():
                if k[0]!=k[1]:
                    couplings.append(k)
            embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            
            sample_set_str=''
            with open(f'../exp/e3/dwave_scaling/{g}/p{vars}/p{vars}_{i}_bqm_emb.txt', 'r') as f:
                lines = f.readlines()
                for r in range(4,len(lines)-1):
                    if r==4:
                        sample_set_str+=lines[r][55:]
                    elif r==len(lines)-2:
                        line=lines[r].replace(" ", "")
                        sample_set_str+=line[:-2]
                    else:
                        sample_set_str+=lines[r]
                sample_set=eval(sample_set_str)
                
                start_index = lines[-1].find("Variables(")
                end_index = lines[-1].find(")", start_index)
                variables = eval(lines[-1][start_index + len("Variables("):end_index])
                
                scale_factor = eval(lines[0].split()[-1])
                
            o_distr = {}
            for sample in sample_set:
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
                
                o = int(utils.count_unsatisfied_clauses(assignment, instance.clauses))
                if o not in o_distr.keys():
                    o_distr[o] = {sample[-3]: sample[-2]}
                else:
                    if sample[-3] not in o_distr[o].keys():
                        o_distr[o][sample[-3]] = sample[-2]
                    else:
                        o_distr[o][sample[-3]] += sample[-2]
            o_distr = dict(sorted(o_distr.items()))
            print('Energy distribution for each optimum:\n', o_distr)
            o_distr_pond_mean = {}
            for o, distr in o_distr.items():
                sols = sum(value for value in distr.values())
                o_distr_pond_mean[o] = sum(key*value/sols for key, value in distr.items())
            o_distr_pond_mean = dict(sorted(o_distr_pond_mean.items()))
            print('Energy ponderate mean for each optimum:\n', o_distr_pond_mean)

# %%

# %%

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
optimos = list(data.keys())
energias = list(data[optimos[0]].keys())
ocurrencias = [[data[opt][energia] for energia in energias] for opt in optimos]

# Crear la gráfica
plt.figure(figsize=(10, 6))

# Dibujar los puntos
for i, opt in enumerate(optimos):
    for j, energia in enumerate(energias):
        plt.scatter(i, j, s=ocurrencias[i][j]*10, c=ocurrencias[i][j], cmap='viridis', alpha=0.6)

# Configurar los ejes
plt.xticks(range(len(optimos)), optimos)
plt.yticks(range(len(energias)), energias)
plt.xlabel('Optimos')
plt.ylabel('Energias')
plt.title('Gráfico de Oportunidades')

# Añadir una barra de colores
cbar = plt.colorbar()
cbar.set_label('Ocurrencia')

plt.tight_layout()
plt.show()

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ### Compute scaling factor (embedding recicled)

# %%
gadgets = {
    'Nuesslein1': Nuesslein1,
    'Nuesslein2': Nuesslein2,
    'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[5,10,12,20,50]
num_instances = 1
for g, g_module in gadgets.items():
    auto_scale_v = {}
    print('For gadget: ', g)
    for vars in n_vars:
        auto_scale = []
        count_i=0
        count_f=0
        print('For n_vars: ', vars)
        for i in range(num_instances):
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q
            h, J, e = dimod.qubo_to_ising(Q)
            J_max = max(J.values())
            J_min = min(J.values())
            h_max = max(h.values())
            h_min = min(h.values())
            if J_max>1.0 or J_min<(-2.0) or h_max>4.0 or h_min<(-4.0):
                count_i+=1
            bqm = dimod.BinaryQuadraticModel.from_ising(h,J)
            '''with open(f'../exp/e3/DWave/{g}/p{vars}/p{vars}_{i}_dwave.txt', 'r') as f:
                lines=f.readlines()
                line=lines[-1]
                if g=='Nuesslein1' or g=='Nuesslein2':
                    line = eval(line[line.find('{'):-12])
                    embedding = line['embedding_context']['embedding']
                else:
                    line = eval(line[line.find('{'):-10])
                    embedding = {}
                    for k, v in line['embedding_context']['embedding'].items():
                        if k<=vars:
                            embedding[k-1]=v
                        else:
                            embedding[k-2]=v
            #embedding = utils.get_embedding(Q)'''
            couplings=[]
            for k in Q.keys():
                if k[0]!=k[1]:
                    couplings.append(k)
            embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
            bqm_emb = dwave.embedding.embed_bqm(source_bqm=bqm, embedding=embedding, target_adjacency=DWaveSampler().adjacency)
            J_max = max(bqm_emb.quadratic.values())
            J_min = min(bqm_emb.quadratic.values())
            h_max = max(bqm_emb.linear.values())
            h_min = min(bqm_emb.linear.values())
            if J_max>1.0 or J_min<(-2.0) or h_max>4.0 or h_min<(-4.0):
                count_f+=1
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
            auto_scale.append(max(max(max(bqm_emb.linear.values())/4,0),max(min(bqm_emb.linear.values())/(-4),0),max(max(bqm_emb.quadratic.values())/1,0),max(min(bqm_emb.quadratic.values())/(-2),0),coupling_limit))
        
        print('Range problem pre-embed: ', count_i)
        print('Range problem post-embed: ', count_f)
        print('Autoscale needed: ', auto_scale)

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ### Solve and unembed bqm from given parameters

# %%
gadgets = {
    #'Nuesslein1': Nuesslein1,
    #'Nuesslein2': Nuesslein2,
    #'CJ1': CJ1,
    'CJ2': CJ2,
}
#Information for experiment 1
n_vars=[12] #n_vars=[5,10,12,20,50]
num_instances = 1 #20
for g, g_module in gadgets.items():
    for vars in n_vars:
        for i in range(num_instances):
            with open(f'../exp/e3/embedding/{g}/p{vars}/p{vars}_{i}_embedding.txt', 'r') as f:
                lines = f.readlines()
                bqm_emb = eval(lines[2])
                embedding = eval(lines[1])
                auto_scale = float(lines[0].split()[-1])
            bqm_emb.scale(1/(auto_scale))
            response = DWaveSampler(token='Your token').sample(bqm_emb, num_reads=100, annealing_time=100, reduce_intersample_correlation=True, auto_scale=False)
            file_path = f'../exp/e3/problems/p{vars}-{i}.cnf'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q
            h, J, e = dimod.qubo_to_ising(Q)
            bqm = dimod.BinaryQuadraticModel.from_ising(h,J)
            sample_set = dwave.embedding.unembed_sampleset(response, embedding, source_bqm=bqm)
            print(sample_set.samples)
            dir = f'../exp/e3/sols_scaled/{g}/p{vars}'
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            with open(f'{dir}/p{vars}_{i}_dwave_scaled.txt', 'w') as f:
                f.write("o "+str(utils.count_unsatisfied_clauses(sample_set.first.sample, instance.clauses))+"\n")
                f.write("e "+str(sample_set.first.energy)+"\n")
                f.write("v "+str(sample_set.first.sample)+"\n")
                f.write(str(sample_set.info['embedding_context'])+"\n")
                f.write(str(sample_set.samples)+"\n")

# %%

# %%

# %%

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
import dwave.embedding
from dwave.embedding.chain_strength import uniform_torque_compensation
from dwave.system import DWaveSampler
import dwave_networkx as dnx
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite

import itertools

from warnings import warn

import dimod
import minorminer
import functools

from dwave.embedding import (target_to_source, unembed_sampleset, embed_bqm,
                             chain_to_quadratic, EmbeddedStructure)
from dwave.system.warnings import WarningHandler, WarningAction


# %% [markdown]
# ### Edited Embedding Composite

# %%
class EmbeddingCompositeEdited(dimod.ComposedSampler):
    """Maps problems to a structured sampler.

    Automatically minor-embeds a problem into a structured sampler such as a
    D-Wave system. A new minor-embedding is calculated each time one of its
    sampling methods is called.

    Args:
        child_sampler (:class:`dimod.Sampler`):
            A dimod sampler, such as a :obj:`.DWaveSampler`, that accepts
            only binary quadratic models of a particular structure.

        find_embedding (function, optional):
            A function `find_embedding(S, T, **kwargs)` where `S` and `T`
            are edgelists. The function can accept additional keyword arguments.
            Defaults to :func:`minorminer.find_embedding`.

        embedding_parameters (dict, optional):
            If provided, parameters are passed to the embedding method as
            keyword arguments.

        scale_aware (bool, optional, default=False):
            Pass chain interactions to child samplers that accept an `ignored_interactions`
            parameter.

        child_structure_search (function, optional):
            A function `child_structure_search(sampler)` that accepts a sampler
            and returns the :attr:`dimod.Structured.structure`.
            Defaults to :func:`dimod.child_structure_dfs`.

    Examples:

       >>> from dwave.system import DWaveSampler, EmbeddingComposite
       ...
       >>> sampler = EmbeddingComposite(DWaveSampler())
       >>> h = {'a': -1., 'b': 2}
       >>> J = {('a', 'b'): 1.5}
       >>> sampleset = sampler.sample_ising(h, J, num_reads=100)
       >>> sampleset.first.energy
       -4.5


    """
    def __init__(self, child_sampler,
                 find_embedding=minorminer.find_embedding,
                 embedding_parameters=None,
                 scale_aware=False,
                 child_structure_search=dimod.child_structure_dfs):

        self.children = [child_sampler]

        # keep any embedding parameters around until later, because we might
        # want to overwrite them
        self.embedding_parameters = embedding_parameters or {}
        self.find_embedding = find_embedding

        # set the parameters
        self.parameters = parameters = child_sampler.parameters.copy()
        parameters.update(chain_strength=[],
                          chain_break_method=[],
                          chain_break_fraction=[],
                          embedding_parameters=[],
                          return_embedding=[],
                          warnings=[],
                          )

        # set the properties
        self.properties = dict(child_properties=child_sampler.properties.copy())

        # track the child's structure. We use a dfs in case intermediate
        # composites are not structured. We could expose multiple different
        # searches but since (as of 14 june 2019) all composites have single
        # children, just doing dfs seems safe for now.
        self.target_structure = child_structure_search(child_sampler)

        self.scale_aware = bool(scale_aware)

    parameters = None  # overwritten by init
    """dict[str, list]: Parameters in the form of a dict.

    For an instantiated composed sampler, keys are the keyword parameters
    accepted by the child sampler and parameters added by the composite.
    """

    children = None  # overwritten by init
    """list [child_sampler]: List containing the structured sampler."""

    properties = None  # overwritten by init
    """dict: Properties in the form of a dict.

    Contains the properties of the child sampler.
    """

    return_embedding_default = False
    """Defines the default behaviour for :meth:`.sample`'s `return_embedding`
    kwarg.
    """

    warnings_default = WarningAction.IGNORE
    """Defines the default behavior for :meth:`.sample`'s `warnings` kwarg.
    """

    def sample(self, bqm, chain_strength=None,
               chain_break_method=None,
               chain_break_fraction=True,
               embedding_parameters=None,
               return_embedding=None,
               warnings=None,
               **parameters):
        """Sample from the provided binary quadratic model.

        Args:
            bqm (:obj:`dimod.BinaryQuadraticModel`):
                Binary quadratic model to be sampled from.

            chain_strength (float/mapping/callable, optional):
                Sets the coupling strength between qubits representing variables 
                that form a :term:`chain`. Mappings should specify the required 
                chain strength for each variable. Callables should accept the BQM 
                and embedding and return a float or mapping. By default, 
                `chain_strength` is calculated with
                :func:`~dwave.embedding.chain_strength.uniform_torque_compensation`.

            chain_break_method (function/list, optional):
                Method or methods used to resolve chain breaks. If multiple
                methods are given, the results are concatenated and a new field
                called "chain_break_method" specifying the index of the method
                is appended to the sample set.
                See :func:`~dwave.embedding.unembed_sampleset` and
                :mod:`dwave.embedding.chain_breaks`.

            chain_break_fraction (bool, optional, default=True):
                Add a `chain_break_fraction` field to the unembedded response with
                the fraction of chains broken before unembedding.

            embedding_parameters (dict, optional):
                If provided, parameters are passed to the embedding method as
                keyword arguments. Overrides any `embedding_parameters` passed
                to the constructor.

            return_embedding (bool, optional):
                If True, the embedding, chain strength, chain break method and
                embedding parameters are added to :attr:`dimod.SampleSet.info`
                of the returned sample set. The default behaviour is defined
                by :attr:`return_embedding_default`, which itself defaults to
                False.

            warnings (:class:`~dwave.system.warnings.WarningAction`, optional):
                Defines what warning action to take, if any. See
                :mod:`~dwave.system.warnings`. The default behaviour is defined
                by :attr:`warnings_default`, which itself defaults to
                :class:`~dwave.system.warnings.IGNORE`

            **parameters:
                Parameters for the sampling method, specified by the child
                sampler.

        Returns:
            :obj:`dimod.SampleSet`

        Examples:
            See the example in :class:`EmbeddingComposite`.

        """
        if return_embedding is None:
            return_embedding = self.return_embedding_default

        # solve the problem on the child system
        child = self.child

        # apply the embedding to the given problem to map it to the child sampler
        __, target_edgelist, target_adjacency = self.target_structure

        # add self-loops to edgelist to handle singleton variables
        source_edgelist = list(bqm.quadratic) + [(v, v) for v in bqm.linear]

        # get the embedding
        if embedding_parameters is None:
            embedding_parameters = self.embedding_parameters
        else:
            # we want the parameters provided to the constructor, updated with
            # the ones provided to the sample method. To avoid the extra copy
            # we do an update, avoiding the keys that would overwrite the
            # sample-level embedding parameters
            embedding_parameters.update((key, val)
                                        for key, val in self.embedding_parameters
                                        if key not in embedding_parameters)

        embedding = self.find_embedding(source_edgelist, target_edgelist,
                                        **embedding_parameters)

        if bqm and not embedding:
            raise ValueError("no embedding found")

        if not hasattr(embedding, 'embed_bqm'):
            embedding = EmbeddedStructure(target_edgelist, embedding)

        bqm_embedded = embedding.embed_bqm(bqm, chain_strength=chain_strength,
                                           smear_vartype=dimod.SPIN)

        if warnings is None:
            warnings = self.warnings_default
        elif 'warnings' in child.parameters:
            parameters.update(warnings=warnings)

        warninghandler = WarningHandler(warnings)

        warninghandler.chain_strength(bqm, embedding.chain_strength, embedding)
        warninghandler.chain_length(embedding)

        if 'initial_state' in parameters:
            # if initial_state was provided in terms of the source BQM, we want
            # to modify it to now provide the initial state for the target BQM.
            # we do this by spreading the initial state values over the
            # chains
            state = parameters['initial_state']
            parameters['initial_state'] = {u: state[v]
                                           for v, chain in embedding.items()
                                           for u in chain}

        if self.scale_aware and 'ignored_interactions' in child.parameters:

            ignored = []
            for chain in embedding.values():
                # just use 0 as a null value because we don't actually need
                # the biases, just the interactions
                ignored.extend(chain_to_quadratic(chain, target_adjacency, 0))

            parameters['ignored_interactions'] = ignored

        response = child.sample(bqm_embedded, **parameters)

        def async_unembed(response):
            # unembed the sampleset aysnchronously.

            warninghandler.chain_break(response, embedding)

            sampleset = unembed_sampleset(response, embedding, source_bqm=bqm,
                                          chain_break_method=chain_break_method,
                                          chain_break_fraction=chain_break_fraction,
                                          return_embedding=return_embedding)

            if return_embedding:
                sampleset.info['embedding_context'].update(
                    embedding_parameters=embedding_parameters,
                    chain_strength=embedding.chain_strength)

            if chain_break_fraction and len(sampleset):
                warninghandler.issue("All samples have broken chains",
                                     func=lambda: (sampleset.record.chain_break_fraction.all(), None))

            if warninghandler.action is WarningAction.SAVE:
                # we're done with the warning handler so we can just pass the list
                # off, if later we want to pass in a handler or similar we should
                # do a copy
                sampleset.info.setdefault('warnings', []).extend(warninghandler.saved)

            return sampleset

        return response, dimod.SampleSet.from_future(response, async_unembed)


# %%
gadgets = {
    'CJ2': CJ2
}
#Information for experiment 3
n_vars=[20]
iterations=1

for vars in n_vars:
    for i in range(iterations):
        file_path = f"../exp/e3/problems/p{vars}-{i}.cnf"
        
        #Quantum Annealing (D'Wave)
        for g, g_module in gadgets.items():
            dir = f"../exp/e3/DWave_egap/{g}/p{vars}"
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            file_name = f'{dir}/p{vars}-{i}.txt'
            file_name_real = f'{dir}/p{vars}-{i}_real.txt'
            instance = getattr(g_module, g)(file_path)
            instance.fillQ()
            Q = instance.Q 
            bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
            sampler = EmbeddingCompositeEdited(DWaveSampler(token='Your token'))
            response_real, sample_set_scale = sampler.sample(bqm, return_embedding=True, num_reads=100, annealing_time=100, reduce_intersample_correlation=True)
            with open(file_name, 'w') as f:
                f.write(str(sample_set_scale.samples))
            with open(file_name_real, 'w') as f_r:
                f_r.write(str(response_real.samples))
