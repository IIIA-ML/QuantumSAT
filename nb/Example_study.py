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
import SAT_solver
import Splitter
import Gadget
import Joiner
import Solver
import random
from pathlib import Path
import utils
import pickle

# %% [markdown]
# ## Generate problem

# %%
random.seed(901)
num_vars = 100

dir = "../exp/eBeyond/Problems_SAT_solver"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)

p = utils.generate_3sat(num_vars, ratio=4.2)
#p = "c generated problem\np cnf 3 4\n1 3 2 0\n3 -1 2 0\n1 3 -2 0\n1 -3 -2 0\n"
#p = "c generated problem\np cnf 3 2\n1 3 -2 0\n-3 1 -2 0"
#p = "c generated problem\np cnf 3 1\n1 3 2 0"
file_name = "p"+str(num_vars)+".cnf"
file_path = f"{dir}/{file_name}"
with open(file_path,"w") as f:
    f.write(p)

# %% [markdown]
# ## SAT_Solver

# %%
splitter = Splitter.Single_problem()
gadget = [Gadget.CJ2()]
joiner = [Joiner.SAT_variables()]
solver = Solver.D_Wave()
token = 'Your token'

# %%
solution = SAT_solver.SAT_solver(file_path=file_path, splitter=splitter, gadget=gadget, joiner=joiner, solver=solver, token=token, qubit_level=True)
#solution = SAT_solver.SAT_solver(file_path=file_path, splitter=splitter, gadget=gadget, joiner=joiner)

# %%
solution.Solve()

# %% [markdown]
# ## Save pickle

# %%
feature = "CJ2_SATvars_ql"
dir = "../exp/eBeyond/Problems_SAT_solver/Pickles"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)
with open(f"../exp/eBeyond/Problems_SAT_solver/Pickles/p{solution.N}_{feature}.pkl", 'wb') as f:
    pickle.dump(solution, f)

# %%

# %% [markdown]
# ## Open pickle

# %%
with open(f'../exp/eBeyond/Problems_SAT_solver/Pickles/p50_CJ2_SATvars_ql.pkl', 'rb') as f:
    solution = pickle.load(f)

# %% [markdown]
# ## Ising matrix for CJ1 and CJ2

# %%
solution.clauses

# %%
print('CJ2')
solution.print_ISING()

# %%
print('CJ1')
solution.print_ISING()

# %% [markdown]
# ## Study

# %% [markdown]
# #### Optimum found

# %%
print("Solution found: ", dict(sorted(solution.o.items())))
print("Mean: ", sum(key * value for key, value in solution.o.items())/sum(solution.o.values()))
print("Exact solution: ", solution.exact_o)

# %%
print("Solution found: ", dict(sorted(solution.o.items())))
print("Mean: ", sum(key * value for key, value in solution.o.items())/sum(solution.o.values()))
print("Exact solution: ", solution.exact_o)

# %% [markdown]
# #### Non-zero couplings in Q

# %%
import dimod
J = dimod.qubo_to_ising(solution.Q)[1]
non_zero_couplings = 0
for key, val in J.items():
    if val!=0:
        non_zero_couplings += 1
print(non_zero_couplings)

# %% [markdown]
# #### Non-zero couplings in bqm_embedded

# %%
non_zero_couplings = 0
for key, val in solution.solver.bqm_embedded.quadratic.items():
    if val!=0:
        non_zero_couplings += 1
print(non_zero_couplings)

# %% [markdown]
# #### Number of couplings for each variable (total, only with other variables or with other auxiliar variables)

# %%
import dimod
J = dimod.qubo_to_ising(solution.Q)[1]
var_couplings_distr = {}
var_couplings_distr_var = {}
var_couplings_distr_aux = {}
for var in range(solution.N):
    var_couplings = []
    for key, val in J.items():
        if key[0]==var and val!=0:
            var_couplings.append(key[1])
        if key[1]==var and val!=0:
            var_couplings.append(key[0])
    if len(set(var_couplings)) not in var_couplings_distr.keys():
        var_couplings_distr[len(set(var_couplings))] = 1
    else:
        var_couplings_distr[len(set(var_couplings))] += 1
    var_couplings_var = [v for v in var_couplings if solution.Q_encoding.literal_to_variable[v]!=None]
    if len(set(var_couplings_var)) not in var_couplings_distr_var.keys():
        var_couplings_distr_var[len(set(var_couplings_var))] = 1
    else:
        var_couplings_distr_var[len(set(var_couplings_var))] += 1
    var_couplings_aux = [v for v in var_couplings if solution.Q_encoding.literal_to_variable[v]==None]
    if len(set(var_couplings_aux)) not in var_couplings_distr_aux.keys():
        var_couplings_distr_aux[len(set(var_couplings_aux))] = 1
    else:
        var_couplings_distr_aux[len(set(var_couplings_aux))] += 1
print('Distr. couplings for each varaible: ', dict(sorted(var_couplings_distr.items())), '\n   Between SAT_vars: ', dict(sorted(var_couplings_distr_var.items())), '\n   SAT_var with auxiliars: ', dict(sorted((var_couplings_distr_aux.items()))))

# %% [markdown]
# #### Chain_length (mean of all variables, SAT_variables and auxiliar)

# %%
chain_len_distr_var = {}
chain_len_distr_aux = {}
for var, chain in solution.embedding.items():
    if solution.Q_encoding.literal_to_variable[var]!=None:
        if len(chain) not in chain_len_distr_var.keys():
            chain_len_distr_var[len(chain)] = 1
        else:
            chain_len_distr_var[len(chain)] += 1
    else:
        if len(chain) not in chain_len_distr_aux.keys():
            chain_len_distr_aux[len(chain)] = 1
        else:
            chain_len_distr_aux[len(chain)] += 1
print('Chain length distribution (DWave_embedding): Mean: ', (sum(key * value for key, value in chain_len_distr_var.items())+sum(key * value for key, value in chain_len_distr_aux.items()))/(sum(chain_len_distr_var.values())+sum(chain_len_distr_aux.values())), '\n   Variables: ', dict(sorted(chain_len_distr_var.items())), '   Mean: ', sum(key * value for key, value in chain_len_distr_var.items())/sum(chain_len_distr_var.values()), '\n   Auxiliar: ', dict(sorted(chain_len_distr_aux.items())), '   Mean: ', sum(key * value for key, value in chain_len_distr_aux.items())/sum(chain_len_distr_aux.values()))

# %% [markdown]
# #### Max h, J for Q and bqm_embedded

# %%
max_h = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[0].values()])
max_J = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[1].values()])
max_h_emb = max([abs(val) for val in solution.solver.bqm_embedded.linear.values()])
max_J_emb = max([abs(val) for val in solution.solver.bqm_embedded.quadratic.values()])
print('Max |h|: ', max_h, '\nMax |J|: ', max_J, '\nMax |h| in bqm_embedded: ', max_h_emb, '\nMax |J| in bqm_embedded: ', max_J_emb)
if solution.response.info['embedding_context']['chain_strength']:
    print('Chain_strength: ', solution.response.info['embedding_context']['chain_strength'])

# %% [markdown]
# #### Compare all between different problems

# %%
import dimod
import numpy as np

#Possible comparisons: n={5,20,50,70}
n = 100
problems = ['_ql', '_twosub']#, '']
for problem in problems:
    print('-------------------------------------------------------------\n', problem[:-4])
    with open(f'../exp/eBeyond/Problems_SAT_solver/Pickles/p{n}_CJ2_SATvars{problem}.pkl', 'rb') as f:
        solution = pickle.load(f)

        max_h = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[0].values()])
        max_J = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[1].values()])
        if solution.qubit_level==True:
            max_h_emb = max([abs(val) for val in solution.solver.bqm_embedded.linear.values()])
            max_J_emb = max([abs(val) for val in solution.solver.bqm_embedded.quadratic.values()])
        print('Max |h|: ', max_h, '\nMax |J|: ', max_J)
        if solution.qubit_level==True:
            print('Max |h| in bqm_embedded: ', max_h_emb, ' If it is smaller it is because the h is divided with the qubits in the chain\nMax |J| in bqm_embedded: ', max_J_emb)
            print('Chain_strength: ', solution.response.info['embedding_context']['chain_strength'], ' ????????')

            print('Scale factor: ', solution.solver.scale_factor)
            print(len(set([q for qubits in solution.embedding.values() for q in qubits])))
            print(len(set([variables for variables in solution.embedding.keys()])))
        
        J = dimod.qubo_to_ising(solution.Q)[1]
        non_zero_couplings = 0
        for key, val in J.items():
            if val!=0:
                non_zero_couplings += 1
        print('Non 0 couplings in Q(ising): ', non_zero_couplings)

        if solution.qubit_level==True:
            non_zero_couplings = 0
            for key, val in solution.solver.bqm_embedded.quadratic.items():
                if val!=0:
                    non_zero_couplings += 1
            print('Non 0 couplings in bqm_embedded: ', non_zero_couplings)
        
        var_couplings_distr = {}
        var_couplings_distr_var = {}
        var_couplings_distr_aux = {}
        for var in range(solution.N):
            var_couplings = []
            for key, val in J.items():
                if key[0]==var and val!=0:
                    var_couplings.append(key[1])
                if key[1]==var and val!=0:
                    var_couplings.append(key[0])
            if len(set(var_couplings)) not in var_couplings_distr.keys():
                var_couplings_distr[len(set(var_couplings))] = 1
            else:
                var_couplings_distr[len(set(var_couplings))] += 1
            var_couplings_var = [v for v in var_couplings if solution.Q_encoding.literal_to_variable[v]!=None]
            if len(set(var_couplings_var)) not in var_couplings_distr_var.keys():
                var_couplings_distr_var[len(set(var_couplings_var))] = 1
            else:
                var_couplings_distr_var[len(set(var_couplings_var))] += 1
            var_couplings_aux = [v for v in var_couplings if solution.Q_encoding.literal_to_variable[v]==None]
            if len(set(var_couplings_aux)) not in var_couplings_distr_aux.keys():
                var_couplings_distr_aux[len(set(var_couplings_aux))] = 1
            else:
                var_couplings_distr_aux[len(set(var_couplings_aux))] += 1
        print('Distr. couplings for each varaible: ', dict(sorted(var_couplings_distr.items())), '   Mean: ', sum(key * value for key, value in var_couplings_distr.items())/sum(var_couplings_distr.values()), '\n   Between SAT_vars: ', dict(sorted(var_couplings_distr_var.items())), '   Mean: ', sum(key * value for key, value in var_couplings_distr_var.items())/sum(var_couplings_distr_var.values()), '\n   SAT_var with auxiliars: ', dict(sorted(var_couplings_distr_aux.items())), '   Mean: ', sum(key * value for key, value in var_couplings_distr_aux.items())/sum(var_couplings_distr_aux.values()))

        chain_len_distr_var = {}
        chain_len_distr_aux = {}
        for var, chain in solution.embedding.items():
            if solution.Q_encoding.literal_to_variable[var]!=None:
                if len(chain) not in chain_len_distr_var.keys():
                    chain_len_distr_var[len(chain)] = 1
                else:
                    chain_len_distr_var[len(chain)] += 1
            else:
                if len(chain) not in chain_len_distr_aux.keys():
                    chain_len_distr_aux[len(chain)] = 1
                else:
                    chain_len_distr_aux[len(chain)] += 1
        print('Chain length distribution (DWave_embedding): Mean: ', (sum(key * value for key, value in chain_len_distr_var.items())+sum(key * value for key, value in chain_len_distr_aux.items()))/(sum(chain_len_distr_var.values())+sum(chain_len_distr_aux.values())), '\n   Variables: ', dict(sorted(chain_len_distr_var.items())), '   Mean: ', sum(key * value for key, value in chain_len_distr_var.items())/sum(chain_len_distr_var.values()), '\n   Auxiliar: ', dict(sorted(chain_len_distr_aux.items())), '   Mean: ', sum(key * value for key, value in chain_len_distr_aux.items())/sum(chain_len_distr_aux.values()))

        print('Solution distribution: ', dict(sorted(solution.o.items())), '   Mean: ', sum(key * value for key, value in solution.o.items())/sum(solution.o.values()))

# %%
# Possible comparisons: n={5,20,50,55,70}
n = 5

for problem in ["","_ql"]:
    with open(f'../exp/eBeyond/Problems_SAT_solver/Pickles/p{n}_CJ2_SATvars{problem}.pkl', 'rb') as f:
        solution = pickle.load(f)
    print("Qubit_level = ", problem=="_ql")
    print("Solution found: ", dict(sorted(solution.o.items())))
    print("Mean: ", sum(key * value for key, value in solution.o.items())/sum(solution.o.values()))
    print("Exact solution: ", solution.exact_o)
    print("----------------------------------")

# %% [markdown]
# #### Plot with dwave.inspector

# %%
import dwave.inspector
from dwave.system import DWaveSampler
dwave.inspector.show(solution.solver.bqm_embedded, solution.response, DWaveSampler())
#dwave.inspector.show(solution.Q, solution.response, DWaveSampler())

# %% [markdown]
# #### Comparison between gadgets

# %%
import dimod

#Possible comparisons: n={5,20,50}
n = 20

for gadget in ["CJ1","CJ2"]:
    for problem in ['', '_ql']:
        with open(f'../exp/eBeyond/Problems_SAT_solver/Pickles/p{n}_{gadget}_SATvars{problem}.pkl', 'rb') as f:
            solution = pickle.load(f)
        print(gadget, ' ; Qubit_level=', problem=='_ql')
        print("Solution found: ", dict(sorted(solution.o.items())))
        print("Mean: ", sum(key * value for key, value in solution.o.items())/sum(solution.o.values()))
        print("Exact solution: ", solution.exact_o)
        if solution.qubit_level==True:
            max_h = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[0].values()])
            max_J = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[1].values()])
            max_h_emb = max([abs(val) for val in solution.solver.bqm_embedded.linear.values()])
            max_J_emb = max([abs(val) for val in solution.solver.bqm_embedded.quadratic.values()])
            print('Max |h|: ', max_h, '\nMax |J|: ', max_J)
            print('Max |h| in bqm_embedded: ', max_h_emb, '\nMax |J| in bqm_embedded: ', max_J_emb)

            print('Scale factor: ', solution.solver.scale_factor)
        print("----------------------------------")

# %%
import dimod
import numpy as np

#Possible comparisons: n={5,20,50}
n = 50
problems = ['_ql']
gadgets = ['CJ1', 'CJ2']
for problem in problems:
    for gadget in gadgets:
        print('-------------------------------------------------------------\n', gadget)
        with open(f'../exp/eBeyond/Problems_SAT_solver/Pickles/p{n}_{gadget}_SATvars{problem}.pkl', 'rb') as f:
            solution = pickle.load(f)
    
            max_h = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[0].values()])
            max_J = max([abs(val) for val in dimod.qubo_to_ising(solution.Q)[1].values()])
            if solution.qubit_level==True:
                max_h_emb = max([abs(val) for val in solution.solver.bqm_embedded.linear.values()])
                max_J_emb = max([abs(val) for val in solution.solver.bqm_embedded.quadratic.values()])
            print('Max |h|: ', max_h, '\nMax |J|: ', max_J)
            if solution.qubit_level==True:
                print('Max |h| in bqm_embedded: ', max_h_emb, ' If it is smaller it is because the h is divided with the qubits in the chain\nMax |J| in bqm_embedded: ', max_J_emb)
                print('Chain_strength: ', solution.response.info['embedding_context']['chain_strength'], ' ????????')
    
                print('Scale factor: ', solution.solver.scale_factor)
            
            J = dimod.qubo_to_ising(solution.Q)[1]
            non_zero_couplings = 0
            for key, val in J.items():
                if val!=0:
                    non_zero_couplings += 1
            print('Non 0 couplings in Q(ising): ', non_zero_couplings)
    
            if solution.qubit_level==True:
                non_zero_couplings = 0
                for key, val in solution.solver.bqm_embedded.quadratic.items():
                    if val!=0:
                        non_zero_couplings += 1
                print('Non 0 couplings in bqm_embedded: ', non_zero_couplings)
            
            var_couplings_distr = {}
            var_couplings_distr_var = {}
            var_couplings_distr_aux = {}
            for var in range(solution.N):
                var_couplings = []
                for key, val in J.items():
                    if key[0]==var and val!=0:
                        var_couplings.append(key[1])
                    if key[1]==var and val!=0:
                        var_couplings.append(key[0])
                if len(set(var_couplings)) not in var_couplings_distr.keys():
                    var_couplings_distr[len(set(var_couplings))] = 1
                else:
                    var_couplings_distr[len(set(var_couplings))] += 1
                var_couplings_var = [v for v in var_couplings if solution.Q_encoding.literal_to_variable[v]!=None]
                if len(set(var_couplings_var)) not in var_couplings_distr_var.keys():
                    var_couplings_distr_var[len(set(var_couplings_var))] = 1
                else:
                    var_couplings_distr_var[len(set(var_couplings_var))] += 1
                var_couplings_aux = [v for v in var_couplings if solution.Q_encoding.literal_to_variable[v]==None]
                if len(set(var_couplings_aux)) not in var_couplings_distr_aux.keys():
                    var_couplings_distr_aux[len(set(var_couplings_aux))] = 1
                else:
                    var_couplings_distr_aux[len(set(var_couplings_aux))] += 1
            print('Distr. couplings for each varaible: ', dict(sorted(var_couplings_distr.items())), '   Mean: ', sum(key * value for key, value in var_couplings_distr.items())/sum(var_couplings_distr.values()), '\n   Between SAT_vars: ', dict(sorted(var_couplings_distr_var.items())), '   Mean: ', sum(key * value for key, value in var_couplings_distr_var.items())/sum(var_couplings_distr_var.values()), '\n   SAT_var with auxiliars: ', dict(sorted(var_couplings_distr_aux.items())), '   Mean: ', sum(key * value for key, value in var_couplings_distr_aux.items())/sum(var_couplings_distr_aux.values()))
    
            chain_len_distr_var = {}
            chain_len_distr_aux = {}
            for var, chain in solution.embedding.items():
                if solution.Q_encoding.literal_to_variable[var]!=None:
                    if len(chain) not in chain_len_distr_var.keys():
                        chain_len_distr_var[len(chain)] = 1
                    else:
                        chain_len_distr_var[len(chain)] += 1
                else:
                    if len(chain) not in chain_len_distr_aux.keys():
                        chain_len_distr_aux[len(chain)] = 1
                    else:
                        chain_len_distr_aux[len(chain)] += 1
            print('Chain length distribution (DWave_embedding): Mean: ', (sum(key * value for key, value in chain_len_distr_var.items())+sum(key * value for key, value in chain_len_distr_aux.items()))/(sum(chain_len_distr_var.values())+sum(chain_len_distr_aux.values())), '\n   Variables: ', dict(sorted(chain_len_distr_var.items())), '   Mean: ', sum(key * value for key, value in chain_len_distr_var.items())/sum(chain_len_distr_var.values()), '\n   Auxiliar: ', dict(sorted(chain_len_distr_aux.items())), '   Mean: ', sum(key * value for key, value in chain_len_distr_aux.items())/sum(chain_len_distr_aux.values()))
    
            print('Solution distribution: ', dict(sorted(solution.o.items())), '   Mean: ', sum(key * value for key, value in solution.o.items())/sum(solution.o.values()))

# %%
