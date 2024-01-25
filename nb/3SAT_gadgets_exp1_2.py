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
import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt
import random
import dimod
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
import dwave.inspector
from pyqubo import Spin
import neal
import time
import numpy as np
import pandas as pd
from greedy import SteepestDescentSolver
from pyvis.network import Network
import os
import sys
from pathlib import Path

sys.path.append('../src')
import generation
import Nuesslein

# %% [markdown]
# ## Generate the 3SAT problems

# %%
secure = input("Do you really want to generate all the problems? (y/n)")
if secure == 'y' or secure == 'Y':
    # !python3 Generate\ 3sat\ problems.py


# %% [markdown]
# ## Gadgets

# %% [markdown]
# ### Parse a problem and get the clauses

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
# ### Gadget kSAT to 3SAT

# %%
def ksat_3sat(clauses, b):
    new_clauses=[]
    for clause in clauses:
        if len(clause)!=3:
            new_clauses.append(np.array([clause[0],clause[1],b]))
            if len(clause)>3:
                for i in range(2, len(clause)-2):
                    b+=1
                    new_clauses.append(np.array([-(b-1),clause[i],b]))
            new_clauses.append(np.array([-b,clause[len(clause)-2], clause[len(clause)-1]]))
                
        else:
            new_clauses.append(np.array(clause))
    return new_clauses, b+1


# %% [markdown]
# ### (2,2)-gadget from Tseitin and (3,9/2)-gadget with E=1 (Carlos/Jordi paper introduction)

# %%
def tseitin(clauses,b): #Al solucionar 3SAT creamos una nueva variable extra
    h={}
    J={}

    clauses, b = ksat_3sat(clauses, b)
    b_init=b
    
    for clause in clauses:
        sign = [1 if x >= 0 else -1 for x in sorted(clause, key=abs)]
        var = [abs(x) for x in sorted(clause, key=abs)]
        
        for v in var:
            if v not in h.keys():
                h[v]=0.0
        h[b]=0.0
        if (var[0],var[1]) not in J.keys(): J[var[0],var[1]]=0.0
        if (var[1],var[2]) not in J.keys(): J[var[1],var[2]]=0.0
        if (var[0],var[2]) not in J.keys(): J[var[0],var[2]]=0.0
        J[var[0],b]=0.0
        J[var[1],b]=0.0
        J[var[2],b]=0.0
        
        #h y J debido a x1 ∨ x2 ↔ b
        h[var[0]]+= (0.5 if sign[0]>0 else -0.5)
        h[var[1]]+= (0.5 if sign[1]>0 else -0.5)
        h[b]-=1
    
        J[var[0],var[1]]+= (0.5 if ((sign[0]>0 and sign[1]>0) or (sign[0]<0 and sign[1]<0)) else -0.5)
        J[var[0],b]-= (1 if sign[0]>0 else -1)
        J[var[1],b]-= (1 if sign[1]>0 else -1)
        
        #h y J debido a b ∨ x3
        h[b]-=0.5
        h[var[2]]-= (0.5 if sign[2]>0 else -0.5)
        J[var[2],b]+= (0.5 if sign[2]>0 else -0.5)
        
        b+=1
    return h,J,(b-1)-b_init


# %%
def three_sat_max2xor(clauses, b):
    h={}
    J={}

    clauses, b = ksat_3sat(clauses, b)
    b_init=b
    
    for clause in clauses:
        sign = [1 if x >= 0 else -1 for x in sorted(clause, key=abs)]
        var = [abs(x) for x in sorted(clause, key=abs)]
        
        for v in var:
            if v not in h.keys():
                h[v]=0.0
        h[b]=0.0
        if (var[0],var[1]) not in J.keys(): J[var[0],var[1]]=0.0
        if (var[1],var[2]) not in J.keys(): J[var[1],var[2]]=0.0
        if (var[0],var[2]) not in J.keys(): J[var[0],var[2]]=0.0
        J[var[0],b]=0.0
        J[var[1],b]=0.0
        J[var[2],b]=0.0
        
        """h[var[0]]-=0.5
        h[var[2]]-=0.5
        J[var[0],var[2]]+=0.5
        h[var[0]]-=0.5
        h[b]+=0.5
        J[var[0],b]-=0.5
        h[var[2]]-=0.5
        h[b]+=0.5
        J[var[2],b]-=0.5
        h[var[1]]-=1 #MIRAR ESTA PARTE!!!!!!!
        h[b]-=1
        J[var[1],b]+=1"""
        h[b]-=0.5
        h[var[1]]-=sign[1]*0.5
        J[var[0],var[2]]+=sign[0]*sign[2]*0.5
        J[var[0],b]-=sign[0]*0.5
        J[var[2],b]-=sign[2]*0.5
        J[var[1],b]+=sign[1]*0.5

        b+=1

    return h,J,(b-1)-b_init    


# %%
def regular_like(clauses, b): #From ksat to Max2XOR
    h={}
    J={}

    b_init=b
    
    for clause in clauses:
        sign = [1 if x >= 0 else -1 for x in sorted(clause, key=abs)]
        var = [abs(x) for x in sorted(clause, key=abs)]

        v_aux = list(range(b,b+len(var)-1))
        
        for i in range(len(var)):
            """if var[i] not in h.keys():
                h[var[i]]=0.0"""
            for j in range(i+1, len(var)):
                if (var[i],var[j]) not in J.keys():
                    J[var[i],var[j]]=0.0
            #k-2 Auxuliary variables
            for aux in range(b, b+len(var)-1):
                J[var[i],aux]=0.0
                if aux!=b+len(var)-2:
                    J[aux,aux+1]=0.0

        k=0
        J[var[0],var[1]]+=sign[0]*sign[1]*0.5
        J[var[0],v_aux[k]]-=sign[0]*0.5
        J[var[1],v_aux[k]]-=sign[1]*0.5
        if len(var)>2:
            for i in range(2,len(var)):
                k+=1
                J[var[i],v_aux[k-1]]+=sign[i]*0.5
                J[v_aux[k-1],v_aux[k]]-=0.5
                J[var[i],v_aux[k]]-=sign[i]*0.5

        b+=(len(var)-2)
    return h,J,(b-1)-b_init


# %% [markdown]
# ## Solving (Simulated annealing)

# %%
def count_unsatisfied_clauses(assignment, clauses):
    count = 0
    for clause in clauses:
        for literal in clause:
            if (literal > 0 and assignment[literal] == 1) or (literal < 0 and assignment[-literal] == -1):
                count += 1
                break
    return str(len(clauses)-count)


# %%
def solving_SAnnealing(h,J,aux,file_name,clauses,non_zero_couplings):
    import time
    ising_model = dimod.BinaryQuadraticModel(h, J, 0.0, dimod.Vartype.SPIN)
    sampler = neal.SimulatedAnnealingSampler()
    start_time = time.time()
    response = sampler.sample(ising_model, num_reads=1000)
    end_time = time.time()
    
    with open(file_name, 'w') as archivo:
        archivo.write("Time: "+str(round(end_time-start_time, 2))+" s.\n")
        archivo.write("Auxiliary variables: "+str(aux)+"\n")
        archivo.write("Non zero couplings: "+str(non_zero_couplings)+"\n")
        archivo.write("o "+count_unsatisfied_clauses(response.first.sample, clauses)+"\n")
        archivo.write("e "+str(response.first.energy)+"\n")
        archivo.write("v "+str(response.first.sample)+"\n")
        archivo.write(str(response))
    #return response


# %% [markdown]
# #### For problems in files

# %%
secure = input("Do you really want to solve all the problems? (y/n)")
if secure == 'y' or secure == 'Y':
    gadgets_jc=["tseitin", "three_sat_max2xor", "regular_like"]
    gadgets=["chancellor","nuesslein1","nuesslein2","choi"]
    n_vars = np.arange(20,401,20)

    for vars in n_vars:
        dir = "../exp/e1/problems"
        p_dir = Path(dir)
        p_dir.mkdir(parents=True, exist_ok=True)
        file_path = p_dir / ("p"+str(vars)+".cnf")

        clauses, b = parsear_cnf_file(file_path)
        
        for g in gadgets_jc:
            dir = f"../exp/e1/Simulated_Annealing/{g}"
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")

            h, J, aux = globals()[g](clauses, b)
            non_zero_couplings = Nuesslein.get_n_couplings(dimod.utilities.ising_to_qubo(h, J, offset=0.0)[0])
            solving_SAnnealing(h,J,aux,file_name,clauses,non_zero_couplings)
            
        
        for g in gadgets:
            dir = f"../exp/e1/Simulated_Annealing/{g}"
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
            
            module = __import__("Nuesslein")
            instance = getattr(module, g)(clauses, b-1)
            time, opt_pos, response, non_zero_couplings = instance.solve()
            with open(file_name, 'w') as archivo:
                archivo.write("Time: "+str(time)+" s.\n")
                archivo.write("Auxiliary variables: "+"\n") #FALTA PONER EL VALOR
                archivo.write("Non zero couplings: "+str(non_zero_couplings)+"\n")
                archivo.write("o "+str(len(clauses)-opt_pos)+"\n")
                archivo.write("e "+str(response.first.energy)+"\n")
                archivo.write("v "+str(response.first.sample)+"\n")
                archivo.write(str(response))

# %% [markdown]
# ## Study results

# %%
gadgets=["tseitin", "three_sat_max2xor", "regular_like","chancellor","nuesslein1","nuesslein2","choi"]
n_vars = np.arange(20,401,20)

n_couplings = {}
for g in gadgets:
    list=[]
    for vars in n_vars:
        dir = f"../exp/e1/Simulated_Annealing/{g}"
        p_dir = Path(dir)
        file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
    
        with open(file_name, 'r') as archivo:
            l=archivo.readlines()
            list.append(int(l[2].split()[3]))
        n_couplings[g]=list

# %%
for key, values in n_couplings.items():
    plt.plot(n_vars, values, label=key)

plt.xlabel('# Variables (V)')
plt.ylabel('# Non-zero Couplings in Q')
plt.ylim(0,14000)
plt.yticks(np.arange(0,14001,2000))
plt.legend()
plt.show()

# %%
gadgets=["tseitin", "three_sat_max2xor", "regular_like","chancellor","nuesslein1","nuesslein2","choi"]
n_vars = np.arange(20,401,20)

o = {}
for g in gadgets:
    list=[]
    for vars in n_vars:
        dir = f"../exp/e1/Simulated_Annealing/{g}"
        p_dir = Path(dir)
        file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
    
        with open(file_name, 'r') as archivo:
            l=archivo.readlines()
            list.append(int(l[3].split()[1]))
        o[g]=list

# %%
for key, values in o.items():
    plt.plot(n_vars, values, label=key)

plt.xlabel('# Variables (V)')
plt.ylabel('Clauses not satisfied (o)')
plt.legend()
plt.show()

# %%
