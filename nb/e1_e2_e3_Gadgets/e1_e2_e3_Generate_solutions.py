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
import random
import dimod
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
import dwave.inspector
from pyqubo import Spin
import neal
import time
import numpy as np
import os
import sys
from pathlib import Path
from minorminer import find_embedding
import math

sys.path.append('../../src')
import generation
import Nuesslein
import Gadgets_Carlos_Jordi

# %% [markdown]
# ## Generate the 3SAT problems

# %%
secure = input("Do you really want to generate all the problems? (y/n)")
if secure == 'y' or secure == 'Y':
    # !python3 ../Generate\ 3sat\ problems.py


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
    log_qubits = len(response.first.sample.keys())
    
    with open(file_name, 'w') as archivo:
        archivo.write("Time: "+str(round(end_time-start_time, 2))+" s.\n")
        archivo.write("Logical qubits: "+str(log_qubits)+"\n")
        archivo.write("Couplings: "+str([keys for keys, values in J.items() if values!=0.0])+"\n")
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
    gadgets_jc=["tseitin", "three_sat_max2xor", "regular_like", "tree_like", "clique"]
    gadgets=["chancellor","nuesslein1","nuesslein2","choi"]

    #experiment 1
    #n_vars = np.arange(20,401,20)
    #experiment 2
    n_vars = np.arange(15,28,3)
    num_instances = 20
    
    for vars in n_vars:
        #experiment 2
        for i in range(num_instances):
            dir = "../../exp/e2/problems" #Change when experiment change
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            #experiment 1
            #file_path = p_dir / ("p"+str(vars)+".cnf")
            #experiment 2
            file_path = p_dir / ("p"+str(vars)+"-"+str(i)+".cnf")
    
            clauses, b = parsear_cnf_file(file_path)
            
            for g in gadgets_jc:
                dir = f"../../exp/e2/Simulated_Annealing/{g}" #Change when experiment change
                p_dir = Path(dir)
                p_dir.mkdir(parents=True, exist_ok=True)
                #experiment 1
                #file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
                #experiment 2
                file_name = p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt")
    
                h, J, aux = getattr(Gadgets_Carlos_Jordi, g)(clauses, b)
                non_zero_couplings = Nuesslein.get_n_couplings(dimod.utilities.ising_to_qubo(h, J, offset=0.0)[0])
                solving_SAnnealing(h,J,aux,file_name,clauses,non_zero_couplings)
                
            
            for g in gadgets:
                dir = f"../../exp/e2/Simulated_Annealing/{g}" #Change when experiment change
                p_dir = Path(dir)
                p_dir.mkdir(parents=True, exist_ok=True)
                #experiment 1
                #file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
                #experiment 2
                file_name = p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt")
                
                module = __import__("Nuesslein")
                instance = getattr(module, g)(clauses, b-1)
                time, opt_pos, response, non_zero_couplings, log_qubits = instance.solve()
                with open(file_name, 'w') as archivo:
                    archivo.write("Time: "+str(time)+" s.\n")
                    archivo.write("Logical qubits: "+str(log_qubits)+"\n")
                    archivo.write("Couplings: "+str([keys for keys, values in instance.Q.items() if (values!=0.0 and keys[0]!=keys[1])])+"\n")
                    archivo.write("Non zero couplings: "+str(non_zero_couplings)+"\n")
                    archivo.write("o "+str(len(clauses)-opt_pos)+"\n")
                    archivo.write("e "+str(response.first.energy)+"\n")
                    archivo.write("v "+str(response.first.sample)+"\n")
                    archivo.write(str(response))


# %% [markdown]
# ### Exact solver (Akmaxsat)

# %%
def exact_ak(file_path, file_name):
    sample_set = !./../../Maxsatz/maxsatz {file_path}

    with open(file_name, 'w') as archivo:
        for line in sample_set:
            archivo.write(line + '\n')


# %% [markdown]
# ### Approx solver (NuWLS)

# %%
def NuWLS(file_path, file_name):
    command=f""
    sample_set = !./../../NuWLS-c-main/bin/starexec_run_default-runsolver {file_path}

    with open(file_name, 'w') as archivo:
        for line in sample_set:
            archivo.write(line + '\n')


# %%
secure = input("Do you really want to solve all the problems? (y/n)")
if secure == 'y' or secure == 'Y':
    solver = input("What solver you want? (Akmaxsat/NuWLS)")

    #experiment 1
    #n_vars = np.arange(380,401,20) #np.arange(20,401,20)
    #experiment 2
    #n_vars = np.arange(15,28,3)
    #num_instances = 20
    #experiment 3
    n_vars = [5,10,12]
    num_instances = 20
    
    for vars in n_vars:
        #experiment 2/3
        for i in range(num_instances):
            dir = "../../exp/e3/problems" #Change when experiment change
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            #experiment 1
            #file_path = p_dir / ("p"+str(vars)+".cnf")
            #experiment 2/3
            file_path = p_dir / ("p"+str(vars)+"-"+str(i)+".cnf")                

            if solver == "Akmaxsat":
                dir = f"../../exp/e3/Akmaxsat_(exact)" #Change when experiment change
            if solver == "NuWLS":
                dir = f"../../exp/e3/NuWLS_(approx)" #Change when experiment change
            p_dir = Path(dir)
            p_dir.mkdir(parents=True, exist_ok=True)
            #experiment 1
            #file_name = p_dir / ("p"+str(vars)+".txt")
            #experiment 2/3
            file_name = p_dir / ("p"+str(vars)+"-"+str(i)+".txt")

            if solver == "Akmaxsat":
                exact_ak(file_path, file_name)
            if solver == "NuWLS":
                NuWLS(file_path, file_name)


# %% [markdown]
# ## Generate embedding

# %%
def embed(exp):
    gadgets=["tseitin", "three_sat_max2xor", "regular_like", "tree_like", "clique", "chancellor","nuesslein1","nuesslein2","choi"]
    embedding = {}
    
    for g in gadgets:
        
        if exp == "e1":
            n_vars = np.arange(20,401,20)
            for vars in n_vars:
                dir = f"../../exp/e1/Simulated_Annealing/{g}"
                p_dir = Path(dir)
                file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
                try:
                    with open(file_name, 'r') as archivo:
                        l = archivo.readlines()
                        embed = find_embedding(eval(l[2][11:]), DWaveSampler().edgelist, random_seed=10)
                except Exception as e:
                    print(f"Error while embedding for {g}, vars={vars}, instance={i}: {e}")
                    embed = None
                dir = f"../../exp/e1/Embedding/{g}"
                p_dir = Path(dir)
                p_dir.mkdir(parents=True, exist_ok=True)
                with open(p_dir / ("p"+str(vars)+"_"+g+".txt"), "w") as archivo:
                    archivo.write(str(embed))

        if exp == "e2":
            n_vars = np.arange(15,28,3)
            num_instances = 20
            for vars in n_vars:
                for i in range(num_instances):
                    dir = f"../../exp/e2/Simulated_Annealing/{g}"
                    p_dir = Path(dir)
                    file_name = p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt")
                    try:
                        with open(file_name, 'r') as archivo:
                            l = archivo.readlines()
                            embed = find_embedding(eval(l[2][11:]), DWaveSampler().edgelist, random_seed=10)
                    except Exception as e:
                        print(f"Error while embedding for {g}, vars={vars}, instance={i}: {e}")
                        embed = None
                    dir = f"../../exp/e2/Embedding/{g}"
                    p_dir = Path(dir)
                    p_dir.mkdir(parents=True, exist_ok=True)
                    with open(p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt"), "w") as archivo:
                        archivo.write(str(embed))

        if exp == "e3":
            n_vars = [5,10,12]
            num_instances = 20
            for vars in n_vars:
                for i in range(num_instances):
                    clauses, b = parsear_cnf_file(f"../../exp/e3/problems/p{vars}-{i}.cnf")
                    if g in ["tseitin", "three_sat_max2xor", "regular_like", "tree_like", "clique"]:
                        h, J, aux = getattr(Gadgets_Carlos_Jordi, g)(clauses, b)
                        embed = find_embedding([keys for keys, values in J.items() if values!=0.0], DWaveSampler().edgelist, random_seed=10)
                    else:
                        instance = getattr(Nuesslein, g)(clauses, b-1)
                        instance.fillQ()
                        embed=find_embedding([keys for keys, values in instance.Q.items() if (values!=0.0 and keys[0]!=keys[1])], DWaveSampler().edgelist, random_seed=10)
                    dir = f"../../exp/e3/Embedding/{g}"
                    p_dir = Path(dir)
                    p_dir.mkdir(parents=True, exist_ok=True)
                    with open(p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt"), "w") as archivo:
                        archivo.write(str(embed))


# %%
secure = input("Do you really want to embed all the problems? (y/n)")
if secure == 'y' or secure == 'Y':
    exp = input("Which experiment you want to embed? (e1/e2/e3)")
    embed(exp)


# %% [markdown]
# ### DWave (e3)

# %%
def solving_DWave(h,J_old,embedding,file_name,clauses,num_reads,anneal_time):
    J={key: value for key,value in J_old.items() if value!=0.0}
    
    sampler = FixedEmbeddingComposite(DWaveSampler(token='DEV-291d80af600d6eb433a8019c579070ba37436e9a'), embedding=embedding)
    sample_set = sampler.sample_ising(h, J, num_reads=num_reads, annealing_time=anneal_time, return_embedding=True, reduce_intersample_correlation=True)

    qpu_access_time = sample_set.info['timing']['qpu_access_time']
    #Tiempo usado por la QPU para resolver el problema de forma activa
    qpu_sampling_time = sample_set.info['timing']['qpu_sampling_time']

    with open(file_name, 'w') as archivo:
        archivo.write("QPU access time:"+str(qpu_access_time*10**(-6))+"s."+" ; QPU sampling time:"+str(qpu_sampling_time*10**(-6))+"s.\n")
        archivo.write("o "+count_unsatisfied_clauses(sample_set.first.sample, clauses)+"\n")
        archivo.write("e "+str(sample_set.first.energy)+"\n")
        archivo.write("v "+str(sample_set.first.sample)+"\n")
        archivo.write(str(sample_set))


# %% [markdown]
# #### Study which num_read and anneal_time to use

# %%
secure = input("Do you really want to solve all the problems? (y/n)")
if secure=="y" or secure=="Y":
    gadgets=["tseitin"]
    n_vars=[10]
    iterations=1
    reads=[20,50,100]
    anneal_t=[20,50,100,200]
    for vars in n_vars:
        for i in range(iterations):
            file_path=f"../../exp/e3/problems/p{vars}-{i}.cnf"
            for g in gadgets:
                clauses, b = parsear_cnf_file(file_path)
                h, J, aux = getattr(Gadgets_Carlos_Jordi, g)(clauses, b)
                with open(f"../../exp/e3/Embedding/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
                    embedding = eval(archivo.readlines()[0])
                for num_reads in reads:
                    for anneal_time in anneal_t:
                        dir = f"../../exp/e3/DWave/Test_reads_anneal_t/{g}"
                        p_dir = Path(dir)
                        p_dir.mkdir(parents=True, exist_ok=True)
                        file_name= p_dir / ("p"+str(vars)+"-"+str(i)+"_reads_"+str(num_reads)+"_anneal_t_"+str(anneal_time)+".txt")
                        solving_DWave(h,J,embedding,file_name,clauses,num_reads,anneal_time)

# %% [markdown]
# #### Solve all problems (e3)

# %%
secure = input("Do you really want to solve all the problems? (y/n)")
if secure=="y" or secure=="Y":
    gadgets=["tseitin"] #["tseitin", "three_sat_max2xor", "regular_like", "tree_like"]
    n_vars=[5] #[5,10,12]
    iterations=1 #20
    for vars in n_vars:
        for i in range(iterations):
            file_path=f"../../exp/e3/problems/p{vars}-{i}.cnf"
            for g in gadgets:
                dir = f"../../exp/e3/DWave/{g}"
                p_dir = Path(dir)
                p_dir.mkdir(parents=True, exist_ok=True)
                file_name= p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt")
                clauses, b = parsear_cnf_file(file_path)
                if g in ["tseitin", "three_sat_max2xor", "regular_like", "tree_like", "choi"]:
                    h, J, aux = getattr(Gadgets_Carlos_Jordi, g)(clauses, b)
                else:
                    instance = getattr(Nuesslein, g)(clauses, b-1)
                    instance.fillQ()
                    h, J, e = dimod.qubo_to_ising(instance.Q, 0.0)
                with open(f"../../exp/e3/Embedding/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
                    embedding = eval(archivo.readlines()[0])
                solving_DWave(h,J,embedding,file_name,clauses,100,150)

# %%
