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
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd


# %% [markdown]
# ## Study results

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
def compare_with_exact(exp, g, vars, i, o):
    if exp=="e1":
        file_name=f"../../exp/e1/Akmaxsat_(exact)/p{vars}.txt"
    if exp=="e2":
        file_name=f"../../exp/e2/Akmaxsat_(exact)/p{vars}-{i}.txt"
    if exp=="e3":
        file_name=f"../../exp/e3/Akmaxsat_(exact)/p{vars}-{i}.txt"
    with open(file_name, "r") as archivo:
        lines=archivo.readlines()
        for l in reversed(lines):
            if l.split()[0]=='o':
                optim=int(l.split()[1])
                break

    return o-optim


# %%
gadgets=["tseitin", "three_sat_max2xor", "regular_like", "tree_like", "clique", "chancellor","nuesslein1","nuesslein2","choi"]
exp = input("Which experiment you want to study? (e1/e2)")

n_couplings = {}
o = {}
log_q = {}
t = {}

for g in gadgets:
    list_n_couplings=[]
    list_o=[]
    list_log_q=[]
    list_t=[]
    if exp == "e1":
        n_vars = np.arange(20,361,20) #np.arange(20,401,20)
        for vars in n_vars:
            dir = f"../../exp/e1/Simulated_Annealing/{g}"
            p_dir = Path(dir)
            file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
            with open(file_name, 'r') as archivo:
                l=archivo.readlines()
                list_n_couplings.append(int(l[3].split()[3]))
                list_o.append(compare_with_exact(exp,g,vars,None,int(l[4].split()[1])))
                list_log_q.append(len(set([keys for keys, values in eval(l[6][2:]).items()])))
                list_t.append(float(l[0].split()[1]))
        
    if exp == "e2":
        n_vars = np.arange(15,28,3)
        num_instances = 20
        for vars in n_vars:
            list_n_couplings_no_mean = []
            list_o_no_mean = []
            list_log_q_no_mean = []
            list_t_no_mean = []
            for i in range(num_instances):
                dir = f"../../exp/e2/Simulated_Annealing/{g}"
                p_dir = Path(dir)
                file_name = p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt")
                with open(file_name, 'r') as archivo:
                    l=archivo.readlines()
                    list_n_couplings_no_mean.append(int(l[3].split()[3]))
                    #list_o_no_mean.append(int(l[4].split()[1]))
                    list_o_no_mean.append(compare_with_exact(exp,g,vars,i,int(l[4].split()[1])))
                    list_log_q_no_mean.append(len(set([keys for keys, values in eval(l[6][2:]).items()])))
                    list_t_no_mean.append(float(l[0].split()[1]))
            list_n_couplings.append(sum(list_n_couplings_no_mean)/len(list_n_couplings_no_mean))
            list_o.append(sum(list_o_no_mean)/len(list_o_no_mean))
            list_log_q.append(sum(list_log_q_no_mean)/len(list_log_q_no_mean))
            list_t.append(sum(list_t_no_mean)/len(list_t_no_mean))
    
    n_couplings[g]=list_n_couplings
    o[g]=list_o
    log_q[g]=list_log_q
    t[g]=list_t

# %%
plt.subplot(2, 2, 1)
for key, values in n_couplings.items():
    if key=="nuesslein1" or key=="nuesslein2" or key=="tseitin" or key=="three_sat_max2xor":
        plt.plot(n_vars, values, label=key, marker='+', markeredgewidth=1, linestyle='--')
plt.xlabel('# Variables (V)')
plt.ylabel('# Non-zero Couplings in Q')
plt.grid(linestyle='--')
plt.legend()

plt.subplot(2, 2, 2)
for key, values in o.items():
    if key=="nuesslein1" or key=="nuesslein2" or key=="tseitin" or key=="three_sat_max2xor":
        plt.plot(n_vars, values, label=key, marker='+', markeredgewidth=1, linestyle='-')
plt.xlabel('# Variables (V)')
plt.ylabel('Optim SA - Optim exact (o)')
plt.grid(linestyle='--')
plt.legend()

plt.subplot(2, 2, 3)
for key, values in log_q.items():
    if key=="nuesslein1" or key=="nuesslein2" or key=="tseitin" or key=="three_sat_max2xor":
        plt.plot(n_vars, values, label=key, marker='+', markeredgewidth=1, linestyle='-')
plt.xlabel('# Variables (V)')
plt.ylabel('# Logical qubits')
plt.grid(linestyle='--')
plt.legend()

plt.subplot(2, 2, 4)
for key, values in t.items():
    if key=="nuesslein1" or key=="nuesslein2" or key=="tseitin" or key=="three_sat_max2xor":
        plt.plot(n_vars, values, label=key, marker='+', markeredgewidth=1, linestyle='-')
plt.xlabel('# Variables (V)')
plt.ylabel('Time (s)')
plt.grid(linestyle='--')
plt.legend()


plt.gcf().set_size_inches(10, 10)  # Ancho x Alto
plt.show()


# %% [markdown]
# ### Study embedding

# %%
def study_embed(exp):
    with open(f"../../exp/{exp}/Embedding.txt", 'r') as archivo:
        l=archivo.readlines()
        embed=eval(l[0])
        gadgets=["tseitin", "three_sat_max2xor", "regular_like", "tree_like", "clique", "chancellor","nuesslein1","nuesslein2","choi"]
        phys_q = {}
        chain_len = {}
        for g in gadgets:
            list_embedding=[]
            list_phys_q=[]
            list_chain_len=[]
            
            if exp == "e1":
                n_vars = np.arange(20,401,20)
                for vars in n_vars:
                    list_phys_q.append(len(set([value for sublist in embed[g][vars].values() for value in sublist])))
                    list_chain_len.append(max(len(value) for key, value in embed[g][vars].items()))
                phys_q[g]=list_phys_q
                chain_len[g]=list_chain_len

            if exp == "e2":
                n_vars = np.arange(15,28,3)
                num_instances = 20
                for vars in n_vars:
                    list_phys_q_no_mean=[]
                    list_chain_len_no_mean=[]
                    for i in range(num_instances):
                        list_phys_q_no_mean.append(len(set([value for sublist in embed[g][vars][i].values() for value in sublist])))
                        list_chain_len_no_mean.append(max(len(value) for key, value in embed[g][vars][i].items()))
                    list_phys_q.append(list_phys_q_no_mean)
                    list_chain_len.append(list_chain_len_no_mean)
                phys_q[g]=list_phys_q
                chain_len[g]=list_chain_len
    
    for g in gadgets:
        list_log_q=[]
        dir = f"../../exp/{exp}/Simulated_Annealing/{g}"
        if exp == "e1":
            for vars in n_vars:
                p_dir = Path(dir)
                file_name = p_dir / ("p"+str(vars)+"_"+g+".txt")
                with open(file_name, 'r') as archivo:
                    l=archivo.readlines()
                    list_log_q.append(len(set([keys for keys, values in eval(l[6][2:]).items()])))
            
        if exp == "e2":
            for vars in n_vars:
                list_log_q_no_mean = []
                for i in range(num_instances):
                    p_dir = Path(dir)
                    file_name = p_dir / ("p"+str(vars)+"-"+str(i)+"_"+g+".txt")
                    with open(file_name, 'r') as archivo:
                        l=archivo.readlines()
                        list_log_q_no_mean.append(len(set([keys for keys, values in eval(l[6][2:]).items()])))
                list_log_q.append(list_log_q_no_mean)
    
        log_q[g]=list_log_q
    
    return phys_q, chain_len, log_q


# %%
exp = input("Which experiment you want to study? (e1/e2)")
if exp == "e1" or exp == "e2":
    phys_q, chain_len, log_q = study_embed(exp)

# %%
if exp == "e1":
    n_vars = np.arange(20,401,20)
if exp == "e2":
    n_vars = np.arange(15,28,3)
plt.subplot(2, 2, 1)
for key, values in log_q.items():
    if key=="nuesslein1" or key=="nuesslein2" or key=="tseitin" or key=="three_sat_max2xor":
        means = [np.median(arr) for arr in values]
        plt.plot(n_vars, means, label=key, marker='+', markeredgewidth=1, linestyle='--')
plt.xlabel('# Variables (V)')
plt.ylabel('# Logical Qubits')
plt.grid(linestyle='--')
plt.legend()

plt.subplot(2, 2, 2)
for key, values in phys_q.items():
    if key=="nuesslein1" or key=="nuesslein2" or key=="tseitin" or key=="three_sat_max2xor":
        means = [np.median(arr) for arr in values]
        q1 = [np.percentile(arr, 25) for arr in values]
        q3 = [np.percentile(arr, 75) for arr in values]
        plt.plot(n_vars, means, label=key, marker='+', markeredgewidth=1, linestyle='--')
        plt.fill_between(n_vars, q1, q3, alpha=0.3)
plt.xlabel('# Variables (V)')
plt.ylabel('# Physical Qubits')
plt.grid(linestyle='--')
plt.legend()

plt.subplot(2, 2, 3)
for key, values in chain_len.items():
    if key=="nuesslein1" or key=="nuesslein2" or key=="tseitin" or key=="three_sat_max2xor":
        means = [np.median(arr) for arr in values]
        plt.plot(n_vars, means, label=key, marker='+', markeredgewidth=1, linestyle='--')
plt.xlabel('# Variables (V)')
plt.ylabel('Maximum chain length')
plt.grid(linestyle='--')
plt.legend()

plt.gcf().set_size_inches(10, 10)  # Ancho x Alto
plt.show()

# %% [markdown]
# ### Study of the search of num_reads and anneal_time

# %%
gadgets=["tseitin"]
n_vars=[10]
iterations=1
reads=[20,50,100]
anneal_t=[20,50,100,200]
o=[]
for vars in n_vars:
    for i in range(iterations):
        for g in gadgets:
            with open(f"../../exp/e3/Akmaxsat_(exact)/p{vars}-{i}.txt", "r") as exact:
                lines=exact.readlines()
                for l in reversed(lines):
                    if l[0]=="o":
                        o_exact=int(l.split()[1])
                        break
            for num_reads in reads:
                o_reads=[]
                for anneal_time in anneal_t:
                    dir = f"../../exp/e3/DWave/Test_reads_anneal_t/{g}"
                    p_dir = Path(dir)
                    file_name= p_dir / ("p"+str(vars)+"-"+str(i)+"_reads_"+str(num_reads)+"_anneal_t_"+str(anneal_time)+".txt")
                    with open(file_name, "r") as archivo:
                        o_reads.append(o_exact-int(archivo.readlines()[1].split()[1]))
                o.append(o_reads)

plt.imshow(o, cmap='viridis', extent=[min(anneal_t), max(anneal_t), min(reads), max(reads)],
           aspect='auto', origin='lower', interpolation='none', vmin=min(min(o)), vmax=max(max(o)))
plt.xlabel('Annealing time (microsec)')
plt.ylabel('# Reads')
plt.xticks(anneal_t)
plt.yticks(reads)
cbar = plt.colorbar()
cbar.set_label('O_exact-O_found')
plt.title('Optimum exact- Optimum found (tseitin)')

# %% [markdown]
# ### Study results DWave

# %%
gadgets=["tseitin", "three_sat_max2xor", "regular_like", "tree_like", "nuesslein2"]
n_vars=[5,10,12]
iterations=20
head=["(V=5, C=21)", "(V=10, C=42)", "(V=12, C=50)"]
o={}
o_compare={}
k=0
for vars in n_vars:
    o_g=[]
    o_compare_g=[]
    for g in gadgets:
        o_it=[]
        o_compare_it=[]
        for i in range(iterations):
            with open(f"../../exp/e3/DWave/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
                l=archivo.readlines()
                o_it.append(int(vars*4.2)-int(l[1].split()[1]))
                o_compare_it.append(compare_with_exact("e3",g,vars,i,int(l[1].split()[1])))
        o_g.append(np.mean(o_it))
        o_compare_g.append(np.mean(o_compare_it))
    o[head[k]]=o_g
    o_compare[head[k]]=o_compare_g
    k+=1

# %%
print("OPTIMUM FOUND")
pd.DataFrame(o, gadgets)

# %%
print("OPTIMUM FOUND - OPTIMUM EXACT")
pd.DataFrame(o_compare, gadgets)

# %%
gadgets=["nuesslein2"]
n_vars=[10]
iterations=20
head=["(V=10, C=42)"]
o={}
o_compare={}
k=0
for vars in n_vars:
    o_g=[]
    o_compare_g=[]
    for g in gadgets:
        o_it=[]
        o_compare_it=[]
        for i in range(iterations):
            with open(f"../../exp/e3/QUBO_J/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
                l=archivo.readlines()
                o_it.append(int(vars*4.2)-int(l[1].split()[1]))
                o_compare_it.append(compare_with_exact("e3",g,vars,i,int(l[1].split()[1])))
        o_g.append(np.mean(o_it))
        o_compare_g.append(np.mean(o_compare_it))
    o[head[k]]=o_g
    o_compare[head[k]]=o_compare_g
    k+=1

# %%
print("OPTIMUM FOUND")
pd.DataFrame(o, gadgets)

# %%
gadgets=["tseitin", "three_sat_max2xor", "nuesslein2", "nuesslein1"]
n_vars=[5,10,12,20,50]
iterations=20
head=["(V=5, C=21)", "(V=10, C=42)", "(V=12, C=50)", "(V=20, C=84)", "(V=50, C=210)"]
o={}
k=0
for vars in n_vars:
    o_g={}
    o_ver={}
    for g in gadgets:
        o_list={}
        for i in range(iterations):
            clauses, b = parsear_cnf_file(f"../../exp/e3/problems/p{vars}-{i}.cnf")
            sample_set_str=""
            #with open(f"../../exp/e3/DWave/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
            with open(f"../../exp/Dwave_problems/{g}/p{vars}/p{vars}_{i}_dwave.txt", "r") as archivo:
                l=archivo.readlines()
                for r in range(3,len(l)-1):
                    if r==3:
                        sample_set_str+=l[r][55:]
                    elif r==len(l)-2:
                        line=l[r].replace(" ", "")
                        sample_set_str+=line[:-2]
                    else:
                        sample_set_str+=l[r]
                sample_set=eval(sample_set_str)
                for sample in sample_set:
                    assignment={i+1:sample[0][i] for i in range(len(sample[0]))}
                    o_found = int(count_unsatisfied_clauses(assignment, clauses))
                    o_compare = compare_with_exact("e3", g, vars, i, o_found)
                    if o_compare not in o_list.keys():
                        o_list[o_compare]=sample[-2]
                    else:
                        o_list[o_compare]+=sample[-2]
        o_g[g]={key:value/iterations for key, value in o_list.items()} #value/(iterations*num_reads)*100%
    o[head[k]]=o_g
    k+=1

# %%
fig, axs = plt.subplots(3, 2, figsize=(8, 12))

i = 0
for o_vars in o.values():
    j = 1
    for o_g_name, o_g in o_vars.items():
        x_pos = np.arange(len(o_g.keys()))
        axs[i // 2, i % 2].bar(x_pos + j * 0.2, o_g.values(), label=o_g_name, width=0.2)
        j += 1
    axs[i // 2, i % 2].set_title(list(o.keys())[i])
    axs[i // 2, i % 2].set_xlabel('Optimum_found - Optimum exact')
    axs[i // 2, i % 2].set_ylabel('%')
    axs[i // 2, i % 2].set_xlim([0,6])
    axs[i // 2, i % 2].set_xticks(np.arange(6))
    #axs[i // 2, i % 2].set_xticklabels(o_g.keys())
    axs[i // 2, i % 2].legend()
    i += 1

plt.tight_layout()
plt.show()

# %% [markdown]
# ### Study results SA

# %%
exp = input("What problem do you want to study with SA? (e2/e3)")
gadgets=["tseitin", "three_sat_max2xor", "nuesslein2", "nuesslein1"]
if exp=="e2":
    n_vars = np.arange(15,28,3)
    head=[f"(V={v}, C={int(4.2*v)})" for v in n_vars]
if exp=="e3":
    n_vars=[5,10,12,20]
    head=["(V=5, C=21)", "(V=10, C=42)", "(V=12, C=50)", "(V=20, C=84)"]
iterations=20
o={}
o_compare={}
k=0
for vars in n_vars:
    o_g=[]
    o_compare_g=[]
    for g in gadgets:
        o_it=[]
        o_compare_it=[]
        for i in range(iterations):
            with open(f"../../exp/{exp}/Simulated_Annealing/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
                l=archivo.readlines()
                o_it.append(int(vars*4.2)-int(l[4].split()[1]))
                o_compare_it.append(compare_with_exact(exp,g,vars,i,int(l[4].split()[1])))
        o_g.append(np.mean(o_it))
        o_compare_g.append(np.mean(o_compare_it))
    o[head[k]]=o_g
    o_compare[head[k]]=o_compare_g
    k+=1

# %%
print("OPTIMUM FOUND (SA)")
pd.DataFrame(o, gadgets)

# %%
print("OPTIMUM FOUND (SA) - OPTIMUM EXACT")
pd.DataFrame(o_compare, gadgets)

# %%
gadgets=["tseitin", "three_sat_max2xor", "nuesslein2", "nuesslein1"]
n_vars=[5,10,12,20]
iterations=20
head=["(V=5, C=21)", "(V=10, C=42)", "(V=12, C=50)", "(V=20, C=84)"]
o={}
k=0
for vars in n_vars:
    o_g={}
    for g in gadgets:
        o_list={}
        for i in range(iterations):
            clauses, b = parsear_cnf_file(f"../../exp/e3/problems/p{vars}-{i}.cnf")
            sample_set_str=""
            #with open(f"../../exp/e3/DWave/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
            with open(f"../../exp/e3/Simulated_Annealing/{g}/p{vars}-{i}_{g}.txt", "r") as archivo:
                l=archivo.readlines()
                for r in range(7,len(l)-1):
                    if r==7:
                        sample_set_str+=l[r][55:]
                    elif r==len(l)-2:
                        line=l[r].replace(" ", "")
                        sample_set_str+=line[:-2]
                    else:
                        sample_set_str+=l[r]
                sample_set=eval(sample_set_str)
                for sample in sample_set:
                    assignment={i+1:sample[0][i] for i in range(len(sample[0]))}
                    o_found = int(count_unsatisfied_clauses(assignment, clauses))
                    o_compare = compare_with_exact("e3", g, vars, i, o_found)
                    if o_compare not in o_list.keys():
                        o_list[o_compare]=sample[-1]
                    else:
                        o_list[o_compare]+=sample[-1]
        o_g[g]={key:value/iterations for key, value in o_list.items()} #value/(iterations*num_reads)*100%
    o[head[k]]=o_g
    k+=1

# %%
fig, axs = plt.subplots(2, 2, figsize=(8, 8))

i = 0
for o_vars in o.values():
    j = 1
    for o_g_name, o_g in o_vars.items():
        x_pos = np.arange(len(o_g.keys()))
        axs[i // 2, i % 2].bar(x_pos + j * 0.2, o_g.values(), label=o_g_name, width=0.2)
        j += 1
    axs[i // 2, i % 2].set_title(list(o.keys())[i])
    axs[i // 2, i % 2].set_xlabel('Optimum_found - Optimum exact')
    axs[i // 2, i % 2].set_ylabel('%')
    axs[i // 2, i % 2].set_xlim([0,11])
    axs[i // 2, i % 2].set_xticks(np.arange(11))
    #axs[i // 2, i % 2].set_xticklabels(o_g.keys())
    axs[i // 2, i % 2].legend()
    i += 1

plt.tight_layout()
plt.show()

# %%
