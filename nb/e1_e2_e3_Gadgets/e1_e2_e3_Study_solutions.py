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


# %% [markdown]
# ## Study results

# %%
def compare_with_exact(exp, g, vars, i, o):
    if exp=="e1":
        file_name=f"../../exp/e1/Akmaxsat_(exact)/p{vars}.txt"
    if exp=="e2":
        file_name=f"../../exp/e2/Akmaxsat_(exact)/p{vars}-{i}.txt"
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
plt.subplot(3, 2, 1)
for key, values in n_couplings.items():
    plt.plot(n_vars, values, label=key)
plt.xlabel('# Variables (V)')
plt.ylabel('# Non-zero Couplings in Q')
#plt.ylim(0,14000)
#plt.yticks(np.arange(0,14001,2000))
plt.legend()

plt.subplot(3, 2, 2)
for key, values in o.items():
    plt.plot(n_vars, values, label=key)
plt.xlabel('# Variables (V)')
plt.ylabel('Optim SA - Optim exact (o)')
plt.legend()

plt.subplot(3, 2, 3)
for key, values in log_q.items():
    plt.plot(n_vars, values, label=key)
plt.xlabel('# Variables (V)')
plt.ylabel('# Logical qubits')
plt.legend()

plt.subplot(3, 2, 4)
for key, values in t.items():
    plt.plot(n_vars, values, label=key)
plt.xlabel('# Variables (V)')
plt.ylabel('Time (s)')
plt.legend()

plt.subplot(3, 2, 5)
for key, values in o.items():
    if key=="regular_like" or key=="tree_like":
        plt.plot(n_vars, values, label=key)
plt.xlabel('# Variables (V)')
plt.ylabel('Optim SA - Optim exact (o)')
plt.legend()

plt.gcf().set_size_inches(10, 15)  # Ancho x Alto
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
                    list_chain_len.append(np.median(list_chain_len_no_mean))
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
                list_log_q.append(np.median(list_log_q_no_mean))
    
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
    plt.plot(n_vars, values, label=key)
plt.xlabel('# Variables (V)')
plt.ylabel('# Logical Qubits')
plt.legend()

plt.subplot(2, 2, 2)
for key, values in phys_q.items():
    if key!="choi":
        means = [np.median(arr) for arr in values]
        q1 = [np.percentile(arr, 25) for arr in values]
        q3 = [np.percentile(arr, 75) for arr in values]
        plt.plot(n_vars, means, label=key, marker='o', linestyle='--')
        plt.fill_between(n_vars, q1, q3, alpha=0.3)
plt.xlabel('# Variables (V)')
plt.ylabel('# Physical Qubits')
plt.legend()
plt.grid(linestyle='--')

plt.subplot(2, 2, 3)
for key, values in chain_len.items():
    plt.plot(n_vars, values, label=key)
plt.xlabel('# Variables (V)')
plt.ylabel('Maximum chain length')
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
n_vars=[5,10,12]
iterations=20
for g in gadgets:
    for vars in n_vars:
        for i in range(iterations):
            with open(f"../../exp/e3/DWave/{g}/p{vars}-{i}_{g}.txt"

# %%
