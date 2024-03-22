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
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import sys
from pathlib import Path
import time
sys.path.append('../src')

import utils
import CJ2

from dwave.system import DWaveSampler, EmbeddingComposite
from minorminer import find_embedding

# %%

# %% [markdown]
# ### Generate (large enough) problem

# %%
random.seed(901)
num_vars = 5
p = utils.generate_3sat(num_vars, ratio=4.2)

# %%

# %% [markdown]
# #### Save problem in exp/embedding/problems/

# %%
dir = "../exp/embedding/problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)

with open(p_dir / ("p"+str(num_vars)+".cnf"),"w") as f:
    f.write(p)

# %%

# %%

# %% [markdown]
# ### Perform embedding with different seeds
# ##### When solving with DWave, do the solutions depend strongly on the embedding?

# %%
dwave_token = "Your token"

# %%
instance = CJ2.CJ2("../exp/embedding/problems/p"+str(num_vars)+".cnf")
instance.fillQ()

# %%
for i in range(0,10):
    print(f"Finding Embedding #{i}...")
    sampler = EmbeddingComposite(DWaveSampler(token=dwave_token), embedding_parameters={"random_seed":i})
    print(f"Embedding #{i} found.")
    
    response = sampler.sample_qubo(instance.Q, num_reads=100, annealing_time=100, \
                               reduce_intersample_correlation=True, return_embedding=True)

    print(f"Solved with Embedding #{i}. Appending in txt...")    
    with open(p_dir / ("p"+str(num_vars)+"_solutions.txt"),"a") as f:
        f.write("o "+str(utils.count_unsatisfied_clauses(response.first.sample, instance.clauses))+"\n")
        f.write("e "+str(response.first.energy)+"\n")
        f.write("v "+str(response.first.sample)+"\n")
        f.write(str(response.info['embedding_context']))
        f.write(str(response.samples))

# %%
