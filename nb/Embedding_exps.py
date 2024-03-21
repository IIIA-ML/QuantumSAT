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
num_vars = 150
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
dwave_token = ""

# %%
instance = CJ2.CJ2("../exp/embedding/problems/p"+str(num_vars)+".cnf")
instance.fillQ()

# %%
embedding = utils.get_embedding(instance.Q, token=dwave_token, random_seed=10)

# %%
embedding

# %%

# %%

# %%
