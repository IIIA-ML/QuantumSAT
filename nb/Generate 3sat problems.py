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
from generation import generate_3sat
import random
import numpy as np
from pathlib import Path

# %%
# Generate problems for Fig.1 in Nusslein

random.seed(1345)
n_vars = np.arange(20,401,20)
dir = "../../exp/e1/problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)
for vars in n_vars:
    p = generate_3sat(vars, ratio=4.2)
    with open(p_dir / ("p"+str(vars)+".cnf"),"w") as f:
        f.write(p)

        
    
    
    
    

# %%
# Generate problems for Fig.2 in Nusslein


random.seed(13435)
n_vars = np.arange(15,28,3)
num_instances = 20
dir = "../../exp/e2/problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)
for vars in n_vars:
    for i in range(num_instances):
        p = generate_3sat(vars, ratio=4.2)
        with open(p_dir / ("p"+str(vars)+"-"+str(i)+".cnf"),"w") as f:
            f.write(p)


# %%





# %%
# Generate problems for Table.3 in Nusslein


random.seed(178)
n_vars=[5,10,12]
num_instances = 20
dir = "../../exp/e3/problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)
for vars in n_vars:
    for i in range(num_instances):
        p = generate_3sat(vars, ratio=4.2)
        with open(p_dir / ("p"+str(vars)+"-"+str(i)+".cnf"), "w") as f:
            f.write(p)

# %%