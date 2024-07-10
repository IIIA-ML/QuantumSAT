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

# %%
random.seed(901)
num_vars = 20

dir = "../exp/eBeyond/Bian_study/Problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)

#p = utils.generate_3sat(num_vars, ratio=4.2)
p = "c generated problem\np cnf 3 4\n1 3 2 0\n3 -1 2 0\n1 3 -2 0\n1 -3 -2 0\n"
#p = "c generated problem\np cnf 3 2\n1 3 -2 0\n-3 1 -2 0"
#p = "c generated problem\np cnf 3 1\n1 3 2 0"
file_name = "p"+str(num_vars)+".cnf"
file_path = f"{dir}/{file_name}"
with open(file_path,"w") as f:
    f.write(p)

# %%
splitter = Splitter.Two_subproblems()
gadget = [Gadget.CJ1(), Gadget.CJ2()]
joiner = [Joiner.SAT_variables(), Joiner.Bian()]
#solver = Solver.D_Wave()
#token = 'DEV-291d80af600d6eb433a8019c579070ba37436e9a'

# %%
solution = SAT_solver.SAT_solver(file_path=file_path, splitter=splitter, gadget=gadget, joiner=joiner)

# %%
solution.Solve()

# %%
solution.response

# %%
solution.subproblems

# %%
import dimod
import numpy as np
import matplotlib.pyplot as plt
bqm = dimod.BinaryQuadraticModel.from_qubo(solution.Q, offset=0).change_vartype("SPIN", True)
num_vars = len(solution.Q_encoding.literal_to_variable.keys())
matrix = np.zeros((num_vars, num_vars))
for i, value in bqm.linear.items():
    matrix[i, i] = value
for (i, j), value in bqm.quadratic.items():
    matrix[j, i] = value
fig, ax = plt.subplots()
ax.axis('tight')
ax.axis('off')

label = sorted(list(solution.Q_encoding.literal_to_variable.keys()))

table = ax.table(cellText=matrix, cellLoc='center', loc='center', colLabels=label, rowLabels=label, colWidths=[0.1]*num_vars)
table.scale(1, 1.5)
plt.title('BQM (Ising)')
plt.show()

# %%
