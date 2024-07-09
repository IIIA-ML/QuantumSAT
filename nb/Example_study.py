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

# %%
random.seed(901)
num_vars = 20

dir = "../exp/eBeyond/Bian_study/Problems"
p_dir = Path(dir)
p_dir.mkdir(parents=True, exist_ok=True)

p = utils.generate_3sat(num_vars, ratio=4.2)
#p = "c generated problem\np cnf 3 4\n1 2 3 0\n-1 2 3 0\n1 -2 3 0\n 1 -2 -3 0\n"
#p = "c generated problem\np cnf 3 1\n1 3 2 0\n3 1 2 0"
file_name = "p"+str(num_vars)+".cnf"
file_path = f"{dir}/{file_name}"
with open(file_path,"w") as f:
    f.write(p)

# %%
splitter = Splitter.Single_clause()
gadget = Gadget.CJ1()
joiner = Joiner.SAT_variables()
solver = Solver.D_Wave()
token = 'DEV-291d80af600d6eb433a8019c579070ba37436e9a'

# %%
solution = main.SAT_solver(file_path=file_path, Splitter=splitter, Gadget=gadget, Solver=solver, Joiner=joiner, token=token, qubit_level=True)

# %%
solution.Solve()

# %%
solution.response

# %%
solution.clauses

# %%
solution.N

# %%
a = {1:1, 2:2, 3:3, 4:None}
set(a.values())

# %%
solution.Q

# %%
solution.Q_encoding.literal_to_variable

# %%
solution.Q_encoding.variable_to_literal

# %%
solution.Q_encoding.variables

# %%
