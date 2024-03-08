# QuantumSAT

This repo contains code for SAT solvers on quantum computers.

## Installation / Starting a work session

Currently only Linux is supported (Windows users can use WSL)

1. Clone the repo
2. source bin/init-local.sh
3. dwave setup
4. You can launch jupyter from console.


After 3.: Answer 'Y' to All and then paste the Solver API token (it is in your Dashboard of D-Wave Leap)

## Directories

* src: keeps the python sources
* bin: executables
* nb: Jupyter notebooks
* exp: Stores the experiments

## Coding and git usage etiquette

* Avoid storing large files or a large number of files in git. 
* Use `jupytext` to pair your notebooks (the repo is configured so that any notebook in the `nb` directory is automatically paired) and commit only the .py files. No .ipynb file should be in the repo.

# Implementing 3-SAT Gadgets for Quantum Annealers with random instances
In order to verify the results shown in the paper, the following steps must be followed:

1. (temp) Checkout to the "development" branch.
2. Execute the notebook Solve.ipynb located in the nb folder. This notebook is organized into four sections. The initial section focuses on generating various 3-SAT instances employed in the different experiments. Specifically, it creates three distinct types of problems and stores them as CNF files within folders labeled e1, e2, and e3, respectively, all situated within the overarching directory named exp. The following cells in the notebook compute the number of non-zero couplings of the QUBO matrix using instances in e1, the number of physical qubits using instances in e2, and solve the problems in e3 using DWave and the exact solver MaxSatZ. The results for all these experiments are stored within the e1, e2, and e3 folders, respectively.
3. After having executed the Solve.ipynb notebook, it is necessary to execute the cells in the Plots.ipynb notebook to replicate the figures and tables presented in the paper. This notebook is also located inside the nb folder. The different cells read out the experiment files stored previously and output Fig.1, Table 1 and Table 2 from the paper, among others. 

