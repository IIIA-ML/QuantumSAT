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


