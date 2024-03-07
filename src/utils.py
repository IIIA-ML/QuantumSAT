import random
import numpy as np
from pathlib import Path
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
from minorminer import find_embedding
import neal
from pathlib import Path
import subprocess



def generate_3sat(num_vars, ratio=None, num_clauses=None):
    if ratio is not None:
        if num_clauses is not None:
            raise Exception("Clause ratio and number of clauses cannot be jointly specified")
        else:
            num_clauses = int(num_vars * ratio)
            if num_clauses <= 0:
                num_clauses = 1
    elif num_clauses is None:
        raise Exception("Please specify clause ratio or number of clauses")
    
    problem = "c generated problem\n"
    problem += "p cnf " + str(num_vars) + " " + str(num_clauses) + "\n"
    clauses = {}
    for i in range(num_clauses):  
        while True:
            lit1 = random.randint (1, num_vars)
            if random.randint (0, 1) == 1:
                lit1 = -lit1
            while True:
                lit2 = random.randint (1, num_vars)
                if random.randint (0, 1) == 1:
                    lit2 = -lit2
                if lit1 != lit2 and -lit1 != lit2:
                    break
            while True:
                lit3 = random.randint (1, num_vars)
                if random.randint (0, 1) == 1:
                    lit3 = -lit3
                if lit1 != lit3 and -lit1 != lit3 and lit2 != lit3 and -lit2 != lit3:
                    break
            assert lit1 != lit2
            assert lit1 != -lit2
            assert lit1 != lit3
            assert lit1 != -lit3
            assert lit2 != lit3
            assert lit2 != -lit3
            #sort literals
            if abs(lit2) < abs(lit1):
                tmp = lit1
                lit1 = lit2
                lit2 = tmp
            if abs(lit3) < abs(lit1):
                tmp = lit1
                lit1 = lit3
                lit3 = tmp
            if abs(lit3) < abs(lit2):
                tmp = lit2
                lit2 = lit3
                lit3 = tmp
            assert abs(lit1) < abs(lit2)
            assert abs(lit1) < abs(lit3)
            assert abs(lit2) < abs(lit3)
            clause = str(lit1) + " " + str(lit2) + " " + str(lit3)
            if not clause in clauses:
                clauses[clause] = True
                problem += (clause + " 0\n")
                break
    return problem



def generate_3_sat_problems():
    # Generate problems for Fig.1 in Nusslein
    random.seed(1345)
    n_vars = np.arange(20,401,20)
    num_instances=20
    dir = "../exp/e1/problems"
    p_dir = Path(dir)
    p_dir.mkdir(parents=True, exist_ok=True)
    for vars in n_vars:
        for i in range(num_instances):
            p = generate_3sat(vars, ratio=4.2)
            with open(p_dir / ("p"+str(vars)+"-"+str(i)+".cnf"),"w") as f:
                f.write(p)

    # Generate problems for Table.2
    random.seed(13435)
    n_vars = np.arange(15,28,3)
    num_instances = 20
    dir = "../exp/e2/problems"
    p_dir = Path(dir)
    p_dir.mkdir(parents=True, exist_ok=True)
    for vars in n_vars:
        for i in range(num_instances):
            p = generate_3sat(vars, ratio=4.2)
            with open(p_dir / ("p"+str(vars)+"-"+str(i)+".cnf"),"w") as f:
                f.write(p)

    # Generate problems for Table.1 and Fig.2
    random.seed(178)
    n_vars=[5,10,12,20,50]
    num_instances = 20
    dir = "../exp/e3/problems"
    p_dir = Path(dir)
    p_dir.mkdir(parents=True, exist_ok=True)
    for vars in n_vars:
        for i in range(num_instances):
            p = generate_3sat(vars, ratio=4.2)
            with open(p_dir / ("p"+str(vars)+"-"+str(i)+".cnf"), "w") as f:
                f.write(p)


def parse_cnf_file(file_name):
    clauses=[]
    with open(file_name, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith('c') and not line.startswith('p') and line:
                clauses.append(np.array(list(map(int, line.split()[:-1]))))
            if line.startswith('p'):
                V = int(line.split()[2])
                
    return clauses, V



def get_non_zero_coupings(Q):
    count = 0
    for key, value in Q.items():
        if key[0] != key[1]:
            if value != 0:
                count += 1
    return count



def get_embedding(Q):
    couplings = []
    for key, value in Q.items():
        if key[0] != key[1]:
            if value != 0:
                couplings.append(key)
    try:
        embedding = find_embedding(couplings, DWaveSampler().edgelist, random_seed=10)
    except Exception as e:
        embedding = None
    
    return embedding



def solve_with_DWave(Q):    
    sampler = EmbeddingComposite(DWaveSampler(token='Your_Token'))
    response = sampler.sample_qubo(Q, num_reads=100, annealing_time=100, reduce_intersample_correlation=True)

    return response



def count_unsatisfied_clauses(assignment, clauses):
    count = 0
    for c in clauses:
        for l in c:
            if (l > 0 and assignment[np.abs(l)-1] == 1) or (l < 0 and assignment[np.abs(l)-1] == -1):
                count += 1
                break
                
    return len(clauses)-count



def compare_with_exact(g, vars, i, o):
    file_name = f"../exp/e3/Maxsatz_(exact)/p{vars}-{i}.txt"
    with open(file_name, "r") as archivo:
        lines = archivo.readlines()
        for l in reversed(lines):
            if l.split()[0] == 'o':
                optimum = int(l.split()[1])
                break

    return o-optimum