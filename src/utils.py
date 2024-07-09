import random
import numpy as np
from pathlib import Path


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