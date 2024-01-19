#!/usr/bin/env python

import random

# Problem generation functions

# Adapted from Robert Daniel Brummayer 3satgen code

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
                   
if __name__ == "__main__":
    print(generate_3sat(10,ratio=4))         




