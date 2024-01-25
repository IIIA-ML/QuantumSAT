#from dwave_qbsolv import QBSolv
import neal
import dimod
import numpy as np
import time


# this function downloads a random k-SAT formula
def create_formula(num_vars, num_clauses, k):

    formula = []
    while len(formula) < num_clauses:
        vars = np.random.choice(range(1,num_vars+1), size=k, replace=False)
        signs = np.random.choice([-1,+1], size=k, replace=True)
        formula.append(vars * signs)

    return formula


# this function solves a given QUBO-Matrix Q with Qbsolv
"""def solve_with_qbsolv(Q):
    response = QBSolv().sample_qubo(Q, num_repeats=1000)
    return response.samples()[0]"""

#NEWWWWWW
def solve_with_SA(Q):
    model = dimod.BinaryQuadraticModel.from_qubo(Q)
    sampler = neal.SimulatedAnnealingSampler()
    start_time = time.time()
    response = sampler.sample(model, num_reads=1000) #neal.SimulatedAnnealingSampler().sample_qubo(Q,num_repeats=1000)
    end_time = time.time()
    return response, end_time-start_time #response.samples()[0]

# this function calculates the value of a solution for a given QUBO-Matrix Q
def getValue(Q, solution):
    ones = [x for x in solution.keys() if solution[x] == 1]
    value = 0
    for x in ones:
        for y in ones:
            if (x,y) in Q.keys():
                value += Q[(x,y)]
    return value


# this function prints the first n row/columns of a QUBO-Matrix Q
def printQUBO(Q, n):
    for row in range(n):
        for column in range(n):
            if row > column:
                print("      ", end = '')
                continue
            printing = ""
            if (row,column) in Q.keys() and Q[(row,column)] != 0:
                printing = str(Q[(row,column)])
            printing += "_____"
            printing = printing[:5]
            printing += " "
            print(printing, end = '')
        print("")


# this function checks, whether a given assignment satisfies a given SAT-formula
def check_solution(formula, assignment):
    n = 0
    for c in formula:
        for l in c:
            if l < 0 and assignment[abs(l)-1] == 0:
                n += 1
                break
            elif l > 0 and assignment[abs(l)-1] == 1:
                n += 1
                break
    return n


def get_n_couplings(Q):
    n = 0
    for k in Q.keys():
        if Q[k] != 0:
            n += 1
    return n

class nuesslein2:

    def __init__(self, formula, V):
        # sort the formula (i.e. all negative literals are at the back of the clause)
        self.formula = [sorted(c, reverse=True) for c in formula]
        self.V = V
        self.Q = {}

    # new values are added to the QUBO-Matrix Q via this monitor
    def add(self, x, y, value):
        x = np.abs(x) - 1
        y = np.abs(y) - 1
        if x > y:
            x,y = y,x
        if (x,y) in self.Q.keys():
            self.Q[(x,y)] += value
        else:
            self.Q[(x,y)] = value

    # this function creates the QUBO-Matrix Q
    def fillQ(self):
        for i, c in enumerate(self.formula):
            if list(np.sign(c)) == [1, 1, 1]:
                self.add(c[0], c[1], 2)
                self.add(c[0], self.V + i + 1, -2)
                self.add(c[1], self.V + i + 1, -2)
                self.add(c[2], c[2], -1)
                self.add(c[2], self.V + i + 1, 1)
                self.add(self.V + i + 1, self.V + i + 1, 1)
            elif list(np.sign(c)) == [1, 1, -1]:
                self.add(c[0], c[1], 2)
                self.add(c[0], self.V + i + 1, -2)
                self.add(c[1], self.V + i + 1, -2)
                self.add(c[2], c[2], 1)
                self.add(c[2], self.V + i + 1, -1)
                self.add(self.V + i + 1, self.V + i + 1, 2)
            elif list(np.sign(c)) == [1, -1, -1]:
                self.add(c[0], c[0], 2)
                self.add(c[0], c[1], -2)
                self.add(c[0], self.V + i + 1, -2)
                self.add(c[1], self.V + i + 1, 2)
                self.add(c[2], c[2], 1)
                self.add(c[2], self.V + i + 1, -1)
            else:
                self.add(c[0], c[0], -1)
                self.add(c[0], c[1], 1)
                self.add(c[0], c[2], 1)
                self.add(c[0], self.V + i + 1, 1)
                self.add(c[1], c[1], -1)
                self.add(c[1], c[2], 1)
                self.add(c[1], self.V + i + 1, 1)
                self.add(c[2], c[2], -1)
                self.add(c[2], self.V + i + 1, 1)
                self.add(self.V + i + 1, self.V + i + 1, -1)

    # this function starts creating Q, solving it and interpreting the solution
    # (e.g. deciding whether the formula is satisfiable or not)
    def solve(self):
        self.fillQ()
        response, time = solve_with_SA(self.Q)
        answer = response.first.sample
        assignment = [list(answer.values())[i] for i in range(self.V)]
        non_zero_couplings = get_n_couplings(self.Q)
        return time, check_solution(self.formula, assignment), response, non_zero_couplings


class choi:

    def __init__(self, formula, V):
        # sort the formula (i.e. all negative literals are at the back of the clause)
        self.L = []
        self.formula = formula
        for c in formula:
            self.L.extend(c)
        self.V = V
        self.Q = {}

    # new values are added to the QUBO-Matrix Q via this monitor
    def add(self, x, y, value):
        if x > y:
            x,y = y,x
        if (x,y) in self.Q.keys():
            self.Q[(x,y)] += value
        else:
            self.Q[(x,y)] = value

    # this function creates the QUBO-Matrix Q
    # Explanations can be found in the paper
    def fillQ(self):
        for i in range(len(self.L)):
            for j in range(len(self.L)):
                if i > j:
                    continue
                if i == j:
                    self.add(i, j, -1)
                elif j - i <= 2 and j//3 == i//3:
                    self.add(i, j, 3)
                elif abs(self.L[i]) == abs(self.L[j]) and self.L[i] != self.L[j]:
                    self.add(i, j, 3)

    # this function starts creating Q, solving it and interpreting the solution
    # (e.g. deciding whether the formula is satisfiable or not)
    def solve(self):
        self.fillQ()
        response, time = solve_with_SA(self.Q)
        answer = response.first.sample
        assignment = [0 for _ in range(self.V)]
        for i in range(len(self.L)):
            if list(answer.values())[i] == 1:
                if self.L[i] < 0:
                    assignment[abs(self.L[i])-1] = 0
                else:
                    assignment[abs(self.L[i])-1] = 1
        non_zero_couplings = get_n_couplings(self.Q)
        return time, check_solution(self.formula, assignment), response, non_zero_couplings


class chancellor:

    def __init__(self, formula, V):
        # sort the formula (i.e. all negative literals are at the back of the clause)
        self.formula = [sorted(c, reverse=True) for c in formula]
        self.V = V
        self.Q = {}

    # new values are added to the QUBO-Matrix Q via this monitor
    def add(self, x, y, value):
        x = np.abs(x) - 1
        y = np.abs(y) - 1
        if x > y:
            x,y = y,x
        if (x,y) in self.Q.keys():
            self.Q[(x,y)] += value
        else:
            self.Q[(x,y)] = value

    # this function creates the QUBO-Matrix Q
    def fillQ(self):
        for i, c in enumerate(self.formula):
            if list(np.sign(c)) == [1, 1, 1]:
                self.add(c[0], c[0], -48)
                self.add(c[0], c[1], 24)
                self.add(c[0], c[2], 24)
                self.add(c[0], self.V + i + 1, 40)
                self.add(c[1], c[1], -48)
                self.add(c[1], c[2], 24)
                self.add(c[1], self.V + i + 1, 40)
                self.add(c[2], c[2], -48)
                self.add(c[2], self.V + i + 1, 40)
                self.add(self.V + i + 1, self.V + i + 1, -64)
            elif list(np.sign(c)) == [1, 1, -1]:
                self.add(c[0], c[0], -40)
                self.add(c[0], c[1], 24)
                self.add(c[0], c[2], 16)
                self.add(c[0], self.V + i + 1, 40)
                self.add(c[1], c[1], -40)
                self.add(c[1], c[2], 16)
                self.add(c[1], self.V + i + 1, 40)
                self.add(c[2], c[2], -32)
                self.add(c[2], self.V + i + 1, 40)
                self.add(self.V + i + 1, self.V + i + 1, -56)
            elif list(np.sign(c)) == [1, -1, -1]:
                self.add(c[0], c[0], -40)
                self.add(c[0], c[1], 16)
                self.add(c[0], c[2], 16)
                self.add(c[0], self.V + i + 1, 40)
                self.add(c[1], c[1], -40)
                self.add(c[1], c[2], 24)
                self.add(c[1], self.V + i + 1, 40)
                self.add(c[2], c[2], -40)
                self.add(c[2], self.V + i + 1, 40)
                self.add(self.V + i + 1, self.V + i + 1, -64)
            else:
                self.add(c[0], c[0], -40)
                self.add(c[0], c[1], 24)
                self.add(c[0], c[2], 24)
                self.add(c[0], self.V + i + 1, 40)
                self.add(c[1], c[1], -40)
                self.add(c[1], c[2], 24)
                self.add(c[1], self.V + i + 1, 40)
                self.add(c[2], c[2], -40)
                self.add(c[2], self.V + i + 1, 40)
                self.add(self.V + i + 1, self.V + i + 1, -56)

    # this function starts creating Q, solving it and interpreting the solution
    # (e.g. deciding whether the formula is satisfiable or not)
    def solve(self):
        self.fillQ()
        #answer = utils.solve_with_qbsolv(self.Q)
        response, time = solve_with_SA(self.Q)
        answer = response.first.sample
        assignment = [list(answer.values())[i] for i in range(self.V)]
        non_zero_couplings = get_n_couplings(self.Q)
        return time, check_solution(self.formula, assignment), response, non_zero_couplings


class nuesslein1:

    def __init__(self, formula, V):
        # sort the formula (i.e. all negative literals are at the back of the clause)
        self.formula = [sorted(c, reverse=True) for c in formula]
        self.L = []
        for i in range(V):
            self.L.append(i+1)
            self.L.append(-(i+1))
        self.V = V
        self.Q = {}

    # new values are added to the QUBO-Matrix Q via this monitor
    def add(self, x, y, value):
        if x > y:
            x,y = y,x
        if (x,y) in self.Q.keys():
            self.Q[(x,y)] += value
        else:
            self.Q[(x,y)] = value

    def R1(self, x):
        n = 0
        for c in self.formula:
            if x in c:
                n += 1
        return n

    def R2(self, x, y):
        n = 0
        for c in self.formula:
            if x in c and y in c:
                n += 1
        return n

    # this function creates the QUBO-Matrix Q
    def fillQ(self):
        for i in range(2*self.V + len(self.formula)):
            for j in range(2*self.V + len(self.formula)):
                if i > j:
                    continue
                if i == j and j < 2*self.V:
                    self.add(i, j, -self.R1(self.L[i]))
                elif i == j and j >= 2*self.V:
                    self.add(i, j, 2)
                elif j < 2*self.V and j-i == 1 and i%2 == 0:
                    self.add(i, j, len(self.formula)+1)
                elif i < 2*self.V and j < 2*self.V:
                    self.add(i, j, self.R2(self.L[i], self.L[j]))
                elif j >= 2*self.V and i < 2*self.V and self.L[i] in self.formula[j-2*self.V]:
                    self.add(i, j, -1)

    # this function starts creating Q, solving it and interpreting the solution
    # (e.g. deciding whether the formula is satisfiable or not)
    def solve(self):
        self.fillQ()
        response, time = solve_with_SA(self.Q)
        answer = response.first.sample
        assignment = [list(answer.values())[2*i] for i in range(self.V)]
        non_zero_couplings = get_n_couplings(self.Q)
        return time, check_solution(self.formula, assignment), response, non_zero_couplings