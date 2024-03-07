import utils
import numpy as np


class Nuesslein2:

    def __init__(self, file_path):
        clauses, V = utils.parse_cnf_file(file_path)
        # sort the clauses (i.e. all negative literals are at the back of the clause)
        self.clauses = [sorted(c, reverse=True) for c in clauses]
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
        for i, c in enumerate(self.clauses):
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

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        answer = response.first.sample
        assignment = [list(answer.values())[i] for i in range(self.V)]
        
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']
        
        #Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write("QPU access time:"+str(qpu_access_time*10**(-6))+"s."+" ; QPU sampling time:"+str(qpu_sampling_time*10**(-6))+"s.\n")
            archivo.write("o "+utils.count_unsatisfied_clauses(assignment, self.clauses)+"\n")
            archivo.write("e "+str(response.first.energy)+"\n")
            archivo.write("v "+str(response.first.sample)+"\n")
            archivo.write(str(response.samples))