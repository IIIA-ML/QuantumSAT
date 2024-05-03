import utils
import numpy as np


class CJ1_bian:

    def __init__(self, file_path=None, clauses=None, V=None):
        if file_path is not None:
            self.clauses, self.V = utils.parse_cnf_file(file_path)
        else:
            self.clauses = clauses
            self.V = V
        self.b = int(3*len(self.clauses))
        self.Q = {}
        self.variables = {}
        self.sign = []

    def add(self, a, b, x, y, value):
        if a not in self.variables.keys():
            self.variables[a] = x
        if b not in self.variables.keys():
            self.variables[b] = y
        if a > b:
            a,b = b,a
        if (a,b) in self.Q.keys():
            self.Q[(a,b)] += value
        else:
            self.Q[(a,b)] = value

    def fillQ(self):
        k=0
        for i, c in enumerate(self.clauses):
            s = [1 if l>0 else -1 for l in c]
            var = [abs(l) for l in c]
            self.add(k, k, var[0], var[0], 3*s[0]-s[0]*s[1])
            self.add(k+1, k+1, var[1], var[1], 3*s[1]-s[0]*s[1])
            self.add(k+2, k+2, var[2], var[2], -2*s[2])
            self.add(self.b+i, self.b+i, self.V+i+1, self.V+i+1, -3+2*s[0]+2*s[1]-s[2])
            self.add(k, k+1, var[0], var[1], 2*s[0]*s[1])
            self.add(k, self.b+i, var[0], self.V+i+1, -4*s[0])
            self.add(k+1, self.b+i, var[1], self.V+i+1, -4*s[1])
            self.add(k+2, self.b+i, var[2], self.V+i+1, 2*s[2])
            k+=3
            self.sign = self.sign + s

        for v in range(1,self.V+1):
            repeated_variables = sorted([i for i, variable in self.variables.items() if variable==v])
            '''for j in range(len(repeated_variables)):
                for k in range(j+1, len(repeated_variables)):
                    self.add(repeated_variables[j], repeated_variables[j], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[k]])
                    self.add(repeated_variables[k], repeated_variables[k], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[k]])
                    self.add(repeated_variables[j], repeated_variables[k], v, v, -4*self.sign[repeated_variables[j]]*self.sign[repeated_variables[k]])'''
            for j in range(len(repeated_variables)-1):
                self.add(repeated_variables[j], repeated_variables[j], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                self.add(repeated_variables[j+1], repeated_variables[j+1], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                self.add(repeated_variables[j], repeated_variables[j+1], v, v, -4*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']
        
        #Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write("QPU access time:"+str(qpu_access_time*10**(-6))+"s."+" ; QPU sampling time:"+str(qpu_sampling_time*10**(-6))+"s.\n")
            archivo.write("o "+utils.count_unsatisfied_clauses(response.first.sample, self.clauses)+"\n")
            archivo.write("e "+str(response.first.energy)+"\n")
            archivo.write("v "+str(response.first.sample)+"\n")
            archivo.write(str(response.samples))
