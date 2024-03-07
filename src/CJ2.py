import utils
import numpy as np


class CJ2:

    def __init__(self, file_path):
        self.clauses, self.V = utils.parse_cnf_file(file_path)
        self.Q = {}

    def add(self, x, y, value):
        x = np.abs(x) - 1
        y = np.abs(y) - 1
        if x > y:
            x,y = y,x
        if (x,y) in self.Q.keys():
            self.Q[(x,y)] += value
        else:
            self.Q[(x,y)] = value

    def fillQ(self):
        for i, c in enumerate(self.clauses):
            s = [1 if l>0 else -1 for l in c]
            var = [abs(l) for l in c]
            self.add(var[0], var[0], s[0]-s[0]*s[2])
            self.add(var[1], var[1], -2*s[1])
            self.add(var[2], var[2], s[2]-s[0]*s[2])
            self.add(self.V+i+1, self.V+i+1, -1+s[0]-s[1]+s[2])
            self.add(var[0], var[2], 2*s[0]*s[2])
            self.add(var[0], self.V+i+1, -2*s[0])
            self.add(var[1], self.V+i+1, 2*s[1])
            self.add(var[2], self.V+i+1, -2*s[2])

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