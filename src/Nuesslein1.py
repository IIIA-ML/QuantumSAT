import utils


class Nuesslein1:

    def __init__(self, file_path):
        clauses, V = utils.parse_cnf_file(file_path)
        # sort the clauses (i.e. all negative literals are at the back of the clause)
        self.clauses = [sorted(c, reverse=True) for c in clauses]
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
        for c in self.clauses:
            if x in c:
                n += 1
        return n

    def R2(self, x, y):
        n = 0
        for c in self.clauses:
            if x in c and y in c:
                n += 1
        return n

    # this function creates the QUBO-Matrix Q
    def fillQ(self):
        for i in range(2*self.V + len(self.clauses)):
            for j in range(2*self.V + len(self.clauses)):
                if i > j:
                    continue
                if i == j and j < 2*self.V:
                    self.add(i, j, -self.R1(self.L[i]))
                elif i == j and j >= 2*self.V:
                    self.add(i, j, 2)
                elif j < 2*self.V and j-i == 1 and i%2 == 0:
                    self.add(i, j, len(self.clauses)+1)
                elif i < 2*self.V and j < 2*self.V:
                    self.add(i, j, self.R2(self.L[i], self.L[j]))
                elif j >= 2*self.V and i < 2*self.V and self.L[i] in self.clauses[j-2*self.V]:
                    self.add(i, j, -1)

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        answer = response.first.sample
        assignment = [list(answer.values())[2*i] for i in range(self.V)]
        
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']
        
        #Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write("QPU access time:"+str(qpu_access_time*10**(-6))+"s."+" ; QPU sampling time:"+str(qpu_sampling_time*10**(-6))+"s.\n")
            archivo.write("o "+utils.count_unsatisfied_clauses(assignment, self.clauses)+"\n")
            archivo.write("e "+str(response.first.energy)+"\n")
            archivo.write("v "+str(response.first.sample)+"\n")
            archivo.write(str(response.samples))