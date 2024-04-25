class choi:

    def __init__(self, file_path):
        clauses, V = utils.parse_cnf_file(file_path)
        # sort the clauses (i.e. all negative literals are at the back of the clause)
        self.clauses = clauses
        self.L = []
        for c in clauses:
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

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        answer = response.first.sample
        assignment = [0 for _ in range(self.V)]
        for i in range(len(self.L)):
            if answer[i] == 1:
                if self.L[i] < 0:
                    assignment[abs(self.L[i])-1] = 0
                else:
                    assignment[abs(self.L[i])-1] = 1
        
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']
        
        #Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write("QPU access time:"+str(qpu_access_time*10**(-6))+"s."+" ; QPU sampling time:"+str(qpu_sampling_time*10**(-6))+"s.\n")
            archivo.write("o "+utils.count_unsatisfied_clauses(assignment, self.clauses)+"\n")
            archivo.write("e "+str(response.first.energy)+"\n")
            archivo.write("v "+str(response.first.sample)+"\n")
            archivo.write(str(response.samples))
