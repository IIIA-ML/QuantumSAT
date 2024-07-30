import numpy as np

class Splitter:
    def Parse_cnf_file(self, file_name):
        clauses=[]
        with open(file_name, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if not line.startswith('c') and not line.startswith('p') and line:
                    clauses.append(np.array(list(map(int, line.split()[:-1]))))
                if line.startswith('p'):
                    V = int(line.split()[2])
                    
        return clauses, V
    
    def Split(self, clauses):
        pass
    
    
class Single_problem(Splitter):
    def Split(self, clauses):
        return [clauses]

class Two_subproblems(Splitter):
    def Split(self, clauses):
        return [clauses[:int(len(clauses)/2)], clauses[int(len(clauses)/2):]]

class Multiple_five(Splitter):
    def Split(self, clauses):
        clause_list=[]
        new_clauses=[]
        for i, c in enumerate(clauses):
            if i==len(clauses)-1 and i%21!=0:
                new_clauses.append(clause_list)
                break
            if i!=0 and i%21==0:
                new_clauses.append(clause_list)
                clause_list=[c]
            else:
                clause_list.append(c)
        return new_clauses

    #Separar en subgrups segons el nombre de literals en clausules (per exemple)
    #def k_SAT(self):
    #Si es genera un subgrup d'una sola clausula, aquest subgrup s'ha de definir com [[2,3]], per aixi quan es crida clauses es pot fer for c in clauses