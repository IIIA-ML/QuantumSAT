import utils
import random
import Splitter
import Gadget
import Joiner
import Solver



class Variable_encoding:
    def __init__(self, encoding=None, variables=None, variable_to_literal=None, literal_to_variable=None):
        if variables is not None:
            self.variables = variables
        else:
            self.variables = set([v for v in encoding.values() if v is not None])
        if variable_to_literal is not None:
            self.variable_to_literal = variable_to_literal
        else:
            variable_to_literal = {}
            for variable in self.variables:
                variable_to_literal[variable] = []
                for k, v in encoding.items():
                    if variable==v:
                        variable_to_literal[v].append(k)
            self.variable_to_literal = variable_to_literal
        if literal_to_variable is not None:
            self.literal_to_variable = literal_to_variable
        else:
            self.literal_to_variable = encoding



class SAT_solver:
    def __init__(self, file_path=None, clauses=None, N=None, Splitter=None, Gadget=None, Joiner=None, Solver=None, token=None, qubit_level=False):
        self.file_path = file_path    
        if file_path is not None:
            self.clauses, self.N = Splitter.Parse_cnf_file(file_path)
            for clause in self.clauses:
                random.shuffle(clause)
        else:
            self.clauses = clauses,
            self.N = N
        self.Splitter = Splitter
        self.Gadget = Gadget
        self.Joiner = Joiner
        self.Solver = Solver
        self.token = token
        self.qubit_level = qubit_level
        
    def Solve(self):
        if self.Splitter is not None:
            self.subproblems = self.Splitter.Split(self.clauses)
            if self.Gadget is not None:
                Mappers = []
                for subproblem in self.subproblems:
                    Mappers.append(self.Gadget.FillQ(subproblem)) #[Q, encoding]
                if self.Joiner is not None:
                    self.Q, Q_encoding = self.Joiner.Join(Mappers)
                    self.Q_encoding = Variable_encoding(encoding = Q_encoding)
                    if self.Solver is not None:
                        self.response, self.embedding = self.Solver.Solve(self.Q, self.token, self.qubit_level)
                        embedding = self.response.info['embedding_context']['embedding']



