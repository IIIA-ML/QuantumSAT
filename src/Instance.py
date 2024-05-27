import utils
import numpy as np
import dimod
from minorminer import find_embedding
from dwave.system import DWaveSampler, FixedEmbeddingComposite
import dwave.embedding


class Instance:
    def __init__(self, file_path=None, clauses=None, N=None):
        self.file_path = file_path
        if file_path is not None:
            self.clauses, self.N = utils.parse_cnf_file(file_path)
        else:
            self.clauses = clauses,
            self.N = N
        self.Nuesslein1 = Nuesslein1(self.clauses, self.N)
        self.Nuesslein2 = Nuesslein2(self.clauses, self.N)
        self.CJ1 = CJ1(self.clauses, self.N)
        self.CJ2 = CJ2(self.clauses, self.N)
        self.CJ1_bian = CJ1_bian(self.clauses, self.N)
        self.CJ2_bian = CJ2_bian(self.clauses, self.N)


class Nuesslein1(Instance):

    def __init__(self, clauses, N):
        self.clauses = clauses
        self.N = N
        # sort the clauses (i.e. all negative literals are at the back of the clause)
        self.clauses = [sorted(c, reverse=True) for c in clauses]
        self.L = []
        for i in range(N):
            self.L.append(i + 1)
            self.L.append(-(i + 1))
        self.Q = {}
        self.bqm = None
        self.bqm_embedded = {}
        self.scaled_bqm_embedded = {}
        self.scale_factor = None
        self.response_dwave_from_scaled_bqm_embedded = None
        self.all_solution_space = {}

    # new values are added to the QUBO-Matrix Q via this monitor
    def add(self, x, y, value):
        if x > y:
            x, y = y, x
        if (x, y) in self.Q.keys():
            self.Q[(x, y)] += value
        else:
            self.Q[(x, y)] = value

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
        if self.Q:
            raise ValueError("Q is already filled")
        for i in range(2 * self.N + len(self.clauses)):
            for j in range(2 * self.N + len(self.clauses)):
                if i > j:
                    continue
                if i == j and j < 2 * self.N:
                    self.add(i, j, -self.R1(self.L[i]))
                elif i == j and j >= 2 * self.N:
                    self.add(i, j, 2)
                elif j < 2 * self.N and j - i == 1 and i % 2 == 0:
                    self.add(i, j, len(self.clauses) + 1)
                elif i < 2 * self.N and j < 2 * self.N:
                    self.add(i, j, self.R2(self.L[i], self.L[j]))
                elif j >= 2 * self.N and i < 2 * self.N and self.L[i] in self.clauses[j - 2 * self.N]:
                    self.add(i, j, -1)

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        answer = response.first.sample
        assignment = [list(answer.values())[2 * i] for i in range(self.N)]

        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']

        # Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write("QPU access time:" + str(qpu_access_time * 10 ** (-6)) + "s." + " ; QPU sampling time:" + str(
                qpu_sampling_time * 10 ** (-6)) + "s.\n")
            archivo.write("o " + utils.count_unsatisfied_clauses(assignment, self.clauses) + "\n")
            archivo.write("e " + str(response.first.energy) + "\n")
            archivo.write("v " + str(response.first.sample) + "\n")
            archivo.write(str(response.samples))

    def compute_bqm(self):
        if not self.Q:
            self.fillQ()
        bqm = dimod.BinaryQuadraticModel.from_qubo(self.Q, offset=0).change_vartype("SPIN", True)
        self.bqm = bqm
        return bqm

    def compute_bqm_embedded(self, token):
        couplings = []
        if not self.Q:
            self.fillQ()
        for k in self.Q.keys():
            if k[0] != k[1]:
                couplings.append(k)
        embedding = find_embedding(couplings, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
        if not self.bqm:
            self.compute_bqm()
        bqm_embedded = dwave.embedding.embed_bqm(source_bqm=self.bqm, embedding=embedding,
                                                              target_adjacency=DWaveSampler(token=token).adjacency)
        self.bqm_embedded['bqm_embedded'] = bqm_embedded
        self.bqm_embedded['embedding'] = embedding
        return {'bqm_embedded': bqm_embedded, 'embedding': embedding}

    def compute_scaled_bqm_embedded(self, token):
        if not self.bqm_embedded:
            self.compute_bqm_embedded(token)
        if not self.scaled_bqm_embedded:
            scale_factor = utils.compute_scale_factor(self.bqm_embedded['bqm_embedded'])
            self.scale_factor = scale_factor
            temp_bqm_embedded = self.bqm_embedded['bqm_embedded'].copy()
            if scale_factor > 1:
                temp_bqm_embedded.scale(1 / scale_factor)
            self.scaled_bqm_embedded['scaled_bqm_embedded'] = temp_bqm_embedded
            self.scaled_bqm_embedded['scale_factor'] = scale_factor
            self.scale_factor = scale_factor
        return {'scaled_bqm_embedded': self.scaled_bqm_embedded['scaled_bqm_embedded'], 'scale_factor': self.scale_factor}

    def solve_dwave_from_scaled_bqm_embedded(self, token, num_reads=100, annealing_time=100):
        if not self.scaled_bqm_embedded:
            self.compute_scaled_bqm_embedded(token)
        if not self.response_dwave_from_scaled_bqm_embedded:
            embedding_refactored = utils.refactor_embedding_1_to_1(self.scaled_bqm_embedded['scaled_bqm_embedded'])
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=embedding_refactored)
            response = sampler.sample(self.scaled_bqm_embedded['scaled_bqm_embedded'], num_reads=num_reads, annealing_time=annealing_time,
                                      return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            response_aux = response.to_serializable()
            self.response_dwave_from_scaled_bqm_embedded = dimod.SampleSet.from_serializable(response_aux)
        return self.response_dwave_from_scaled_bqm_embedded


class Nuesslein2(Instance):

    def __init__(self, clauses, N):
        self.clauses = clauses
        self.N = N
        # sort the clauses (i.e. all negative literals are at the back of the clause)
        self.clauses = [sorted(c, reverse=True) for c in clauses]
        self.N = N
        self.Q = {}
        self.bqm = None
        self.bqm_embedded = {}
        self.scaled_bqm_embedded = {}
        self.scale_factor = None
        self.response_dwave_from_scaled_bqm_embedded = None
        self.all_solution_space = {}

    # new values are added to the QUBO-Matrix Q via this monitor
    def add(self, x, y, value):
        x = np.abs(x) - 1
        y = np.abs(y) - 1
        if x > y:
            x, y = y, x
        if (x, y) in self.Q.keys():
            self.Q[(x, y)] += value
        else:
            self.Q[(x, y)] = value

    # this function creates the QUBO-Matrix Q
    def fillQ(self):
        if self.Q:
            raise ValueError("Q is already filled")
        for i, c in enumerate(self.clauses):
            if list(np.sign(c)) == [1, 1, 1]:
                self.add(c[0], c[1], 2)
                self.add(c[0], self.N + i + 1, -2)
                self.add(c[1], self.N + i + 1, -2)
                self.add(c[2], c[2], -1)
                self.add(c[2], self.N + i + 1, 1)
                self.add(self.N + i + 1, self.N + i + 1, 1)
            elif list(np.sign(c)) == [1, 1, -1]:
                self.add(c[0], c[1], 2)
                self.add(c[0], self.N + i + 1, -2)
                self.add(c[1], self.N + i + 1, -2)
                self.add(c[2], c[2], 1)
                self.add(c[2], self.N + i + 1, -1)
                self.add(self.N + i + 1, self.N + i + 1, 2)
            elif list(np.sign(c)) == [1, -1, -1]:
                self.add(c[0], c[0], 2)
                self.add(c[0], c[1], -2)
                self.add(c[0], self.N + i + 1, -2)
                self.add(c[1], self.N + i + 1, 2)
                self.add(c[2], c[2], 1)
                self.add(c[2], self.N + i + 1, -1)
            else:
                self.add(c[0], c[0], -1)
                self.add(c[0], c[1], 1)
                self.add(c[0], c[2], 1)
                self.add(c[0], self.N + i + 1, 1)
                self.add(c[1], c[1], -1)
                self.add(c[1], c[2], 1)
                self.add(c[1], self.N + i + 1, 1)
                self.add(c[2], c[2], -1)
                self.add(c[2], self.N + i + 1, 1)
                self.add(self.N + i + 1, self.N + i + 1, -1)

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        answer = response.first.sample
        assignment = [list(answer.values())[i] for i in range(self.N)]

        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']

        # Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write("QPU access time:" + str(qpu_access_time * 10 ** (-6)) + "s." + " ; QPU sampling time:" + str(
                qpu_sampling_time * 10 ** (-6)) + "s.\n")
            archivo.write("o " + utils.count_unsatisfied_clauses(assignment, self.clauses) + "\n")
            archivo.write("e " + str(response.first.energy) + "\n")
            archivo.write("v " + str(response.first.sample) + "\n")
            archivo.write(str(response.samples))

    def compute_bqm(self):
        if not self.Q:
            self.fillQ()
        bqm = dimod.BinaryQuadraticModel.from_qubo(self.Q, offset=0).change_vartype("SPIN", True)
        self.bqm = bqm
        return bqm

    def compute_bqm_embedded(self, token):
        couplings = []
        if not self.Q:
            self.fillQ()
        for k in self.Q.keys():
            if k[0] != k[1]:
                couplings.append(k)
        embedding = find_embedding(couplings, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
        if not self.bqm:
            self.compute_bqm()
        bqm_embedded = dwave.embedding.embed_bqm(source_bqm=self.bqm, embedding=embedding,
                                                              target_adjacency=DWaveSampler(token=token).adjacency)
        self.bqm_embedded['bqm_embedded'] = bqm_embedded
        self.bqm_embedded['embedding'] = embedding
        return {'bqm_embedded': bqm_embedded, 'embedding': embedding}

    def compute_scaled_bqm_embedded(self, token):
        if not self.bqm_embedded:
            self.compute_bqm_embedded(token)
        if not self.scaled_bqm_embedded:
            scale_factor = utils.compute_scale_factor(self.bqm_embedded['bqm_embedded'])
            self.scale_factor = scale_factor
            temp_bqm_embedded = self.bqm_embedded['bqm_embedded'].copy()
            if scale_factor > 1:
                temp_bqm_embedded.scale(1 / scale_factor)
            self.scaled_bqm_embedded['scaled_bqm_embedded'] = temp_bqm_embedded
            self.scaled_bqm_embedded['scale_factor'] = scale_factor
            self.scale_factor = scale_factor
        return {'scaled_bqm_embedded': self.scaled_bqm_embedded['scaled_bqm_embedded'], 'scale_factor': self.scale_factor}

    def solve_dwave_from_scaled_bqm_embedded(self, token, num_reads=100, annealing_time=100):
        if not self.scaled_bqm_embedded:
            self.compute_scaled_bqm_embedded(token)
        if not self.response_dwave_from_scaled_bqm_embedded:
            embedding_refactored = utils.refactor_embedding_1_to_1(self.scaled_bqm_embedded['scaled_bqm_embedded'])
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=embedding_refactored)
            response = sampler.sample(self.scaled_bqm_embedded['scaled_bqm_embedded'], num_reads=num_reads, annealing_time=annealing_time,
                                      return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            response_aux = response.to_serializable()
            self.response_dwave_from_scaled_bqm_embedded = dimod.SampleSet.from_serializable(response_aux)
        return self.response_dwave_from_scaled_bqm_embedded


class CJ1(Instance):
    def __init__(self, clauses, N):
        self.clauses = clauses
        self.N = N
        self.Q = {}
        self.bqm = None
        self.bqm_embedded = {}
        self.scaled_bqm_embedded = {}
        self.scale_factor = None
        self.response_dwave_from_scaled_bqm_embedded = None
        self.all_solution_space = {}

    def add(self, x, y, value):
        x = np.abs(x) - 1
        y = np.abs(y) - 1
        if x > y:
            x, y = y, x
        if (x, y) in self.Q.keys():
            self.Q[(x, y)] += value
        else:
            self.Q[(x, y)] = value

    def fillQ(self):
        if self.Q:
            raise ValueError("Q is already filled")
        for i, c in enumerate(self.clauses):
            s = [1 if l > 0 else -1 for l in c]
            var = [abs(l) for l in c]
            self.add(var[0], var[0], 3 * s[0] - s[0] * s[1])
            self.add(var[1], var[1], 3 * s[1] - s[0] * s[1])
            self.add(var[2], var[2], -2 * s[2])
            self.add(self.N + i + 1, self.N + i + 1, -3 + 2 * s[0] + 2 * s[1] - s[2])
            self.add(var[0], var[1], 2 * s[0] * s[1])
            self.add(var[0], self.N + i + 1, -4 * s[0])
            self.add(var[1], self.N + i + 1, -4 * s[1])
            self.add(var[2], self.N + i + 1, 2 * s[2])


    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']

        # Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write("QPU access time:" + str(qpu_access_time * 10 ** (-6)) + "s." + " ; QPU sampling time:" + str(
                qpu_sampling_time * 10 ** (-6)) + "s.\n")
            archivo.write("o " + utils.count_unsatisfied_clauses(response.first.sample, self.clauses) + "\n")
            archivo.write("e " + str(response.first.energy) + "\n")
            archivo.write("v " + str(response.first.sample) + "\n")
            archivo.write(str(response.samples))

    def compute_bqm(self):
        if not self.Q:
            self.fillQ()
        bqm = dimod.BinaryQuadraticModel.from_qubo(self.Q, offset=0).change_vartype("SPIN", True)
        self.bqm = bqm
        return bqm

    def compute_bqm_embedded(self, token):
        couplings = []
        if not self.Q:
            self.fillQ()
        for k in self.Q.keys():
            if k[0] != k[1]:
                couplings.append(k)
        embedding = find_embedding(couplings, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
        if not self.bqm:
            self.compute_bqm()
        bqm_embedded = dwave.embedding.embed_bqm(source_bqm=self.bqm, embedding=embedding,
                                                              target_adjacency=DWaveSampler(token=token).adjacency)
        self.bqm_embedded['bqm_embedded'] = bqm_embedded
        self.bqm_embedded['embedding'] = embedding
        return {'bqm_embedded': bqm_embedded, 'embedding': embedding}

    def compute_scaled_bqm_embedded(self, token):
        if not self.bqm_embedded:
            self.compute_bqm_embedded(token)
        if not self.scaled_bqm_embedded:
            scale_factor = utils.compute_scale_factor(self.bqm_embedded['bqm_embedded'])
            self.scale_factor = scale_factor
            temp_bqm_embedded = self.bqm_embedded['bqm_embedded'].copy()
            if scale_factor > 1:
                temp_bqm_embedded.scale(1 / scale_factor)
            self.scaled_bqm_embedded['scaled_bqm_embedded'] = temp_bqm_embedded
            self.scaled_bqm_embedded['scale_factor'] = scale_factor
            self.scale_factor = scale_factor
        return {'scaled_bqm_embedded': self.scaled_bqm_embedded['scaled_bqm_embedded'], 'scale_factor': self.scale_factor}

    def solve_dwave_from_scaled_bqm_embedded(self, token, num_reads=100, annealing_time=100):
        if not self.scaled_bqm_embedded:
            self.compute_scaled_bqm_embedded(token)
        if not self.response_dwave_from_scaled_bqm_embedded:
            embedding_refactored = utils.refactor_embedding_1_to_1(self.scaled_bqm_embedded['scaled_bqm_embedded'])
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=embedding_refactored)
            response = sampler.sample(self.scaled_bqm_embedded['scaled_bqm_embedded'], num_reads=num_reads, annealing_time=annealing_time,
                                      return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            response_aux = response.to_serializable()
            self.response_dwave_from_scaled_bqm_embedded = dimod.SampleSet.from_serializable(response_aux)
        return self.response_dwave_from_scaled_bqm_embedded


class CJ2(Instance):
    def __init__(self, clauses, N):
        self.clauses = clauses
        self.N = N
        self.Q = {}
        self.bqm = None
        self.bqm_embedded = {}
        self.scaled_bqm_embedded = {}
        self.scale_factor = None
        self.response_dwave_from_scaled_bqm_embedded = None
        self.all_solution_space = {}

    def add(self, x, y, value):
        x = np.abs(x) - 1
        y = np.abs(y) - 1
        if x > y:
            x, y = y, x
        if (x, y) in self.Q.keys():
            self.Q[(x, y)] += value
        else:
            self.Q[(x, y)] = value

    def fillQ(self):
        if self.Q:
            raise ValueError("Q is already filled")
        for i, c in enumerate(self.clauses):
            s = [1 if l > 0 else -1 for l in c]
            var = [abs(l) for l in c]
            self.add(var[0], var[0], s[0] - s[0] * s[2])
            self.add(var[1], var[1], -2 * s[1])
            self.add(var[2], var[2], s[2] - s[0] * s[2])
            self.add(self.N + i + 1, self.N + i + 1, -1 + s[0] - s[1] + s[2])
            self.add(var[0], var[2], 2 * s[0] * s[2])
            self.add(var[0], self.N + i + 1, -2 * s[0])
            self.add(var[1], self.N + i + 1, 2 * s[1])
            self.add(var[2], self.N + i + 1, -2 * s[2])

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']

        # Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write(
                "QPU access time:" + str(qpu_access_time * 10 ** (-6)) + "s." + " ; QPU sampling time:" + str(
                    qpu_sampling_time * 10 ** (-6)) + "s.\n")
            archivo.write("o " + utils.count_unsatisfied_clauses(response.first.sample, self.clauses) + "\n")
            archivo.write("e " + str(response.first.energy) + "\n")
            archivo.write("v " + str(response.first.sample) + "\n")
            archivo.write(str(response.samples))

    def compute_bqm(self):
        if not self.Q:
            self.fillQ()
        bqm = dimod.BinaryQuadraticModel.from_qubo(self.Q, offset=0).change_vartype("SPIN", True)
        self.bqm = bqm
        return bqm

    def compute_bqm_embedded(self, token):
        couplings = []
        if not self.Q:
            self.fillQ()
        for k in self.Q.keys():
            if k[0] != k[1]:
                couplings.append(k)
        embedding = find_embedding(couplings, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
        if not self.bqm:
            self.compute_bqm()
        bqm_embedded = dwave.embedding.embed_bqm(source_bqm=self.bqm, embedding=embedding,
                                                              target_adjacency=DWaveSampler(token=token).adjacency)
        self.bqm_embedded['bqm_embedded'] = bqm_embedded
        self.bqm_embedded['embedding'] = embedding
        return {'bqm_embedded': bqm_embedded, 'embedding': embedding}

    def compute_scaled_bqm_embedded(self, token):
        if not self.bqm_embedded:
            self.compute_bqm_embedded(token)
        if not self.scaled_bqm_embedded:
            scale_factor = utils.compute_scale_factor(self.bqm_embedded['bqm_embedded'])
            self.scale_factor = scale_factor
            temp_bqm_embedded = self.bqm_embedded['bqm_embedded'].copy()
            if scale_factor > 1:
                temp_bqm_embedded.scale(1 / scale_factor)
            self.scaled_bqm_embedded['scaled_bqm_embedded'] = temp_bqm_embedded
            self.scaled_bqm_embedded['scale_factor'] = scale_factor
            self.scale_factor = scale_factor
        return {'scaled_bqm_embedded': self.scaled_bqm_embedded['scaled_bqm_embedded'], 'scale_factor': self.scale_factor}

    def solve_dwave_from_scaled_bqm_embedded(self, token, num_reads=100, annealing_time=100):
        if not self.scaled_bqm_embedded:
            self.compute_scaled_bqm_embedded(token)
        if not self.response_dwave_from_scaled_bqm_embedded:
            embedding_refactored = utils.refactor_embedding_1_to_1(self.scaled_bqm_embedded['scaled_bqm_embedded'])
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=embedding_refactored)
            response = sampler.sample(self.scaled_bqm_embedded['scaled_bqm_embedded'], num_reads=num_reads, annealing_time=annealing_time,
                                      return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            response_aux = response.to_serializable()
            self.response_dwave_from_scaled_bqm_embedded = dimod.SampleSet.from_serializable(response_aux)
        return self.response_dwave_from_scaled_bqm_embedded


class CJ1_bian(Instance):
    def __init__(self, clauses, N):
        self.clauses = clauses
        self.N = N
        self.b = int(3*len(self.clauses))
        self.Q = {}
        self.bian_embedding = {}
        self.sign = []
        self.bqm = None
        self.bqm_embedded = {}
        self.scaled_bqm_embedded = {}
        self.scale_factor = None
        self.response_dwave_from_scaled_bqm_embedded = None
        self.all_solution_space = {}

    def add(self, a, b, x, y, value):
        if a not in self.bian_embedding.keys():
            self.bian_embedding[a] = x
        if b not in self.bian_embedding.keys():
            self.bian_embedding[b] = y
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
            self.add(self.b+i, self.b+i, self.N+i+1, self.N+i+1, -3+2*s[0]+2*s[1]-s[2])
            self.add(k, k+1, var[0], var[1], 2*s[0]*s[1])
            self.add(k, self.b+i, var[0], self.N+i+1, -4*s[0])
            self.add(k+1, self.b+i, var[1], self.N+i+1, -4*s[1])
            self.add(k+2, self.b+i, var[2], self.N+i+1, 2*s[2])
            k+=3
            self.sign = self.sign + s

        for v in range(1,self.N+1):
            repeated_variables = sorted([i for i, variable in self.bian_embedding.items() if variable==v])
            for j in range(len(repeated_variables)-1):
                self.add(repeated_variables[j], repeated_variables[j], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                self.add(repeated_variables[j+1], repeated_variables[j+1], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                self.add(repeated_variables[j], repeated_variables[j+1], v, v, -4*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']

        # Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write(
                "QPU access time:" + str(qpu_access_time * 10 ** (-6)) + "s." + " ; QPU sampling time:" + str(
                    qpu_sampling_time * 10 ** (-6)) + "s.\n")
            archivo.write("o " + utils.count_unsatisfied_clauses(response.first.sample, self.clauses) + "\n")
            archivo.write("e " + str(response.first.energy) + "\n")
            archivo.write("v " + str(response.first.sample) + "\n")
            archivo.write(str(response.samples))

    def compute_bqm(self):
        if not self.Q:
            self.fillQ()
        bqm = dimod.BinaryQuadraticModel.from_qubo(self.Q, offset=0).change_vartype("SPIN", True)
        self.bqm = bqm
        return bqm

    def compute_bqm_embedded(self, token):
        couplings = []
        if not self.Q:
            self.fillQ()
        for k in self.Q.keys():
            if k[0] != k[1]:
                couplings.append(k)
        embedding = find_embedding(couplings, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
        if not self.bqm:
            self.compute_bqm()
        bqm_embedded = dwave.embedding.embed_bqm(source_bqm=self.bqm, embedding=embedding,
                                                              target_adjacency=DWaveSampler(token=token).adjacency, chain_strength=2)
        self.bqm_embedded['bqm_embedded'] = bqm_embedded
        self.bqm_embedded['embedding'] = embedding
        return {'bqm_embedded': bqm_embedded, 'embedding': embedding}

    def compute_scaled_bqm_embedded(self, token):
        if not self.bqm_embedded:
            self.compute_bqm_embedded(token)
        if not self.scaled_bqm_embedded:
            scale_factor = utils.compute_scale_factor(self.bqm_embedded['bqm_embedded'])
            self.scale_factor = scale_factor
            temp_bqm_embedded = self.bqm_embedded['bqm_embedded'].copy()
            if scale_factor > 1:
                temp_bqm_embedded.scale(1 / scale_factor)
            self.scaled_bqm_embedded['scaled_bqm_embedded'] = temp_bqm_embedded
            self.scaled_bqm_embedded['scale_factor'] = scale_factor
            self.scale_factor = scale_factor
        return {'scaled_bqm_embedded': self.scaled_bqm_embedded['scaled_bqm_embedded'], 'scale_factor': self.scale_factor}

    def solve_dwave_from_scaled_bqm_embedded(self, token, num_reads=100, annealing_time=100):
        if not self.scaled_bqm_embedded:
            self.compute_scaled_bqm_embedded(token)
        if not self.response_dwave_from_scaled_bqm_embedded:
            embedding_refactored = utils.refactor_embedding_1_to_1(self.scaled_bqm_embedded['scaled_bqm_embedded'])
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=embedding_refactored)
            response = sampler.sample(self.scaled_bqm_embedded['scaled_bqm_embedded'], num_reads=num_reads, annealing_time=annealing_time,
                                      return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            response_aux = response.to_serializable()
            self.response_dwave_from_scaled_bqm_embedded = dimod.SampleSet.from_serializable(response_aux)
        return self.response_dwave_from_scaled_bqm_embedded


class CJ2_bian(Instance):
    def __init__(self, clauses, N):
        self.clauses = clauses
        self.N = N
        self.b = int(3*len(self.clauses))
        self.Q = {}
        self.bian_embedding = {}
        self.sign = []
        self.bqm = None
        self.bqm_embedded = {}
        self.scaled_bqm_embedded = {}
        self.scale_factor = None
        self.response_dwave_from_scaled_bqm_embedded = None
        self.all_solution_space = {}

    def add(self, a, b, x, y, value):
        if a not in self.bian_embedding.keys():
            self.bian_embedding[a] = x
        if b not in self.bian_embedding.keys():
            self.bian_embedding[b] = y
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
            self.add(k, k, var[0], var[0], s[0]-s[0]*s[2])
            self.add(k+1, k+1, var[1], var[1], -2*s[1])
            self.add(k+2, k+2, var[2], var[2], s[2]-s[0]*s[2])
            self.add(self.b+i, self.b+i, self.N+i+1, self.N+i+1, -1+s[0]-s[1]+s[2])
            self.add(k, k+2, var[0], var[2], 2*s[0]*s[2])
            self.add(k, self.b+i, var[0], self.N+i+1, -2*s[0])
            self.add(k+1, self.b+i, var[1], self.N+i+1, 2*s[1])
            self.add(k+2, self.b+i, var[2], self.N+i+1, -2*s[2])
            k+=3
            self.sign = self.sign + s

        for v in range(1,self.N+1):
            repeated_variables = sorted([i for i, variable in self.bian_embedding.items() if variable==v])
            for j in range(len(repeated_variables)-1):
                self.add(repeated_variables[j], repeated_variables[j], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                self.add(repeated_variables[j+1], repeated_variables[j+1], v, v, 2*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                self.add(repeated_variables[j], repeated_variables[j+1], v, v, -4*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])

    def solve(self, file_name):
        self.fillQ()
        response = utils.solve_with_DWave(self.Q)
        qpu_access_time = response.info['timing']['qpu_access_time']
        qpu_sampling_time = response.info['timing']['qpu_sampling_time']

        # Write the response to a file
        with open(file_name, 'w') as archivo:
            archivo.write(
                "QPU access time:" + str(qpu_access_time * 10 ** (-6)) + "s." + " ; QPU sampling time:" + str(
                    qpu_sampling_time * 10 ** (-6)) + "s.\n")
            archivo.write("o " + utils.count_unsatisfied_clauses(response.first.sample, self.clauses) + "\n")
            archivo.write("e " + str(response.first.energy) + "\n")
            archivo.write("v " + str(response.first.sample) + "\n")
            archivo.write(str(response.samples))

    def compute_bqm(self):
        if not self.Q:
            self.fillQ()
        bqm = dimod.BinaryQuadraticModel.from_qubo(self.Q, offset=0).change_vartype("SPIN", True)
        self.bqm = bqm
        return bqm

    def compute_bqm_embedded(self, token):
        couplings = []
        if not self.Q:
            self.fillQ()
        for k in self.Q.keys():
            if k[0] != k[1]:
                couplings.append(k)
        embedding = find_embedding(couplings, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
        if not self.bqm:
            self.compute_bqm()
        bqm_embedded = dwave.embedding.embed_bqm(source_bqm=self.bqm, embedding=embedding,
                                                              target_adjacency=DWaveSampler(token=token).adjacency, chain_strength=2)
        self.bqm_embedded['bqm_embedded'] = bqm_embedded
        self.bqm_embedded['embedding'] = embedding
        return {'bqm_embedded': bqm_embedded, 'embedding': embedding}

    def compute_scaled_bqm_embedded(self, token):
        if not self.bqm_embedded:
            self.compute_bqm_embedded(token)
        if not self.scaled_bqm_embedded:
            scale_factor = utils.compute_scale_factor(self.bqm_embedded['bqm_embedded'])
            self.scale_factor = scale_factor
            temp_bqm_embedded = self.bqm_embedded['bqm_embedded'].copy()
            if scale_factor > 1:
                temp_bqm_embedded.scale(1 / scale_factor)
            self.scaled_bqm_embedded['scaled_bqm_embedded'] = temp_bqm_embedded
            self.scaled_bqm_embedded['scale_factor'] = scale_factor
            self.scale_factor = scale_factor
        return {'scaled_bqm_embedded': self.scaled_bqm_embedded['scaled_bqm_embedded'], 'scale_factor': self.scale_factor}

    def solve_dwave_from_scaled_bqm_embedded(self, token, num_reads=100, annealing_time=100):
        if not self.scaled_bqm_embedded:
            self.compute_scaled_bqm_embedded(token)
        if not self.response_dwave_from_scaled_bqm_embedded:
            embedding_refactored = utils.refactor_embedding_1_to_1(self.scaled_bqm_embedded['scaled_bqm_embedded'])
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=embedding_refactored)
            response = sampler.sample(self.scaled_bqm_embedded['scaled_bqm_embedded'], num_reads=num_reads, annealing_time=annealing_time,
                                      return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            response_aux = response.to_serializable()
            self.response_dwave_from_scaled_bqm_embedded = dimod.SampleSet.from_serializable(response_aux)
        return self.response_dwave_from_scaled_bqm_embedded