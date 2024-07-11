import dimod
from dwave.embedding import embed_bqm
from minorminer import find_embedding
from dwave.system import DWaveSampler, FixedEmbeddingComposite, EmbeddingComposite
from dwave.samplers import SimulatedAnnealingSampler
import utils


class Solver:
    def __init__(self):
        self.response = []
    
    def Solve(self):
        pass

class D_Wave(Solver):
    def compute_scale_factor(self, bqm_embedded):
        J_max = max(bqm_embedded.quadratic.values())
        J_min = min(bqm_embedded.quadratic.values())
        h_max = max(bqm_embedded.linear.values())
        h_min = min(bqm_embedded.linear.values())
        J_per_qubit = {}
        for key, value in bqm_embedded.quadratic.items():
            if key[0] in J_per_qubit.keys():
                J_per_qubit[key[0]]+=value
            else:
                J_per_qubit[key[0]]=value
            if key[1] in J_per_qubit.keys():
                J_per_qubit[key[1]]+=value
            else:
                J_per_qubit[key[1]]=value
        coupling_limit = max(max(max(J_per_qubit.values())/15,0),max(min(J_per_qubit)/(-18),0))
        scale_factor = max(max(h_max/4,0),max(h_min/(-4),0),max(J_max/1,0),max(J_min/(-2),0),coupling_limit)
        if scale_factor <= 1:
            scale_factor = 1
        return scale_factor

    def refactor_embedding_1_to_1(self, bqm):
        embedding_bqm = {}
        for k in bqm.linear.keys():
            if k not in embedding_bqm.keys():
                embedding_bqm[k] = [k]
        for k in bqm.quadratic.keys():
            if k[0] not in embedding_bqm.keys():
                embedding_bqm[k[0]] = [k[0]]
            if k[1] not in embedding_bqm.keys():
                embedding_bqm[k[1]] = [k[1]]
        return embedding_bqm

    
    def compute_bqm_embedded(self, Q, token):
        self.bqm = dimod.BinaryQuadraticModel.from_qubo(Q, offset=0).change_vartype("SPIN", True)
        couplings = []
        for k in list(self.bqm.quadratic.keys()):
            if k[0] != k[1]:
                couplings.append(k)
        self.embedding = find_embedding(couplings, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
        self.bqm_embedded = embed_bqm(source_bqm=self.bqm, embedding=self.embedding, target_adjacency=DWaveSampler(token=token).adjacency)

    def compute_scaled_bqm_embedded(self, token):
        scale_factor = self.compute_scale_factor(self.bqm_embedded)
        temp_bqm_embedded = self.bqm_embedded.copy()
        if scale_factor > 1:
            temp_bqm_embedded.scale(1 / scale_factor)
        self.scaled_bqm_embedded = temp_bqm_embedded
        self.scale_factor = scale_factor
        
        
    def Solve(self, Q, token, qubit_level):
        self.qubit_level = qubit_level
        if token==None:
            raise Exception("No token entered for solving with D-Wave QPU")
        if self.qubit_level==True:
            self.compute_bqm_embedded(Q, token)
            self.compute_scaled_bqm_embedded(token)
            embedding_refactored = self.refactor_embedding_1_to_1(self.scaled_bqm_embedded)
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=embedding_refactored)
            response = sampler.sample(self.scaled_bqm_embedded, num_reads=100, annealing_time=100, return_embedding=True, reduce_intersample_correlation=True, auto_scale=False)
            response_aux = response.to_serializable()
            self.response = dimod.SampleSet.from_serializable(response_aux)
        else:
            self.embedding = find_embedding(Q, DWaveSampler(token=token).edgelist, random_seed=42, verbose=0)
            sampler = FixedEmbeddingComposite(DWaveSampler(token=token), embedding=self.embedding)
            response = sampler.sample_qubo(Q, num_reads=100, annealing_time=100, return_embedding=True, reduce_intersample_correlation=True)
            response_aux = response.to_serializable()
            self.response = dimod.SampleSet.from_serializable(response_aux)
        
        return self.response, self.embedding

    def SAT_solution(self, real_encoding, clauses):
        solutions_bqm_embedded_df_aux = utils.unchain_dwave_solutions(self.response, real_encoding)
        solutions_bqm_embedded_df = solutions_bqm_embedded_df_aux.loc[solutions_bqm_embedded_df_aux.index.repeat(solutions_bqm_embedded_df_aux.Occurrences)].reset_index(drop=True).drop(columns=['Occurrences'])
        self.o = {}
        for _, row in solutions_bqm_embedded_df.iterrows():
            assigment = dict(row)
            o = utils.count_unsatisfied_clauses(assigment, clauses)
            if o not in self.o.keys():
                self.o[o] = 1
            else:
                self.o[o] += 1
        return self.o
        

    

class Simulated_annealing(Solver):
    def Solve(self, Q, token, qubit_level):
        sampler = SimulatedAnnealingSampler()
        self.response = sampler.sample_qubo(Q, num_reads=100)

        return self.response, {}




    