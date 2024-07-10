class Joiner:
    def __init__(self):
        self.final_encoding = {}
        self.final_Q = {}

    def get_variables(self):
        variables = []
        for mapper in self.mappers: #mapper = [Q's, encoding's]
            for key, var in mapper[1].items():
                if var!=None and var not in variables:
                    variables.append(var)
        return sorted(variables)
        
    def Join(self, mappers):
        self.mappers = mappers
        pass
        
class Bian(Joiner):
    def get_index_aux(self):
        index_aux = 0
        for mapper in self.mappers: #mapper = [Q's, encoding's]
            for k, v in mapper[1].items():
                if v!=None:
                    index_aux += 1
        return index_aux #Index necessary for the encoding of bian
        
    def Join(self, mappers):
        self.mappers = mappers
        b = self.get_index_aux()
        literal = 0
        final_encoding = {}
        final_Q = {}
        for mapper in self.mappers:
            Q = mapper[0]
            encoding = mapper[1]
            len_literals = len(set([l for l,value in encoding.items() if value!=None]))
            for k, v in Q.items():
                key = [0,0]
                if encoding[k[0]]!=None:
                    key[0] = literal+k[0]
                else:
                    key[0] = b+k[0]-len_literals
                if encoding[k[1]]!=None:
                    key[1] = literal+k[1]
                else:
                    key[1] = b+k[1]-len_literals

                if key[0] not in final_encoding.keys():
                    final_encoding[key[0]] = encoding[k[0]]
                if key[1] not in final_encoding.keys():
                    final_encoding[key[1]] = encoding[k[1]]

                final_Q[(key[0], key[1])] = v
                
            literal += len_literals
            b += (len(encoding)-len_literals)
            
        variables = self.get_variables()
        for v in variables:
            repeated_variables = [index for index, variable in final_encoding.items() if variable==v]
            for j in range(len(repeated_variables)-1):
                if (repeated_variables[j], repeated_variables[j+1]) not in final_Q.keys(): #Only add coupling if those literals are not already coupled
                    if (repeated_variables[j], repeated_variables[j]) not in final_Q.keys():
                        final_Q[(repeated_variables[j], repeated_variables[j])] = 4 #*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                    else:
                        final_Q[(repeated_variables[j], repeated_variables[j])] += 4 #*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                    if (repeated_variables[j+1], repeated_variables[j+1]) not in final_Q.keys():
                        final_Q[(repeated_variables[j+1], repeated_variables[j+1])] = 4 #*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                    else:
                        final_Q[(repeated_variables[j+1], repeated_variables[j+1])] += 4 #*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])
                    if (repeated_variables[j], repeated_variables[j+1]) not in final_Q.keys():
                        final_Q[(repeated_variables[j], repeated_variables[j+1])] = -8 #*self.sign[repeated_variables[j]]*self.sign[repeated_variables[j+1]])

        return final_Q, final_encoding
        

class SAT_variables(Joiner):
    def get_index_aux(self):
        index_aux = 0
        var_list = []
        for mapper in self.mappers: #mapper = [Q's, encoding's]
            list_v = [v for v in mapper[1].values() if v is not None]
            if max(list_v) > index_aux:
                index_aux = max(list_v)
            #for k, v in mapper[1].items():
            #    if v!=None and v not in var_list:
            #        index_aux += 1
            #        var_list.append(v)
        return index_aux #Index necessary for the encoding of bian
    
    def Join(self, mappers):
        self.mappers = mappers
        b = self.get_index_aux() #key after studying all the mappers submitted
        final_encoding = {}
        final_Q = {}
        for mapper in self.mappers:
            Q = mapper[0]
            encoding = mapper[1]
            len_literals = len(set([l for l,value in encoding.items() if value!=None]))
            for k, v in Q.items():
                key = [0,0]
                if encoding[k[0]]!=None:
                    key[0] = encoding[k[0]]-1
                else:
                    key[0] = b+k[0]-len_literals
                if encoding[k[1]]!=None:
                    key[1] = encoding[k[1]]-1
                else:
                    key[1] = b+k[1]-len_literals

                if key[0] not in final_encoding.keys():
                    final_encoding[key[0]] = encoding[k[0]]
                if key[1] not in final_encoding.keys():
                    final_encoding[key[1]] = encoding[k[1]]

                if (key[0],key[1]) not in final_Q.keys():
                    final_Q[(key[0], key[1])] = v
                else:
                    final_Q[(key[0], key[1])] += v
            b += (len(encoding)-len_literals)

        return final_Q, final_encoding #final_encoding = {Q_variable: SAT_variable}
        