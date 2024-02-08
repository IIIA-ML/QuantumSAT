# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
import random


def ksat_3sat(clauses, b):
    new_clauses=[]
    for clause in clauses:
        if len(clause)!=3:
            new_clauses.append(np.array([clause[0],clause[1],b]))
            if len(clause)>3:
                for i in range(2, len(clause)-2):
                    b+=1
                    new_clauses.append(np.array([-(b-1),clause[i],b]))
            new_clauses.append(np.array([-b,clause[len(clause)-2], clause[len(clause)-1]]))
                
        else:
            new_clauses.append(np.array(clause))
    return new_clauses, b+1


def tseitin(clauses,b): #Al solucionar 3SAT creamos una nueva variable extra
    h={}
    J={}

    clauses, b = ksat_3sat(clauses, b)
    b_init=b
    
    for clause in clauses:
        sign = [1 if x >= 0 else -1 for x in sorted(clause, key=abs)]
        var = [abs(x) for x in clause]
        
        for v in var:
            if v not in h.keys():
                h[v]=0.0
        h[b]=0.0
        if tuple(sorted((var[0],var[1]))) not in J.keys(): J[tuple(sorted((var[0],var[1])))]=0.0
        if tuple(sorted((var[1],var[2]))) not in J.keys(): J[tuple(sorted((var[1],var[2])))]=0.0
        if tuple(sorted((var[0],var[2]))) not in J.keys(): J[tuple(sorted((var[0],var[2])))]=0.0
        J[var[0],b]=0.0
        J[var[1],b]=0.0
        J[var[2],b]=0.0
        
        #h y J debido a x1 ∨ x2 ↔ b
        h[var[0]]+= (0.5 if sign[0]>0 else -0.5)
        h[var[1]]+= (0.5 if sign[1]>0 else -0.5)
        h[b]-=1
    
        J[tuple(sorted((var[0],var[1])))]+= (0.5 if ((sign[0]>0 and sign[1]>0) or (sign[0]<0 and sign[1]<0)) else -0.5)
        J[var[0],b]-= (1 if sign[0]>0 else -1)
        J[var[1],b]-= (1 if sign[1]>0 else -1)
        
        #h y J debido a b ∨ x3
        h[b]-=0.5
        h[var[2]]-= (0.5 if sign[2]>0 else -0.5)
        J[var[2],b]+= (0.5 if sign[2]>0 else -0.5)
        
        b+=1
    return h,J,(b-1)-b_init


def three_sat_max2xor(clauses, b):
    h={}
    J={}

    clauses, b = ksat_3sat(clauses, b)
    b_init=b
    
    for clause in clauses:
        sign = [1 if x >= 0 else -1 for x in sorted(clause, key=abs)]
        var = [abs(x) for x in clause]
        
        for v in var:
            if v not in h.keys():
                h[v]=0.0
        h[b]=0.0
        if tuple(sorted((var[0],var[1]))) not in J.keys(): J[tuple(sorted((var[0],var[1])))]=0.0
        if tuple(sorted((var[1],var[2]))) not in J.keys(): J[tuple(sorted((var[1],var[2])))]=0.0
        if tuple(sorted((var[0],var[2]))) not in J.keys(): J[tuple(sorted((var[0],var[2])))]=0.0
        J[var[0],b]=0.0
        J[var[1],b]=0.0
        J[var[2],b]=0.0
        
        """h[var[0]]-=0.5
        h[var[2]]-=0.5
        J[var[0],var[2]]+=0.5
        h[var[0]]-=0.5
        h[b]+=0.5
        J[var[0],b]-=0.5
        h[var[2]]-=0.5
        h[b]+=0.5
        J[var[2],b]-=0.5
        h[var[1]]-=1 #MIRAR ESTA PARTE!!!!!!!
        h[b]-=1
        J[var[1],b]+=1"""
        h[b]-=0.5
        h[var[1]]-=sign[1]*0.5
        J[tuple(sorted((var[0],var[2])))]+=sign[0]*sign[2]*0.5
        J[var[0],b]-=sign[0]*0.5
        J[var[2],b]-=sign[2]*0.5
        J[var[1],b]+=sign[1]*0.5

        b+=1

    return h,J,(b-1)-b_init    


def regular_like(clauses, b): #From ksat to Max2XOR
    h={}
    J={}

    b_init=b
    random.seed(42)
    
    for clause in clauses:
        sign = [1 if x >= 0 else -1 for x in clause]
        var = random.sample([abs(x) for x in clause], len(clause))
        
        v_aux = list(range(b,b+len(var)-1))
        
        for i in range(len(var)):
            """if var[i] not in h.keys():
                h[var[i]]=0.0"""
            for j in range(i+1, len(var)):
                if tuple(sorted((var[i],var[j]))) not in J.keys():
                    J[tuple(sorted((var[i],var[j])))]=0.0
            #k-2 Auxuliary variables
            for aux in range(b, b+len(var)-1):
                J[var[i],aux]=0.0
                if aux!=b+len(var)-2:
                    J[aux,aux+1]=0.0

        k=0
        J[tuple(sorted((var[0],var[1])))]+=sign[0]*sign[1]*0.5
        J[var[0],v_aux[k]]-=sign[0]*0.5
        J[var[1],v_aux[k]]-=sign[1]*0.5
        if len(var)>2:
            for i in range(2,len(var)):
                k+=1
                J[var[i],v_aux[k-1]]+=sign[i]*0.5
                J[v_aux[k-1],v_aux[k]]-=0.5
                J[var[i],v_aux[k]]-=sign[i]*0.5

        b+=(len(v_aux))
    return h,J,(b-1)-b_init


def tree_like(clauses, b): #Using only the 1st and 2nd cases
    h={}
    J={}

    b_init=b

    #Sort literals by studying the relation between them in all clauses
    weights = {}
    for clause in clauses:
        for i in range(len(clause)):
            for j in range(i+1,len(clause)):
                if tuple(sorted((clause[i],clause[j]))) not in weights.keys():
                    weights[tuple(sorted((abs(clause[i]),abs(clause[j]))))]=1/len(clause)
                else:
                    weights[tuple(sorted((abs(clause[i]),abs(clause[j]))))]+=1/len(clause)
    for clause in clauses:
        sign = {}
        for x in clause:
            sign[abs(x)] = (1 if x >= 0 else -1)
        var = [abs(x) for x in clause]
        coup=[]
        for i in range(len(var)):
            for j in range(i+1,len(var)):
                coup.append(tuple(sorted((var[i],var[j]))))
        sorted_coup = sorted(coup, key=lambda x: weights[x], reverse=True)
        new_list=[sorted_coup[0][0]*sign[sorted_coup[0][0]], sorted_coup[0][1]*sign[sorted_coup[0][1]]]
        for x in var:
            if x*sign[x] not in new_list:
                new_list.append(x*sign[x])

        sign = [1 if x >= 0 else -1 for x in new_list]
        var = [abs(x) for x in new_list]
        
        v_aux = list(range(b,b+len(var)))
        
        for i in range(len(var)):
            for j in range(i+1, len(var)):
                if tuple(sorted((var[i],var[j]))) not in J.keys():
                    J[tuple(sorted((var[i],var[j])))]=0.0
            #Auxuliary variables
            for aux in range(b, b+len(var)-1):
                J[var[i],aux]=0.0
                if aux!=b+len(var)-2:
                    J[aux,aux+1]=0.0

        if len(var)==2:
            J[var[0],var[1]]+=0.5*sign[0]*sign[1]
            J[var[0],b]-=0.5*sign[0]
            J[var[1],b]-=0.5*sign[1]
            b+=1
        if len(var)>=3:
            J[var[0],var[1]]+=0.5*sign[0]*sign[1]
            J[var[0],b]-=0.5*sign[0]
            J[var[1],b]-=0.5*sign[1]
            b+=1
            for i in range(2,len(var)):
                J[var[i], b-1]+=0.5*sign[i]
                J[b-1,b]-=0.5
                J[var[i],b]-=0.5*sign[i]
                b+=1

        b+=(len(v_aux))
    
    return h,J,(b-1)-b_init


def clique(clauses, b): #Only for kSAT clauses where k%2
    h={}
    J={}

    b_init=b
    
    for clause in clauses:
        sign = [1 if x >= 0 else -1 for x in clause]
        var = [abs(x) for x in clause]
        
        v_aux = list(range(b,b+int(math.log2(len(var)))))
        
        for i in range(len(var)):
            if var[i] not in h.keys():
                h[var[i]]=-0.5*sign[i]
            else:
                h[var[i]]-=0.5*sign[i]
            for j in range(i+1, len(var)):
                if tuple(sorted((var[i],var[j]))) not in J.keys():
                    J[tuple(sorted((var[i],var[j])))]=0.5*sign[i]*sign[j]
                else:
                    J[tuple(sorted((var[i],var[j])))]+=0.5*sign[i]*sign[j]
            #k-2 Auxuliary variables
            for k in range(len(v_aux)):
                h[v_aux[k]]=-2**(k-1)
                J[var[i],v_aux[k]]=2**(k-1)*sign[i]
                for k_2 in range(k+1, len(v_aux)):
                    J[v_aux[k],v_aux[k_2]]=2**(k+k_2-1)

        b+=(len(v_aux))    
    
    return h,J,(b-1)-b_init