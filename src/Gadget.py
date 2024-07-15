class Gadget:
    def __init__(self):
        self.Q = {}
        self.encoding = {}
        
    def get_index_auxiliar(self, clauses):
        num_literals = 0
        for c in clauses:
            num_literals += len(c)
        self.b = num_literals
    
    def add(self, a, b, x, y, value):
        if a not in self.encoding.keys():
            self.encoding[a] = x
        if b not in self.encoding.keys():
            self.encoding[b] = y
        if a > b:
            a,b = b,a
        if (a,b) in self.Q.keys():
            self.Q[(a,b)] += value
        else:
            self.Q[(a,b)] = value

    def FillQ(self, clauses):
        pass
    
class CJ1(Gadget):
    def FillQ(self, clauses):
        self.get_index_auxiliar(clauses) #auxiliar will be the num_literals+1
        literal = 0
        for c in clauses:
            s = [1 if l>0 else -1 for l in c] #sign of the literal (negated or not)
            var = [abs(l) for l in c]
            self.add(literal, literal, var[0], var[0], 3*s[0]-s[0]*s[1])
            self.add(literal+1, literal+1, var[1], var[1], 3*s[1]-s[0]*s[1])
            self.add(literal+2, literal+2, var[2], var[2], -2*s[2])
            self.add(self.b, self.b, None, None, -3+2*s[0]+2*s[1]-s[2])
            self.add(literal, literal+1, var[0], var[1], 2*s[0]*s[1])
            self.add(literal, self.b, var[0], None, -4*s[0])
            self.add(literal+1, self.b, var[1], None, -4*s[1])
            self.add(literal+2, self.b, var[2], None, 2*s[2])
            
            literal += 3
            self.b += 1

        return [self.Q, self.encoding]

class CJ2(Gadget):
    def FillQ(self, clauses):
        self.get_index_auxiliar(clauses) #auxiliar will be the num_literals+1
        literal = 0
        for c in clauses:
            s = [1 if l>0 else -1 for l in c] #sign of the literal (negated or not)
            var = [abs(l) for l in c]
            self.add(literal, literal, var[0], var[0], s[0] - s[0] * s[2])
            self.add(literal+1, literal+1, var[1], var[1], -2 * s[1])
            self.add(literal+2, literal+2, var[2], var[2], s[2] - s[0] * s[2])
            self.add(self.b, self.b, None, None, -1 + s[0] - s[1] + s[2])
            self.add(literal, literal+2, var[0], var[2], 2 * s[0] * s[2])
            self.add(literal, self.b, var[0], None, -2 * s[0])
            self.add(literal+1, self.b, var[1], None, 2 * s[1])
            self.add(literal+2, self.b, var[2], None, -2 * s[2])
            
            literal += 3
            self.b += 1

        return [self.Q, self.encoding]

class Choi(Gadget):
    def FillQ(self, clauses):
        L = []
        for c in clauses:
            L.extend(c)
        for i in range(len(L)):
            for j in range(len(L)):
                if i > j:
                    continue
                if i == j:
                    self.add(i, j, abs(L[i]), abs(L[j]), -1)
                elif j - i <= 2 and j//3 == i//3:
                    self.add(i, j, abs(L[i]), abs(L[j]), 3)
                elif abs(L[i]) == abs(L[j]) and L[i] != L[j]:
                    self.add(i, j, abs(L[i]), abs(L[j]), 3)