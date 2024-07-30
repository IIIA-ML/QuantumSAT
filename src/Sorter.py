class Sorter:
    def Sort(self, clases):
        pass

class Absolute(Sorter):
    def Sort(self, clauses):
        new_clauses = []
        for c in clauses:
            new_clauses.append(sorted(c, key=abs))
        return new_clauses

class CJ2(Sorter):
    def Sort(self, clauses):
        new_clauses = []
        couplings = []
        for c in clauses:
            c = sorted(c, key=abs)
            if (abs(c[0]),abs(c[2])) not in couplings and (abs(c[1]),abs(c[2])) in couplings:
                new_clauses.append([c[1], c[0], c[2]])
            elif (abs(c[0]),abs(c[2])) not in couplings and (abs(c[0]),abs(c[1])) in couplings:
                new_clauses.append([c[0], c[2], c[1]])
            else:
                new_clauses.append(c)
                if (abs(c[0]), abs(c[2])) not in couplings:
                    couplings.append((abs(c[0]), abs(c[2])))
        return new_clauses