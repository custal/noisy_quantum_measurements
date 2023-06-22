import sympy as sp
import numpy as np
from sympy import sqrt, Symbol
from sympy import nsimplify as simp
from sympy.matrices import Matrix
from sympy.physics.quantum import TensorProduct, Dagger
from itertools import combinations
import qutip as qp
from multiprocessing.pool import Pool
import tqdm
import datetime


class Qubit:
    qubit_total = 0

    @classmethod
    def count_qubit(cls):
        cls.qubit_total += 1

    def __init__(self, a, b):
        self.qubit_num = self.qubit_total
        self.count_qubit()

        self.a = sp.nsimplify(a)
        self.b = sp.nsimplify(b)

        self.state = Matrix([self.a, self.b])

    @property
    def dag(self):
        return Dagger(self.state)

    #     def __mul__(self, other):
    #         return Qubitself.state * other

    def __repr__(self):
        return repr(self.state)


def matrix_to_braket(matrix, qubit_num=None):
    cols = matrix.cols
    rows = matrix.rows

    def binary(x):
        if qubit_num is None:
            return x
        else:
            num = f"{x:0>b}"
            return (qubit_num - len(num)) * '0' + num

    col_symbols = [sp.Symbol(f"|{binary(i)}>", commutative=False) for i in range(cols)]
    row_symbols = [sp.Symbol(f"<{binary(i)}|", commutative=False) for i in range(rows)]

    expression = 0

    for i, element in enumerate(matrix):
        row = i // cols
        col = i % cols

        expression += element * col_symbols[col] * row_symbols[row]

    return expression


def search_all_combinations(probabilities):
    """Find the summation of all possible combinations of probabilities"""
    max_len = int(len(probabilities) / 2)
    # Add empty set edge case
    group_probabilities = {frozenset(): 0, frozenset(probabilities.keys()): 1}

    for i in range(max_len):
        combs = combinations(probabilities, i + 1)
        for selection in combs:
            duplicate_probs = dict(probabilities)
            for prob in selection:
                duplicate_probs.pop(prob)
            list1 = [probabilities.get(key) for key in selection]
            list2 = list(duplicate_probs.values())

            group_probabilities[frozenset(selection)] = sum(list1)
            group_probabilities[frozenset(duplicate_probs.keys())] = sum(list2)

    return group_probabilities


class basic_measure:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __call__(self, calc_prob):
        key = calc_prob[0]
        item = calc_prob[1]
        prob = item - self.a * sp.conjugate(self.a)
        prob = prob.subs(self.b, sqrt(1 - self.a * sp.conjugate(self.a))).evalf()
        measure = sqrt(sp.Integral(prob ** 2, (self.a, 0, 1))).evalf(chop=True)

        return key, measure


class root_mean_square_measure:
    def __init__(self, a, b):
        self.a = a
        self.b = b

        self.theta = sp.Symbol("θ", real=True)
        self.phi = sp.Symbol("𝜙", real=True)

    def __call__(self, calc_prob):
        key = calc_prob[0]
        item = calc_prob[1]

        prob = item - self.a * sp.conjugate(self.a)
        prob = prob.subs({self.a: sp.cos(self.theta / 2), self.b: sp.E ** (sp.I * self.phi) * sp.sin(self.theta / 2)})
        measure = sp.simplify(
            sqrt(sp.Integral(sp.sin(self.theta) / (4 * sp.pi) * prob ** 2, (self.theta, 0, sp.pi), (self.phi, 0, 2 *
                                                                                                    sp.pi)))).evalf(
            chop=True)

        return key, measure

def measure_all_probabilities(group_probabilities, a, b, measure_func=basic_measure):
    print(f"Calculating {len(groups)} measures...")

    measures = {}
    for i, (key, item) in enumerate(group_probabilities.items()):
        if i % 10 == 0 and i != 0:
            print(f"Calculated {i} measures...")
        measure = measure_func(item, a, b)
        measures[key] = measure

    return measures


def lowest_measure(group_measures):
    lowest_key = set()
    lowest_val = np.inf
    for key, item in group_measures.items():
        if item <= lowest_val:
            if item == lowest_val:
                lowest_key.add(key)
            else:
                lowest_key = {key}
            lowest_val = item

    return lowest_key, lowest_val


def make_tensor(povm):
    tensor = {}
    for i, p1 in enumerate(povm):
        for j, p2 in enumerate(povm):
            tensor[f"P({i},{j})"] = TensorProduct(p1, p2)

    return tensor


def calc_probabilities(tensor, a, b):
    psi = Matrix([0] * tensor.cols)
    psi[0] = a
    psi[-1] = b

    return (Dagger(psi) * tensor * psi)[0]


if __name__ == "__main__":
    phi0 = Qubit(1, 0)
    phi1 = Qubit(1 / sqrt(3), sqrt(2 / 3))
    phi2 = Qubit(1 / sqrt(3), sp.E ** (sp.I * 2 * sp.pi / 3) * sqrt(2 / 3))
    phi3 = Qubit(1 / sqrt(3), sp.E ** (-sp.I * 2 * sp.pi / 3) * sqrt(2 / 3))

    M0 = simp(1 / 2) * phi0.state * phi0.dag
    M1 = simp(1 / 2) * phi1.state * phi1.dag
    M2 = simp(1 / 2) * phi2.state * phi2.dag
    M3 = simp(1 / 2) * phi3.state * phi3.dag

    tetra = [M0, M1, M2, M3]

    trine_tensor = make_tensor(tetra)

    a = sp.Symbol("a")
    b = sp.Symbol("b")

    probabilities = {}
    for key, tri in trine_tensor.items():
        probabilities[key] = calc_probabilities(tri, a, b)

    groups = search_all_combinations(probabilities)

    print(f"Calculating {len(groups)} measures...")

    all_measures = {}
    with Pool(processes=4) as pool:
        for output in tqdm.tqdm(pool.imap_unordered(root_mean_square_measure(a, b), groups.items()), total=len(groups)):
            all_measures[output[0]] = output[1]

    group, measure = lowest_measure(all_measures)

    with open(f"two_tetrahedron_povm_measurements_{datetime.datetime.now().strftime('%Y_%m_%d-%H_%M_%S')}.txt", "w") as f:
        f.write(str(all_measures))

    print(all_measures)
    print(group, measure)