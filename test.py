from sympy import sqrt
from sympy.physics.quantum import Bra, Ket
from sympy import Matrix
import sympy as sp
import qutip as qp
from multiprocessing.pool import Pool
from time import sleep
from itertools import combinations
from sympy import nsimplify as simp
import tqdm
from wolframclient.evaluation import WolframLanguageAsyncSession
import sympy as sp
import numpy as np
from sympy import sqrt, Symbol, mathematica_code
from sympy import nsimplify as simp
from sympy.matrices import Matrix
from sympy.physics.quantum import TensorProduct, Dagger
from itertools import combinations
import qutip as qp
from multiprocessing.pool import Pool
from wolframclient.evaluation import WolframLanguageSession, WolframCloudSession
from wolframclient.language import wlexpr
import datetime
import time
from sympy.parsing import mathematica

def test(key):
    sleep(2)
    print(1)

    return (key[0], key[1])


def search_all_combinations(group):
    """Find the summation of all possible combinations of groups"""
    max_len = int(len(group) / 2)
    # Add empty set edge case
    sub_groups = {frozenset(): 0, frozenset(group.keys()): 1}

    for i in range(max_len):
        combs = combinations(group, i + 1)
        for selection in combs:
            duplicate_group = dict(group)
            for prob in selection:
                duplicate_group.pop(prob)
            list1 = [group.get(key) for key in selection]
            list2 = list(duplicate_group.values())

            print(i)
            sub_groups[frozenset(selection)] = sum(list1, sp.zeros(4))
            sub_groups[frozenset(duplicate_group.keys())] = sum(list2, sp.zeros(4))

    return sub_groups


def basic_measure(calc_prob, a, b):
    prob = calc_prob - a * sp.conjugate(a)
    prob = prob.subs(b, sqrt(1 - a * sp.conjugate(a))).evalf()
    measure = sqrt(sp.Integral(prob ** 2, (a, 0, 1)))

    return measure


def measure_all_probabilities(group_probabilities, a, b, measure_func=basic_measure):
    print(f"Calculating {len(groups)} measures...")

    measures = {}
    for i, (key, item) in enumerate(group_probabilities.items()):
        if i % 10 == 0 and i != 0:
            print(f"Calculated {i} measures...")
        measure = sp.simplify(measure_func(item, a, b)).evalf(chop=True)
        measures[key] = measure

    return measures


def lowest_measure(group_measures):
    lowest_val = np.inf
    for key, item in group_measures.items():
        if item <= lowest_val:
            if item == lowest_val:
                lowest_key.add(key)
            else:
                lowest_key = {key}
            lowest_val = item

    return lowest_key, lowest_val


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

def make_tensor(povm):
    tensor = {}
    for i, p1 in enumerate(povm):
        for j, p2 in enumerate(povm):
            tensor[f"({i},{j})"] = TensorProduct(p1, p2)

    return tensor


def calc_probabilities(tensor, a, b):
    psi = Matrix([0] * tensor.cols)
    psi[0] = a
    psi[-1] = b

    return (Dagger(psi) * tensor * psi)[0]

if __name__ == "__main__":
    with WolframCloudSession('C:\\Program Files\\Wolfram Research\\Wolfram Engine\\13.2\\WolframKernel.exe') as session:
        command = '2+2'
        x = session.evaluate(wlexpr(command))

    print(x)

    phi0 = Qubit(1, 0)
    phi1 = Qubit(1 / sqrt(3), sqrt(2 / 3))
    phi2 = Qubit(1 / sqrt(3), sp.E ** (sp.I * 2 * sp.pi / 3) * sqrt(2 / 3))
    phi3 = Qubit(1 / sqrt(3), sp.E ** (-sp.I * 2 * sp.pi / 3) * sqrt(2 / 3))

    M0 = simp(1 / 2) * phi0.state * phi0.dag
    M1 = simp(1 / 2) * phi1.state * phi1.dag
    M2 = simp(1 / 2) * phi2.state * phi2.dag
    M3 = simp(1 / 2) * phi3.state * phi3.dag

    tetra = [M0, M1, M2, M3]

    tetra_tensor = make_tensor(tetra)

    groups = search_all_combinations(tetra_tensor)
    # d = {"1": 1, "2": 2, "3": 3, "4": 4, "5": 5}
    # result = {}
    # with Pool(processes=4) as pool:
    #     for output in tqdm.tqdm(pool.imap_unordered(test, d.items()), total=len(d)):
    #         result[output[0]] = output[1]
    #
    # print(result)
    #
    # sp.mathematica_code()
