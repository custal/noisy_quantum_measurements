import numpy as np
from dataclasses import dataclass
from sympy.matrices import Matrix
from sympy import Symbol, mathematica_code, sqrt
import sympy as sp
from pathlib import Path
from wolframclient.language import wlexpr
import inspect
import pickle
import base64

wolfram_kernel_path = Path('C:\\Program Files\\Wolfram Research\\Mathematica\\13.2\\WolframKernel.exe')
dir_data = Path(__file__).parents[0]/"data"

@dataclass
class ProtocolResult:
    """Class for keeping track of protocol results """
    protocol: frozenset
    measure_function: str
    measure_value: float

@dataclass
class EigenValueProtocolResult(ProtocolResult):
    eigen_values: tuple
    x_star_y_star: tuple

@dataclass
class POVMClone:
    """ Class for keeping track of povms once cloned and associated probabilities """
    povm: Matrix
    probability: Symbol

def encode_object(obj):
    return base64.b64encode(pickle.dumps(obj)).decode("ascii")

def decode_object(obj):
    return pickle.loads(base64.b64decode(obj))

def dataclass_to_dict(dclass):
    dclass = dclass.__dict__
    for key, item in dclass.items():
        if isinstance(item, Matrix):
            dclass[key] = matrix_to__list(item)
        elif isinstance(item, Symbol):
            dclass[key] = str(item)

    return dclass

def matrix_to__list(mat: Matrix):
    return mat.__array__().tolist()

def number_to_base(n, b, pack):
    """ n: Number in decimal
        b: Convert number to base b
        pack: Pack the number with zeros up to this value
    """
    if pack < n:
        raise ValueError(f"pack ({pack}) must be greater than n ({n})")

    num_of_digits = int(np.ceil(np.emath.logn(b, pack)))
    digits = [0]*num_of_digits
    i = num_of_digits-1
    while n:
        digits[i] = int(n % b)
        n //= b
        i -= 1
    return digits

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


#######################################################################################################################

def optimal_eig_new_measure(max_eig, min_eig):
    """ New measure of goodness paul sent to me in an email """
    if max_eig <= 1 - (np.sqrt(2)+1) * min_eig:
        (x, y) = (max_eig, (np.sqrt(2) - 1) * (1 - max_eig))
    elif 1 - (np.sqrt(2)+1) * min_eig <= max_eig and max_eig <= 1 - (np.sqrt(2)-1) * min_eig:
        (x, y) = (max_eig, min_eig)
    elif 1 - (np.sqrt(2)-1) * min_eig <= max_eig:
        (x, y) = (1 - (np.sqrt(2) - 1) * min_eig, min_eig)
    else:
        raise ValueError(f"Protocol failed.")

    return ((1 - x)**2 + y**2) / (1 - x + y), x, y, inspect.currentframe().f_code.co_name

def optimal_eig_old_measure(max_eig, min_eig):
    """ measure of goodness from the paper """
    if max_eig <= 1 - 2 * min_eig:
        (x, y) = (max_eig, (1 - max_eig) / 2)
    elif 1 - 2 * min_eig <= max_eig and max_eig <= 1 - min_eig / 2:
        (x, y) = (max_eig, min_eig)
    elif 1 - min_eig / 2 <= max_eig:
        (x, y) = ((1 - min_eig) / 2, min_eig)
    else:
        raise ValueError(f"Protocol failed.")

    return np.sqrt(1 / 3 * ((1 - x - y) ** 2 + (1 - x) * y)), x, y, inspect.currentframe().f_code.co_name

def eigen_value_measure_mathematica(protocol, povm_calculator, session):
    """ Calculate the eigen values of a protocol and use those to find the measure of goodness """
    permutations = povm_calculator.permutations
    dim = list(permutations.values())[0].povm.rows
    povm_sum = sp.zeros(dim)
    for permutation in protocol:
        povm_sum += permutations[permutation].povm

    command = 'Chop[N[Eigenvalues[' + mathematica_code(povm_sum) + ']]]'
    eigen_vals = session.evaluate(wlexpr(command))

    max_eig = max(eigen_vals)
    min_eig = min(eigen_vals)

    measure, x, y, func_name = optimal_eig_old_measure(max_eig, min_eig)

    return EigenValueProtocolResult(protocol, inspect.currentframe().f_code.co_name + "_" + func_name, measure,
                                    eigen_vals, (x, y))


def integral_rms_measure_mathematica(protocol, povm_calculator, session):
    """ Same as above but takes much longer and is not optimal as you have to evaluate the integral """
    a = povm_calculator.a
    b = povm_calculator.b
    permutations = povm_calculator.permutations

    probability = 0
    for permutation in protocol:
        probability += permutations[permutation].probability

    theta = sp.Symbol("θ", real=True)
    phi = sp.Symbol("𝜙", real=True)

    prob = probability - a * sp.conjugate(a)
    prob = prob.subs({a: sp.cos(theta / 2), b: sp.E ** (sp.I * phi) * sp.sin(theta / 2)})

    command = sqrt(sp.Integral(sp.sin(theta) / (4 * sp.pi) * prob ** 2, (theta, 0, sp.pi), (phi, 0, 2 * sp.pi)))
    command = 'N[ReleaseHold[' + mathematica_code(command) + ']]'
    measure = session.evaluate(wlexpr(command))

    return ProtocolResult(protocol, inspect.currentframe().f_code.co_name, measure)

def simplified_root_mean_square_measure(protocol, povm_calculator):
    """ Using the probabilities to calculate the goodness measure, used to check results are consistent between
    probability calculation of goodness measures above """
    a = povm_calculator.a
    b = povm_calculator.b
    permutations = povm_calculator.permutations

    probability = 0
    for permutation in protocol:
        probability += permutations[permutation].probability

    x = probability.subs({a: 1, b: 0}).evalf()
    y = probability.subs({a: 0, b: 1}).evalf()

    p_plus = probability.subs({a: 1 / sqrt(2), b: 1 / sqrt(2)}).evalf()
    p_minus = probability.subs({a: 1 / sqrt(2), b: -1 / sqrt(2)}).evalf()
    p_iplus = probability.subs({a: 1 / sqrt(2), b: sp.I / sqrt(2)}).evalf()
    p_iminus = probability.subs({a: 1 / sqrt(2), b: -sp.I / sqrt(2)}).evalf()
    z = 0.5 * ((p_plus - p_minus) - sp.I * (p_iplus - p_iminus)).evalf()

    measure = sqrt(1 / 3 * ((1 - x - y) ** 2 + (1 - x) * y + z * sp.conjugate(z)))
    measure = sp.simplify(measure)

    return ProtocolResult(protocol, inspect.currentframe().f_code.co_name, measure.evalf(chop=True))

def basic_measure(protocol, povm_calculator):
    """ Very simple measure initially used due to lack of understanding of the actual measure.
    Will need to set a and b sympy symbols to real if you want this to solve in any reasonable time scale """
    a = povm_calculator.a
    b = povm_calculator.b
    permutations = povm_calculator.permutations

    probability = 0
    for permutation in protocol:
        probability += permutations[permutation].probability

    prob = probability - a * sp.conjugate(a)
    prob = sp.simplify(prob.subs(b, sqrt(1 - a * sp.conjugate(a))).evalf())
    measure = sqrt(sp.Integral(prob ** 2, (a, 0, 1)))

    return ProtocolResult(protocol, inspect.currentframe().f_code.co_name, measure.evalf(chop=True))


#######################################################################################################################

def cyclic_symmetry(protocol):
    return type(protocol)([tuple([perm[i-1] for i in range(len(perm))]) for perm in protocol])
