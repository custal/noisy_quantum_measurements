import numpy as np
import sympy as sp
from sympy import Symbol
from sympy.matrices import Matrix
from sympy.physics.quantum import TensorProduct, Dagger
from itertools import combinations
from wolframclient.evaluation import WolframLanguageSession
from tqdm import tqdm
from datetime import datetime
from uuid import uuid1
from copy import deepcopy

from src.utils import number_to_base, POVMClone, wolfram_kernel_path, basic_measure, eigen_value_measure_mathematica,\
    dir_data, decode_object, encode_object


class POVMProtocolCalculator:
    a = Symbol("a")
    b = Symbol("b")
    time_format = '%Y_%m_%d-%H_%M_%S'

    def __init__(self, number_of_clones: int, *povms: Matrix):
        """ Set real to True if calculating using the basic measure """
        self.id_ = None
        self.date = None

        self.number_of_measures = len(povms)
        self.number_of_clones = number_of_clones
        self.number_of_permutations = self.number_of_measures**self.number_of_clones
        self.number_of_protocols = 2**self.number_of_permutations
        self.number_of_lowest_measures = None
        self.povms = povms
        self.calculate_permutations()
        self.calculate_all_protocols()
        self.lowest_measures = None

    def calculate_permutations(self):
        """ Calculate all permutations of the possible measurement outcomes for the given number of clones along with
        the associated povm operator and probabilities for a general state. """

        self.permutations = {}

        for i in range(self.number_of_permutations):
            permutation = tuple(number_to_base(i, self.number_of_measures, self.number_of_permutations))
            tensor_povm = []
            for j in permutation:
                tensor_povm.append(self.povms[j])

            tensor_povm = TensorProduct(*tensor_povm)
            probability = self.calculate_probability(tensor_povm)
            self.permutations[permutation] = POVMClone(tensor_povm, probability)

    def calculate_all_protocols(self):
        """ Find all possible combinations of the permutations """

        # We only need half because the final half is identitcal to the first half.
        # Eg. for 9 perms, taking groups of 2 and 7 is the same as taking groups of 7 and 2
        max_len = int(self.number_of_permutations / 2)

        permutations = list(self.permutations.keys())
        # Add empty set edge case
        self.protocols = {frozenset(): None, frozenset(permutations): None}

        for i in range(max_len):
            combs = combinations(permutations, i + 1)
            for selection in combs:
                duplicate_permutations = set(permutations)
                for outcome in selection:
                    duplicate_permutations.remove(outcome)

                self.protocols[frozenset(selection)] = None
                self.protocols[frozenset(duplicate_permutations)] = None

    def calculate_probability(self, povm):
        """ Find the probability of a generalised state. Using the column number of the povm is a bit hacky but should
        be fine"""
        psi = Matrix([0] * povm.cols)
        psi[0] = self.a
        psi[-1] = self.b

        return (Dagger(psi) * povm * psi)[0]

    def calculate_measures_sympy(self, measure_func=basic_measure):
        print(f"Calculating {self.number_of_protocols} measures...")

        for key, item in tqdm(self.protocols.items()):
            measure = measure_func(key, self)
            self.protocols[key] = measure

        print("Done!")

    def calculate_measures_mathematica(self, measure_func=eigen_value_measure_mathematica):
        print(f"Calculating {self.number_of_protocols} measures...")

        with WolframLanguageSession(wolfram_kernel_path) as session:
            for key, item in tqdm(self.protocols.items()):
                measure = measure_func(key, self, session)
                self.protocols[key] = measure

        print("Done!")

    def calculate_lowest_measures(self):
        lowest_val = np.inf
        for key, item in self.protocols.items():
            measure_value = item.measure_value
            try:
                # More sympy skulduggery. I think I should convert to doing all calculation within mathematica as
                # sympy just isn't sufficient
                check = bool(measure_value <= lowest_val)
            except TypeError:
                measure_value = sp.N(measure_value)

            if measure_value <= lowest_val:
                if measure_value == lowest_val:
                    lowest_measures.append(item)
                else:
                    lowest_measures = [item]
                lowest_val = measure_value

        self.lowest_measures = lowest_measures
        self.number_of_lowest_measures = len(lowest_measures)

    def apply_symmetry(self, symmetry):
        protocols = [prot.protocol for prot in self.lowest_measures]
        symmetry_groups = []

        for prot in protocols:
            found = False
            for group in symmetry_groups:
                if symmetry(prot) in group:
                    group.append(prot)
                    found = True
                    break
            if not found:
                symmetry_groups.append([prot])

        return symmetry_groups

    def save(self, file_name: str = None, dir: str = dir_data, mode: str = "a"):
        """ Save the data to a file. We don't pickle the entire object as it isn't
            conducive to long term storage ie. if we ever want to update this class. JSON is better for this.
        """
        # Generate a new id and date for each save as any changes made should be stored as a unique circuit
        self.id_ = str(uuid1())
        self.date = datetime.now()

        if file_name is None:
            file_name = f"{self.date.strftime(self.time_format)}-{self.id_}.txt"

        save_dict = deepcopy(self.__dict__)

        with open(dir/file_name, mode) as f:
            f.write(str(save_dict))
            f.write("\n")
            f.write(str(encode_object(self)))
            f.write("\n")

    @classmethod
    def load(cls, file_name: str, dir: str = dir_data):
        """
        Create a RandomCircuit object from a json string generated by the save method"""
        with open(dir/file_name, "r") as f:
            povm_pickle = f.readlines()[-1]

        return decode_object(povm_pickle)

if __name__ == "__main__":
    from utils import eigen_value_measure_new_mathematica

    # phi0 = Matrix([1, 0])
    # phi1 = Matrix([sp.nsimplify(1/2), sp.nsimplify(sp.sqrt(3)/2)])
    # phi2 = Matrix([sp.nsimplify(1/2), sp.nsimplify(-sp.sqrt(3)/2)])
    #
    # M0 = sp.nsimplify(2 / 3) * phi0 * Dagger(phi0)
    # M1 = sp.nsimplify(2 / 3) * phi1 * Dagger(phi1)
    # M2 = sp.nsimplify(2 / 3) * phi2 * Dagger(phi2)
    #
    # povm = POVMProtocolCalculator(2, M0, M1, M2)
    #
    # povm.calculate_measures_mathematica(eigen_value_measure_new_mathematica)
    #
    # povm.save("test.txt")
    # povm = POVMProtocolCalculator.load("test.txt")

    # phi0 = Matrix([1, 0])
    # phi1 = Matrix([sp.nsimplify(1 / sp.sqrt(3)), sp.nsimplify(sp.sqrt(2 / 3))])
    # phi2 = Matrix([sp.nsimplify(1 / sp.sqrt(3)), sp.E ** (sp.I * 2 * sp.pi / 3) * sp.nsimplify(sp.sqrt(2 / 3))])
    # phi3 = Matrix([sp.nsimplify(1 / sp.sqrt(3)), sp.E ** (-sp.I * 2 * sp.pi / 3) * sp.nsimplify(sp.sqrt(2 / 3))])
    #
    # M0 = sp.nsimplify(1 / 2) * phi0 * Dagger(phi0)
    # M1 = sp.nsimplify(1 / 2) * phi1 * Dagger(phi1)
    # M2 = sp.nsimplify(1 / 2) * phi2 * Dagger(phi2)
    # M3 = sp.nsimplify(1 / 2) * phi3 * Dagger(phi3)
    #
    # povm_tetra = POVMProtocolCalculator(2, M0, M1, M2, M3)
    #
    # povm_tetra.calculate_measures_mathematica()

    povm_tetra = POVMProtocolCalculator.load("2023_07_24-15_24_43-d577f36e-2a2d-11ee-954c-a4bb6d9bb15d.txt")

    povm_tetra.calculate_lowest_measures()

    print("")