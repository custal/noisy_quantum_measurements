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
import json
from copy import deepcopy

from utils import number_to_base, POVMClone, wolfram_kernel_path, basic_measure, eigen_value_measure_mathematica,\
    dir_data


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

        self.calculate_lowest_measures()
        print("Done!")

    def calculate_measures_mathematica(self, measure_func=eigen_value_measure_mathematica):
        print(f"Calculating {self.number_of_protocols} measures...")

        with WolframLanguageSession(wolfram_kernel_path) as session:
            for key, item in tqdm(self.protocols.items()):
                measure = measure_func(key, self, session)
                self.protocols[key] = measure

        self.calculate_lowest_measures()
        print("Done!")

    def calculate_lowest_measures(self):
        lowest_val = np.inf
        for key, item in self.protocols.items():
            if item.measure_value <= lowest_val:
                if item.measure_value == lowest_val:
                    lowest_measures.append(item)
                else:
                    lowest_measures = [item]
                lowest_val = item.measure_value

        self.lowest_measures = lowest_measures
        self.number_of_lowest_measures = len(lowest_measures)

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

    @classmethod
    def load(cls, json_string: str):
        """
        WORK IN PROGRESS
        Create a RandomCircuit object from a json string generated by the save method"""
        load_dict = json.loads(json_string)
        # Decode objects which couldn't be JSONified
        load_dict["date"] = datetime.strptime(load_dict["date"], POVMProtocolCalculator.time_format)
        load_dict["state"] = list(map(complex, load_dict["state"]))
        load_dict["circuit"] = load_dict["circuit"]

        return cls(**load_dict)


if __name__ == "__main__":
    phi0 = Matrix([1, 0])
    phi1 = Matrix([sp.nsimplify(1/2), sp.nsimplify(sp.sqrt(3)/2)])
    phi2 = Matrix([sp.nsimplify(1/2), sp.nsimplify(-sp.sqrt(3)/2)])

    M0 = sp.nsimplify(2 / 3) * phi0 * Dagger(phi0)
    M1 = sp.nsimplify(2 / 3) * phi1 * Dagger(phi1)
    M2 = sp.nsimplify(2 / 3) * phi2 * Dagger(phi2)

    povm = POVMProtocolCalculator(2, M0, M1, M2)

    povm.calculate_measures_sympy()

    povm.save("test.txt")