from itertools import tee, islice, zip_longest
from typing import List, Tuple

from core.molecule import Molecule
from core.util import tokenize

EquationSide = List[Tuple[float, Molecule]]


class Equation:
    def __init__(self, reactants: EquationSide, products: EquationSide):
        pass

    @staticmethod
    def from_str(eq: str):
        toks: List[str] = tokenize(eq)
        reserve: List[str] = []
        reactants: EquationSide = []
        products: EquationSide = []

        def get_next(it, window=1):
            items, nexts = tee(it, 2)
            nexts = islice(nexts, window, None)
            return zip_longest(items, nexts)

        def merge_res(start: int = 0):
            nonlocal reserve
            if len(reserve) > 1:
                merged = ''.join(reserve[start:])
                reserve = reserve[0:start] + [merged]

        ls = reactants
        for tok, ahead in get_next(toks):
            if (tok == '+' or tok == '=') and (ahead.isnumeric() or ahead.isalpha()):
                if tok == '=':
                    ls = products
                merge_res(1 if reserve[0].isnumeric() else 0)
                ls.append((1.0 if len(reserve) == 1 else float(reserve[0]),
                                  Molecule.complete_formula(reserve[1])))
                reserve = []
            else:
                reserve += tok

        return reactants, products


# print(Molecule.complete_formula('H2[g]'))
print(Equation.from_str('2H2[g] + O2[g] = 2H2O[g]'))