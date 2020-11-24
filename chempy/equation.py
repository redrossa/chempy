import math
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import tee, islice, zip_longest
from typing import List, Union, Dict

from sympy import Matrix

from chempy.molecule import Molecule
from chempy.util import tokenize


def _parse_coeff(src: str = ''):
    if not src:
        frac = Fraction(1)
    elif src.isnumeric():
        frac = Fraction(int(src))
    elif '.' in src:
        ratio = float(src).as_integer_ratio()
        frac = Fraction(ratio[0], ratio[1])
    elif '/' in src:
        ratio = src.split('/')
        if len(ratio) != 2:
            raise ValueError('Invalid coefficient numerical value')
        frac = Fraction(int(ratio[0]), int(ratio[1]))
    else:
        raise ValueError('Invalid coefficient numerical value')
    return frac


class Species:
    def __init__(self, molecule: Molecule, coeff: Union[int, float, str, Fraction] = Fraction(1)):
        num, denom = float(coeff).as_integer_ratio() if isinstance(coeff, float) else (0, 1)
        self._coeff = coeff if isinstance(coeff, Fraction) \
            else _parse_coeff(coeff) if isinstance(coeff, str) \
            else Fraction(num, denom) if isinstance(coeff, float) \
            else Fraction(coeff)
        self._molecule = molecule
        self._atoms = self._molecule.elements.copy()
        for key in self._atoms.keys():
            self._atoms[key] *= self._coeff

    @property
    def coeff(self):
        return self._coeff

    @property
    def molecule(self):
        return self._molecule

    @property
    def atoms(self):
        return self._atoms

    def as_dict(self):
        return {self._molecule: self._coeff}

    def __mul__(self, other):
        return Species(self._molecule, self._coeff * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return Species(self._molecule, self._coeff / other)

    def __rtruediv__(self, other):
        return self.__truediv__(other)

    def __str__(self):
        return (str(self._coeff) if self._coeff.denominator > 1
                else str(int(self._coeff)) if self._coeff.numerator != 1
                else '') + str(self._molecule)

    def __repr__(self):
        return self.__class__.__name__ + '(coeff=' + str(self._coeff) + ', molecule=' + repr(self._molecule) + ')'


class Equation:
    def __init__(self, reactants: List[Species], products: List[Species]):
        self._reactants = reactants
        self._products = products

    @staticmethod
    def from_str(eq: str):
        toks: List[str] = tokenize(eq)
        reactants: List[Species] = []
        products: List[Species] = []

        def get_next(it, window=1):
            items, nexts = tee(it, 2)
            nexts = islice(nexts, window, None)
            return zip_longest(items, nexts)

        def eval_species(start: int, end: int):
            mol_start = start
            while not toks[mol_start].isalpha() and mol_start < end:
                mol_start += 1
            if mol_start == end:
                raise ValueError('Invalid equation')

            coeff = ''.join(toks[start:mol_start])
            mol_toks = toks[mol_start:end]
            mol = Molecule.complete_formula(''.join(mol_toks))
            return Species(mol, coeff=coeff)

        ls = reactants
        prev_plus_index = -1
        for i, items in enumerate(get_next(toks)):
            tok, ahead = items
            try:
                if tok == '=' or tok == '+' and (ahead.isnumeric() or ahead.isalpha()):
                    ls += [eval_species(prev_plus_index + 1, i)]
                    prev_plus_index = i
                    if tok == '=':  # switch to products if '=' is encountered
                        if ls is reactants:
                            ls = products
                        else:
                            raise ValueError('Invalid equation')
            except AttributeError:  # ahead == NoneType has no isnumeric or isalpha
                break
        ls += [eval_species(prev_plus_index + 1, len(toks))]

        return Equation(reactants, products)

    @staticmethod
    def from_dict(reactants: Dict[Molecule, Fraction], products: Dict[Molecule, Fraction]):
        rs = [Species(r, reactants[r]) for r in reactants]
        ps = [Species(p, products[p]) for p in products]
        return Equation(rs, ps)

    def simplify(self):
        species = self._reactants + self._products

        # Turn all fractions into whole numbers
        max_denom = max([sp.coeff.denominator for sp in species])
        if max_denom > 1:
            species = [sp * max_denom for sp in species]

        # Reduce all coefficients to lowest terms
        gcd = reduce(math.gcd, [sp.coeff.numerator for sp in species])
        if gcd > 1:
            species = [sp / gcd for sp in species]

        # Reduce species found on both sides
        reactants, products = Counter(), Counter()
        for reactant in species[:len(self._reactants)]:
            reactants.update(reactant.as_dict())
        for product in species[len(self._reactants):]:
            products.update(product.as_dict())
        union = reactants & products
        for mol in union:
            reactants[mol] -= union[mol]
            products[mol] -= union[mol]
        reactants = {r: reactants[r] for r in reactants if reactants[r] > 0}
        products = {p: products[p] for p in products if products[p] > 0}

        return Equation.from_dict(reactants, products)

    def balance(self):
        simplified = self.simplify()
        species = simplified._reactants + simplified._products

        # Get the number of elements in each species and all elements present in Equation to produce the matrix's rows
        elm_counters = [sp.molecule.elements for sp in species]
        elements = {}
        for counter in elm_counters:
            elms = counter.keys()
            elements |= elms
        elements = dict.fromkeys(elements, Fraction())

        # Use sympy.Matrix to find null space of matrix whose values are to be used as balanced coefficients
        mat_list = [list({**elements, **counter}.values()) for counter in elm_counters]
        eye = Matrix.eye(len(mat_list))
        mat_list = [row + eye.row(i).tolist()[0] for i, row in enumerate(mat_list)]
        mat = Matrix(mat_list).rref()[0]
        sols = [Fraction(abs(coeff.numerator()), coeff.denominator())
                for coeff in mat.row(len(mat_list) - 1).tolist()[0] if abs(coeff) > 0]

        species = [Species(sp.molecule, coeff=sol) for sp, sol in zip(species, sols)]
        return Equation(species[:len(simplified._reactants)], species[len(simplified._reactants):])

    def __add__(self, other: 'Equation'):
        if not isinstance(other, Equation):
            raise TypeError("Expecting type 'Equation', got '" + str(type(other).__name__) + "' instead")
        return Equation(self._reactants + other._reactants, self._products + other._products)

    def __sub__(self, other: 'Equation'):
        if not isinstance(other, Equation):
            raise TypeError("Expecting type 'Equation', got '" + str(type(other).__name__) + "' instead")
        return Equation(self._reactants + other._products, self._products + other._reactants)

    def __mul__(self, other: Union[int, float, str, Fraction]):
        if not isinstance(other, (int, float, str, Fraction)):
            raise TypeError("Expecting type 'int', 'float', 'str', or 'Fraction', got '" +
                            str(type(other).__name__) + "' instead")

        num, denom = float(other).as_integer_ratio() if isinstance(other, float) else (0, 1)
        scalar = other if isinstance(other, Fraction) \
            else _parse_coeff(other) if isinstance(other, str) \
            else Fraction(num, denom) if isinstance(other, float) \
            else Fraction(other)

        reactants = [scalar * r for r in self._reactants]
        products = [scalar * p for p in self._products]
        return Equation(reactants, products)

    def __rmul__(self, other: Union[int, float, str, Fraction]):
        return self.__mul__(other)

    def __truediv__(self, other: Union[int, float, str, Fraction]):
        if not isinstance(other, (int, float, str, Fraction)):
            raise TypeError("Expecting type 'int', 'float', 'str', or 'Fraction', got '" +
                            str(type(other).__name__) + "' instead")

        num, denom = float(other).as_integer_ratio() if isinstance(other, float) else (0, 1)
        scalar = other if isinstance(other, Fraction) \
            else _parse_coeff(other) if isinstance(other, str) \
            else Fraction(num, denom) if isinstance(other, float) \
            else Fraction(other)

        reactants = [scalar / r for r in self._reactants]
        products = [scalar / p for p in self._products]
        return Equation(reactants, products)

    def __rtruediv__(self, other: Union[int, float]):
        return self.__truediv__(other)
    
    def __neg__(self):
        return Equation(self._products, self._reactants)
    
    def __pos__(self):
        return Equation(self._reactants, self._products)
    
    @property
    def reactants(self):
        return self._reactants

    @property
    def products(self):
        return self._products

    def __str__(self):
        return ' + '.join([str(species) for species in self._reactants]) + \
               ' = ' + \
               ' + '.join([str(species) for species in self._products])
