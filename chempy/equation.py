import re
from typing import NamedTuple, List, Tuple

from sympy import eye
from sympy import Rational
from sympy import Matrix
from chempy.molecule import molecule
from chempy.molecule import Molecule

phases_p = re.compile(r'\((s|l|g|aq)[^)]*\)')


class Equation(NamedTuple):
    reactants: List[Tuple[float, Molecule]]
    products: List[Tuple[float, Molecule]]

    def __str__(self):
        return " + ".join([(str(int(reactant[0])) if reactant[0] != 1 else '')
                           + str(reactant[1]) for reactant in self.reactants]) + \
               " = " + \
               " + ".join([(str(int(product[0])) if product[0] != 1 else '')
                           + str(product[1]) for product in self.products])

    def __add__(self, other: 'Equation'):
        return balance(self.reactants + other.products, self.products + other.reactants)

    def __sub__(self, other: 'Equation'):
        pass

    def __mul__(self, scalar):
        pass

    def __rmul__(self, scalar):
        pass

    def balance(self):
        return balance(self.reactants, self.products)


def reverse(reactants: List[Tuple[float, Molecule]], products: List[Tuple[float, Molecule]]):
    return balance(products, reactants)


def simplify(reactants: List[Tuple[float, Molecule]], products: List[Tuple[float, Molecule]]):
    pass


def balance(reactants: List[Tuple[float, Molecule]], products: List[Tuple[float, Molecule]]):
    species = reactants + products
    elements = list(set([e for elements in [[*sp[1].elements] for sp in species] for e in elements]))
    cols = [[sp[0] * sp[1].elements[elements[i]] for i in range(len(elements))] for sp in species]

    m = Matrix(cols)
    if m.rows == m.cols:
        raise TypeError('The equation is an impossible reaction')

    identity = eye(m.rows)
    for i in range(identity.cols):
        m = m.col_insert(m.cols + i, identity.col(i))

    # Tuple returned from Matrix.rref()
    m = m.rref()[0]

    coeffs = [abs(Rational(n)) for n in list(m.row(m.rows - 1)) if n]
    denoms = [sol.q for sol in coeffs]
    max_den = max(denoms)
    coeffs = [float(sol * max_den) for sol in coeffs]

    # TODO: molecules on both sides of equation should cancel

    reactants = [(coeff, reactant[1]) for coeff, reactant in zip(coeffs[:len(reactants)], reactants)]
    products = [(coeff, product[1]) for coeff, product in zip(coeffs[len(reactants):], products)]

    return Equation(reactants, products)


def parse_expression(exp: str):
    exp = "".join(exp.split())
    species_seps = [match.end() for match in re.finditer(phases_p, exp)]
    entities = [exp[0:species_seps[0]]]

    for behind, sep in zip(species_seps[0:len(species_seps) - 1], species_seps[1:]):
        entities.append(exp[behind + 1:sep])

    species = []

    for entity in entities:
        coeff_end = 0
        while entity[coeff_end].isnumeric():
            coeff_end += 1
        coeff = entity[:coeff_end]
        coeff = 1.0 if not coeff_end else float(coeff)
        mol = molecule(entity[coeff_end:])
        species.append((coeff, mol))

    return species


def equation(eq: str):
    eq_exps = eq.split('=')
    if len(eq_exps) != 2:
        raise ValueError('Invalid chemical equation')

    reactants_exp = eq_exps[0]
    products_exp = eq_exps[1]

    reactants = parse_expression(reactants_exp)
    products = parse_expression(products_exp)

    return Equation(reactants, products)



