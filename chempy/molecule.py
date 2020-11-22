import collections
from typing import List, Dict, NamedTuple

valid_elements = [
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]


class Molecule(NamedTuple):
    formula: str
    elements: Dict[str, float]
    charge: str
    phases: List[str]

    def __str__(self):
        return self.formula

    def __repr__(self):
        return self.__str__()


def tokenize(formula: str) -> List[str]:
    if '!' in formula:
        raise ValueError('Invalid token in "' + formula + "' at index " + str(formula.find('!')))

    tokens = []
    store = ""

    def add_tok(tok: str):
        if tok != '':
            tokens.append(tok)

    def retrieve_store():
        nonlocal store
        tmp = store
        store = ''
        return tmp

    for i, c in enumerate(formula + '!'):
        if c == '(' or c == ')':
            add_tok(retrieve_store())
            add_tok(c)
        elif c.isupper():
            add_tok(retrieve_store())
            store += c
        elif c.islower():
            add_tok(retrieve_store() + c)
        elif c.isnumeric():
            if not store.isnumeric():
                add_tok(retrieve_store())
            store += c
        elif c == '!':
            add_tok(retrieve_store())
        else:
            raise ValueError('Invalid token in "' + formula + "' at index " + str(i))
    return tokens


def parse_molecule(formula: str):
    op_stack = []
    result = []
    reserve = []
    tokens = tokenize(formula)

    def eval_op():
        if len(op_stack) == 0 or len(reserve) == 0:
            return

        op = op_stack.pop()

        if op == ')':
            op_stack.append(op)
        elif op.isnumeric():
            multiplier = int(op)
            reserve[-1] *= multiplier
        else:
            raise ValueError('Invalid operator "' + op + '"')

    def merge_res():
        nonlocal reserve
        if len(reserve) == 2:
            reserve = [reserve[0] + reserve[1]]

    for tok in reversed(tokens):
        if tok.isnumeric() or tok == ')':
            op_stack.append(tok)
        elif tok == '(':
            while op_stack[-1] != ')':
                eval_op()
            op_stack.pop()  # pop matching ')'
            eval_op()
        elif tok in valid_elements:
            reserve.append([tok])
            eval_op()
        else:
            raise ValueError('Invalid element "' + tok + '" in "' + formula + '"')

        merge_res()
        if len(op_stack) == 0:
            result += reserve.pop()

    return result


def molecule(formula: str):
    formula = "".join(formula.split())

    if formula[-1] != ')':
        raise ValueError('Molecule phases not listed')

    phase_index = formula.rfind('(')
    charge_index = formula.rfind('+')
    if charge_index == -1:
        charge_index = formula.rfind('-')

    phase = formula[phase_index + 1:len(formula) - 1].split(',')
    charge = "0" if charge_index == -1 else formula[charge_index:phase_index]
    elements_formula = formula[0:phase_index if charge_index == -1 else charge_index]

    parsed_formula = parse_molecule(elements_formula)

    return Molecule(formula, collections.Counter(parsed_formula), charge, phase)
