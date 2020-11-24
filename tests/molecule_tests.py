import unittest
from collections import Counter

from chempy import Molecule
from chempy import molecule
from chempy import util


case_inputs = [
    ('e', '-', []),
    ('He', 0, ['g']),
    ('H2O', '0', ['l']),
    ('SO4', -2, ['aq']),
    ('HgS', 0, ['s', 'red']),
    ('Fe2O3', 0, ['s']),
    ('FePO4', 0, ['s']),
    ('Ca(NO3)2', 0, ['s']),
    ('(NH4)2SO4', 0, ['s']),
    ('(CH2)2(NH2)2H', '+1', ['aq'])
]


def adjust_charge(charge):
    return charge if isinstance(charge, int) else charge + ('1' if charge == '+' or charge == '-' else '')

def molecule_as_str(formula, charge, states):
    charge = str(charge) if str(charge) != '0' else ''
    tmp = formula + charge + '(' + ', '.join(states) + ')'
    return tmp


class MoleculeParseTests(unittest.TestCase):
    def test_parse(self):
        cases = [
            (molecule._parse(util.tokenize(case_inputs[0][0])), ['e']),
            (molecule._parse(util.tokenize(case_inputs[1][0])), ['He']),
            (molecule._parse(util.tokenize(case_inputs[2][0])), ['O', 'H', 'H']),
            (molecule._parse(util.tokenize(case_inputs[3][0])), ['O', 'O', 'O', 'O', 'S']),
            (molecule._parse(util.tokenize(case_inputs[4][0])), ['S', 'Hg']),
            (molecule._parse(util.tokenize(case_inputs[5][0])), ['O', 'O', 'O', 'Fe', 'Fe']),
            (molecule._parse(util.tokenize(case_inputs[6][0])), ['O', 'O', 'O', 'O', 'P', 'Fe']),
            (molecule._parse(util.tokenize(case_inputs[7][0])), ['O', 'O', 'O', 'N', 'O', 'O', 'O', 'N', 'Ca']),
            (molecule._parse(util.tokenize(case_inputs[8][0])), ['O', 'O', 'O', 'O', 'S', 'H', 'H', 'H', 'H', 'N', 'H', 'H', 'H', 'H', 'N']),
            (molecule._parse(util.tokenize(case_inputs[9][0])), ['H', 'H', 'H', 'N', 'H', 'H', 'N', 'H', 'H', 'C', 'H', 'H', 'C'])
        ]

        for actual, expected in cases:
            self.assertEqual(expected, actual)


class MoleculeTests(unittest.TestCase):
    def test_Molecule_init(self):
        cases = [(Molecule(util.tokenize(formula), adjust_charge(charge), states),
                  formula,
                  Counter(molecule._parse(util.tokenize(formula))),
                  int(adjust_charge(charge)),
                  states)
                 for formula, charge, states in case_inputs]

        for actual, exp_formula, exp_elements, exp_charge, exp_states in cases:
            self.assertEqual(exp_formula, actual.formula)
            self.assertEqual(exp_elements, actual.elements)
            self.assertEqual(exp_charge, actual.charge)
            self.assertEqual(exp_states, actual.states)

    def test_Molecule_complete_formula(self):
        cases = [(Molecule.complete_formula(molecule_as_str(formula, charge, states)),
                  formula,
                  Counter(molecule._parse(util.tokenize(formula))),
                  int(adjust_charge(charge)),
                  states)
                 for formula, charge, states in case_inputs]

        for actual, exp_formula, exp_elements, exp_charge, exp_states in cases:
            self.assertEqual(exp_formula, actual.formula)
            self.assertEqual(exp_elements, actual.elements)
            self.assertEqual(exp_charge, actual.charge)
            self.assertEqual(exp_states, actual.states)


if __name__ == '__main__':
    unittest.main()
