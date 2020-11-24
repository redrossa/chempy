import unittest
from fractions import Fraction

from chempy import equation, Molecule


class CoefficientParseTest(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(Fraction(1), equation._parse_coeff(''))
        self.assertEqual(Fraction(1), equation._parse_coeff('1'))
        self.assertEqual(Fraction(12), equation._parse_coeff('12'))
        self.assertEqual(Fraction(3, 2), equation._parse_coeff('   1  .   5 '))
        self.assertEqual(Fraction(1, 5), equation._parse_coeff('1/5'))
        self.assertEqual(Fraction(1, 5), equation._parse_coeff('   1  /   5 '))
        self.assertRaises(ValueError, lambda: equation._parse_coeff('-2'))


molecules = [
    Molecule.complete_formula('e-'),
    Molecule.complete_formula('He(g_'),
    Molecule.complete_formula('H2O(l)'),
    Molecule.complete_formula('SO4-2(aq)'),
    Molecule.complete_formula('HgS(s, red)'),
    Molecule.complete_formula('Fe2O3(s)'),
    Molecule.complete_formula('FePO4(s)'),
    Molecule.complete_formula('Ca(NO3)2(s)'),
    Molecule.complete_formula('(NH4)2SO4(s)'),
    Molecule.complete_formula('(CH2)2(NH2)2H+(aq)')
]

coefficients = [
    [0, 0.0, '0', '0.0', '0/1', Fraction()],
    [1, 1.0, '1', '1.0', '1/1', Fraction(1)],
    [3, 3.0, '3', '3.0', '3/1', Fraction(3)],
    [   1.5,      '1.5', '3/2', Fraction(3, 2)]
]


class SpeciesTest(unittest.TestCase):
    def test_init(self):
        # default coeff == 1
        molecule = molecules[0]
        sp = equation.Species(molecule)
        self.assertEqual(sp.coeff, Fraction(1))
        self.assertEqual(sp.molecule, molecule)
        self.assertEqual(sp.atoms, molecule.elements)

        # passing coeff
        for coeffs in coefficients:
            for coeff, molecule in zip(coeffs, molecules):
                sp = equation.Species(molecule, coeff)
                self.assertEqual(sp.coeff, coeffs[len(coeffs) - 1])
                self.assertEqual(sp.molecule, molecule)
                atoms = molecule.elements
                for key in atoms:
                    atoms[key] *= coeffs[len(coeffs) - 1]
                self.assertEqual(sp.atoms, atoms)

    def test_as_dict(self):
        for coeffs in coefficients:
            for coeff, molecule in zip(coeffs, molecules):
                sp = equation.Species(molecule, coeff)
                self.assertEqual(sp.as_dict(), {molecule: coeffs[len(coeffs) - 1]})

    def test_mul(self):
        for coeffs in coefficients:
            for coeff, molecule in zip(coeffs, molecules):
                sp = equation.Species(molecule) * coeffs[len(coeffs) - 1]
                self.assertEqual(sp.coeff, coeffs[len(coeffs) - 1])
                self.assertEqual(sp.molecule, molecule)
                atoms = molecule.elements
                for key in atoms:
                    atoms[key] *= coeffs[len(coeffs) - 1]
                self.assertEqual(sp.atoms, atoms)

    def test_div(self):
        for coeffs in coefficients[1:]:
            for coeff, molecule in zip(coeffs, molecules):
                sp = equation.Species(molecule, coeff) / coeffs[len(coeffs) - 1]
                self.assertEqual(sp.coeff, Fraction(1))
                self.assertEqual(sp.molecule, molecule)
                self.assertEqual(sp.atoms, molecule.elements)


class EquationTests(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
