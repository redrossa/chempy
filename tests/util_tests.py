import unittest

import chempy.util as util


class UtilTests(unittest.TestCase):
    def test_tokenize(self):
        self.assertEqual(util.tokenize('1'), ['1'])
        self.assertEqual(util.tokenize('12'), ['12'])
        self.assertEqual(util.tokenize('1/2'), ['1', '/', '2'])
        self.assertEqual(util.tokenize('1.2'), ['1', '.', '2'])
        self.assertEqual(util.tokenize('H'), ['H'])
        self.assertEqual(util.tokenize('He'), ['He'])
        self.assertEqual(util.tokenize('Helium'), ['Helium'])
        self.assertEqual(util.tokenize('1H'), ['1', 'H'])
        self.assertEqual(util.tokenize('12He'), ['12', 'He'])
        self.assertEqual(util.tokenize('12Helium'), ['12', 'Helium'])
        self.assertEqual(util.tokenize('1/2Helium'), ['1', '/', '2', 'Helium'])
        self.assertEqual(util.tokenize('1.2Helium'), ['1', '.', '2', 'Helium'])
        self.assertEqual(util.tokenize('1.2H+'), ['1', '.', '2', 'H', '+'])
        self.assertEqual(util.tokenize('1.2He+'), ['1', '.', '2', 'He', '+'])
        self.assertEqual(util.tokenize('1.2Helium+'), ['1', '.', '2', 'Helium', '+'])
        self.assertEqual(util.tokenize('1.2e-'), ['1', '.', '2', 'e', '-'])
        self.assertEqual(util.tokenize('1.2electron-'), ['1', '.', '2', 'electron', '-'])
        self.assertEqual(util.tokenize('1.2electron- + heat'), ['1', '.', '2', 'electron', '-', '+', 'heat'])
        self.assertEqual(util.tokenize('1.2camelCase'), ['1', '.', '2', 'camel', 'Case'])
        self.assertEqual(util.tokenize('1.2camelCase(s,g)'), ['1', '.', '2', 'camel', 'Case', '(', 's', ',', 'g', ')'])
        self.assertEqual(util.tokenize('1.2camCas+(s,g)'), ['1', '.', '2', 'cam', 'Cas', '+', '(', 's', ',', 'g', ')'])
        self.assertEqual(util.tokenize('1.2ionIon+(s)='), ['1', '.', '2', 'ion', 'Ion', '+', '(', 's', ')', '='])
        self.assertEqual(util.tokenize('2NH + Ofc3'), ['2', 'N', 'H', '+', 'Ofc', '3'])
        self.assertEqual(util.tokenize('2NH + Ofc3/2'), ['2', 'N', 'H', '+', 'Ofc', '3', '/', '2'])
        self.assertEqual(util.tokenize('23.3/3syms: ! # &^asd ( |][? $% = *% { ;Tok? $ ,)'),
                         ['23', '.', '3', '/', '3', 'syms', ':', '!', '#', '&', '^', 'asd', '(', '|', ']', '[', '?',
                          '$', '%', '=', '*', '%', '{', ';', 'Tok', '?', '$', ',', ')'])


if __name__ == '__main__':
    unittest.main()
