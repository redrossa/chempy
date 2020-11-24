from collections import Counter
from typing import Union, List

from chempy.util import tokenize


def _parse(toks: List[str]):
    op_stack: List[str] = []
    result: List[str] = []
    reserve: List[List[str]] = []

    def eval_op():
        if len(op_stack) == 0 or len(reserve) == 0:
            return

        op = op_stack.pop()

        if op == ')':
            op_stack.append(op)
        elif op.isnumeric():
            multiplier = int(op)
            reserve[-1] *= multiplier

    def merge_res():
        nonlocal reserve
        if len(reserve) > 1:
            merged = []
            for item in reserve:
                merged += item
            reserve = [merged]

    for tok in reversed(toks):
        if tok.isnumeric() or tok == ')':  # the valid operators
            op_stack.append(tok)
        elif tok == '(':
            while op_stack[-1] != ')':
                eval_op()
            op_stack.pop()  # pop matching ')'
            eval_op()
        elif tok.isalpha():
            reserve.append([tok])
            eval_op()
        else:
            raise ValueError('Invalid character while parsing molecule formula: ' + tok)

        merge_res()
        if len(op_stack) == 0:
            result += reserve.pop()

    return result


class Molecule:
    def __init__(self, formula_toks: List[str], charge: Union[str, int] = 0, states: List = None):
        self._formula = ''.join(formula_toks)
        self._elements = Counter(_parse(formula_toks))
        self._charge = int(charge)
        self._states = states if states else []
        self._charge_sign = '+' if self._charge >= 0 else '-'
        self._charge_mag = abs(self._charge)

    @staticmethod
    def complete_formula(complete_formula: str):
        toks: List[str] = tokenize(complete_formula)

        try:
            charge_start = toks.index('+')
        except ValueError:
            try:
                charge_start = toks.index('-')
            except ValueError:
                charge_start = 0  # charge cannot be first character

        if charge_start > 0:
            # sign, number, ... OR sign, ...
            try:
                charge_is_one = not toks[charge_start + 1].isnumeric()
            except IndexError:
                charge_is_one = True
            states_start = charge_start + (1 if charge_is_one else 2)
            charge = ''.join(toks[charge_start:states_start] + (['1'] if charge_is_one else []))
        else:
            charge = 0
            try:
                states_start = toks.index('(')
            except ValueError:
                states_start = 0  # states cannot be first character

        # states must take up the rest of the list
        states = toks[states_start + 1:len(toks) - 1] if states_start else None
        formula_toks = toks[:charge_start] if charge_start else toks[:states_start] if states_start else toks

        return Molecule(formula_toks, charge, states)

    @property
    def formula(self):
        return self._formula

    @property
    def elements(self):
        return self._elements

    @property
    def charge(self):
        return self._charge

    @property
    def states(self):
        return self._states.copy()

    @property
    def charge_sign(self):
        return self._charge_sign

    @property
    def charge_mag(self):
        return self._charge_mag

    def __hash__(self):
        return hash((self._formula, self._charge, *self._states))

    def __eq__(self, other: 'Molecule'):
        return isinstance(other, Molecule) \
               and self._formula == other._formula \
               and self._charge == other._charge \
               and self._states == other._states

    def __str__(self):
        return self._formula + \
               (self._charge_sign + (str(self._charge_mag) if self._charge_mag != 1 else '')
                if self._charge_mag != 0 else '') + \
               ('(' + ', '.join(state for state in self._states) + ')' if self._states else '')

    def __repr__(self):
        return self.__class__.__name__ + '(' \
               'formula=' + self._formula + ', ' \
               'charge='+ self._charge_sign + str(self._charge_mag) + ', ' + \
               'states=' + str(self._states) + \
               ')'
