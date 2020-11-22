import collections
from typing import Union, List


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


from core import tokenize


print(_parse(tokenize('AliefOfc3')))


class Molecule:
    __slots__ = ('_formula', '_elements', '_charge', '_states')

    def __init__(self, formula_toks: List[str], charge: Union[str, int] = 0, states: List = None):
        self._formula = ''.join(formula_toks)
        self._elements = collections.Counter(_parse(formula_toks))
        self._charge = charge
        self._states = states if not states else []

    @staticmethod
    def complete_formula(complete_formula: str):
        pass

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
        return self._states
