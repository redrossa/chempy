# from typing import NamedTuple, List
#
# import chempy.config as _config
#
#
# class Element(NamedTuple):
#     protons: int
#     symbol: str
#     name: str
#     mass: str
#     config: str
#     radius: str
#     ionization: str
#     affinity: str
#     phase: str
#     mp: str
#     bp: str
#     density: str
#     series: str
#     elneg: str
#     oxidations: List[str]
#
#     def __str__(self):
#         return self.symbol
#
#     def __repr__(self):
#         return self.__str__()
#
#
# def element_protons(protons: int):
#     doc = _config.collections['elements'].find_one({'protons': protons})
#     return Element(*list(doc.values())[1:])
#
#
# def element_symbol(symbol: str):
#     doc = _config.collections['elements'].find_one({'symbol': symbol})
#     return Element(*list(doc.values())[1:])
#
#
# def element_name(name: str):
#     doc = _config.collections['elements'].find_one({'name': name.lower()})
#     return Element(*list(doc.values())[1:])
#
#
# def elements_with(**kwargs):
#     docs = _config.collections['elements'].find(dict(kwargs))
#     return [Element(*list(doc.values())[1:]) for doc in docs]
