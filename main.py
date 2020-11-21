from pymongo import MongoClient

import chempy

collection_name = 'thermodynamics'

# client = MongoClient("mongodb+srv://editor:IDskK7y6xTnE48qA@chem-res.nrvkb.mongodb.net/chem-values-db?retryWrites=true&w=majority")
# db = client.get_default_database()
# collection = db.get_collection(collection_name)
#
# element = ptable.element_symbol('H')

eq = chempy.equation("H2(g) + O2(g) = H2O(g) + H2(g)")
print(eq.balance())