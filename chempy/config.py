from pymongo import MongoClient

default_db = MongoClient(
    "mongodb+srv://user:E15XrvXBlTnVBpiM@chem-res.nrvkb.mongodb.net/chem-values-db?retryWrites=true&w=majority")\
    .get_default_database()

collections = {
    'elements': default_db['elements'],
    'thermodynamics': default_db['thermodynamics'],
    'bonds': default_db['bonds']
}
