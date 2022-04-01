import json
import logging
import re
import sys

import regex
from biocluster.config import Config

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)


from mbio.packages.ref_rna_v2.annotation.go_annotation2 import Terms




def update_go(GO_db, version=1.0):
    client = Config().get_mongo_client(mtype="ref_rna", ref=True)
    ref_db = client[Config().get_mongo_dbname("ref_rna", ref=True)]
    go = ref_db.GO
    cursor = go.find({})
    go_term_dict = {document['go_id']: document for document in cursor}
    for go_id, term in GO_db.items():
        if go_id in go_term_dict:
            pass
        else:
            if hasattr(term, "syn"):
                syn_list = [x.split("\"")[1] for x in term.syn]
            else:
                syn_list = []

            if hasattr(term, "is_a"):
                is_a_list = [x[0].id for x in term.is_a]
            else:
                is_a_list = []

            term_dict = {
                "is_a": is_a_list,
                "definition": term.defi.split("\"")[1],
                "depth": max(term.levels) if len(term.levels) !=0 else 0,
                "level": min(term.levels)  if len(term.levels) !=0 else 0,
                "version": version,
                "name": term.name,
                "synonym":  syn_list,
                "ontology": term.namespace,
                "go_id": go_id
            }

            print term_dict
            go.insert_one(term_dict)



if __name__ == '__main__':
    GO_db = Terms(obo_fp=sys.argv[2])
    version = float(sys.argv[1])
    update_go(GO_db, version)
