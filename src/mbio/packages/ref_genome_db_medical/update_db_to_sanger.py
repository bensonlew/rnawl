# -*- coding: utf-8 -*-
# __author__ : 'shicaiping'
# __date__: 20200317

import pymongo


def main(auto, genome_ids, test):
    tsg_client = pymongo.MongoClient("mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
    sanger_client = pymongo.MongoClient("mongodb://rna:y6a3t5n1w0y7@10.100.1.10/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
    tsg_db = tsg_client["sanger_ref_rna_v2"]
    sanger_db = sanger_client["sanger_ref_rna_v2"]
    tsg_db_collection = tsg_db['sg_genome_db']
    sanger_db_collection = sanger_db['sg_genome_db']
    if auto:
        ## find the latest genome_id
        genome_id_latest = "GM0001"
        sanger_docs = sanger_db_collection.find()
        for doc in sanger_docs:
            genome_id = doc['genome_id']
            if int(genome_id[2:]) > int(genome_id_latest[2:]):
                genome_id_latest = genome_id
        tsg_docs = tsg_db_collection.find()
        ## insert recently updated records
        insert_genome_id = dict()
        for tsg_doc in tsg_docs:
            tsg_genome_id = tsg_doc['genome_id']
            if int(tsg_genome_id[2:]) > int(genome_id_latest[2:]):
                del tsg_doc['_id']
                insert_genome_id[tsg_genome_id] = tsg_doc
        if insert_genome_id:
            for key in sorted(insert_genome_id):
                sanger_doc = sanger_db_collection.find_one({'genome_id': key})
                if not sanger_doc:
                    if test:
                        print key
                    else:
                        print ("Update genome_id is: {}".format(key))
                        sanger_db_collection.insert_one(insert_genome_id[key])
        else:
            print ("There is nothing to be updated.")
    if genome_ids:
        genome_ids = genome_ids.strip().split(",")
        for id in sorted(genome_ids):
            tsg_doc = tsg_db_collection.find_one({'genome_id': id})
            sanger_doc = sanger_db_collection.find_one({'genome_id': id})
            if not tsg_doc:
                print ("Wrong input genome_id: {}".format(id))
            if sanger_doc:
                print ("genome_id {} is already in sanger mongodb".format(id))
            if tsg_doc and not sanger_doc:
                del tsg_doc["_id"]
                if test:
                    print id
                else:
                    print ("Update genome_id is: {}".format(id))
                    sanger_db_collection.insert_one(tsg_doc)
            else:
                print ("There is nothing to be updated.")


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='This script is used to insert sg_genome_db info to sanger ref_rna_v2_db.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--auto', type=bool, default=False,
                        help='automatic update tsg sg_genome_db info to sanger sg_genome_db of refrna database, True of False, deault is False')
    group.add_argument('--genome_ids', type=str, default="",
                        help='genome id to be updated, if multiple should be separated by comma')
    parser.add_argument('--test', type=bool, default=False,
                        help='True print update info, False update sanger mongodb directly')

    args = parser.parse_args()

    main(args.auto, args.genome_ids, args.test)
