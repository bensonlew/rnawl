# -*- coding: utf-8 -*-
# __author__ = 'sanger'
from __future__ import division
import os
from biocluster.config import Config
from bson.objectid import ObjectId
import types
import json
import re
from types import StringTypes
import gridfs
from collections import OrderedDict
import pandas as pd
from itertools import chain

# project_type = 'labelfree'
# db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
# if project_type == 'labelfree' or project_type == 'denovo_rna_v2':
#     main_table = 'sg_express'
#     detail_table = 'sg_express_detail'
# else:
#     main_table = 'sg_exp'
#     detail_table = 'sg_exp_detail'

def export_exp_matrix_labelfree(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    conn = db['sg_express_detail']
    group_dict = bind_obj.sheet.option('protein_group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    con_main = db['sg_express']
    express_id = con_main.find_one({"task_id": data, "type": "ratio"})["main_id"]
    exp_records = conn.find({"express_id": express_id}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('accession_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    return output

def export_exp_matrix_denovo(data, option_name, dir_path, bind_obj=None):
    project_type = 'denovo_rna_v2'
    #add by fwy 20210128
    db = selecet_db(project_type, "sg_exp_detail", "exp_id",data, bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('rna_group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    exp_records = conn.find({"exp_id": ObjectId(data)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('seq_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_exp_matrix_refrnav1(data,option_name,dir_path,bind_obj=None):
    """
    type????????????gene/transcript?????????workflow????????????????????????
    express_level????????????fpkm/tpm?????????workflow????????????????????????
    """
    # db = Config().mongo_client[Config().MONGODB + "_ref_rna"]
    # add by fwy 20210128
    db = selecet_db("ref_rna", "sg_express_detail", 'express_id',data, bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
    fpkm_path = os.path.join(dir_path, "%s_fpkm.matrix" % option_name)
    count_path = os.path.join(dir_path, "%s_count.matrix" % option_name)
    bind_obj.logger.debug("????????????????????????:%s???fpkm??????:%s" % (count_path, fpkm_path))
    collection = db['sg_express_detail']
    my_collection = db['sg_express']
    type = "gene"
    bind_obj.logger.debug(type)
    level = bind_obj.sheet.option("express_level")
    group_dict = bind_obj.sheet.option('rna_group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    #sample_group = bind_obj.sheet.option("sample_group")
    results = collection.find({'$and': [{'express_id': ObjectId(data)}, {'type': '{}'.format(type)},{"sample_group":"sample"},{"value_type":level}]})
    # count_results = collection.find({'$and': [{'express_id': ObjectId(data)}, {'type': '{}'.format(type)},{"sample_group":"sample"},{"value_type":"count"}]})
    my_result = my_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????express_id:{}???sg_express???????????????".format(ObjectId(data)))
    # samples = my_result['specimen']
    def write_file(path, collcetion_results):
        with open(path, "wb") as f:
            head = '\t'.join(samples)
            f.write('\t' + head + '\n')
            for result in collcetion_results:
                #bind_obj.logger.debug(result)
                gene_id = result['seq_id']
                fpkm_write = '{}'.format(gene_id)
                for sam in samples:
                    fpkm = sam
                    try:
                        fpkm_write += '\t{}'.format(result[fpkm])
                        #count_write += '\t{}'.format(result[count])
                    except Exception:
                        pass
                        #bind_obj.set_error("{}??????".format(result[fpkm]))
                fpkm_write += '\n'
                #count_write += '\n'
                f.write(fpkm_write)
                #c.write(count_write)
    write_file(fpkm_path, results)
    # write_file(count_path, count_results)
    # paths = ','.join([fpkm_path, count_path])
    return fpkm_path

def export_exp_matrix_refrnav2(data, option_name, dir_path, bind_obj=None):
    project_type = 'ref_rna_v2'
    # add by fwy 20210128
    db = selecet_db(project_type, "sg_exp_detail","exp_id", data, bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('rna_group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    exp_records = conn.find({"exp_id": ObjectId(data)}, target_cols)
    if 'type' in bind_obj.sheet.options():
        print("parameter 'type' found, and we will decide ref or new ")
        new_exp_records = list()
        if bind_obj.sheet.option('type') == "ref":
            for record in exp_records:
                if not record['seq_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                    new_exp_records.append(record)
            exp_records = new_exp_records
        elif bind_obj.sheet.option('type') == 'new':
            for record in exp_records:
                if record['seq_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                    new_exp_records.append(record)
            exp_records = new_exp_records
        else:
            pass
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('seq_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_exp_matrix_medical_transcriptome(data, option_name, dir_path, bind_obj=None):
    project_type = 'medical_transcriptome'
    # add by fwy 20210128
    db = selecet_db(project_type, "sg_exp_detail","exp_id", data, bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('rna_group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(transcript_id=1,gene_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    exp_records = conn.find({"exp_id": ObjectId(data)}, target_cols)
    if 'type' in bind_obj.sheet.options():
        print("parameter 'type' found, and we will decide ref or new ")
        new_exp_records = list()
        if bind_obj.sheet.option('type') == "ref":
            for record in exp_records:
                try:
                    if not record['transcript_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                        new_exp_records.append(record)
                    good_id = "transcript_id"
                except:
                    if not record['gene_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                        new_exp_records.append(record)
                    good_id = "gene_id"
            exp_records = new_exp_records
        elif bind_obj.sheet.option('type') == 'new':
            for record in exp_records:
                try:
                    if record['transcript_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                        new_exp_records.append(record)
                    good_id = "transcript_id"
                except:
                    if record['gene_id'].startswith(('MSTRG', 'TCONS', 'XLOC')):
                        new_exp_records.append(record)
                    good_id = "gene_id"
            exp_records = new_exp_records
        else:
            pass
    exp_matrix = pd.DataFrame(list(exp_records))
    if "transcript_id" in list(exp_matrix.columns):
        good_id = "transcript_id"
    elif "gene_id" in list(exp_matrix.columns):
        good_id = "gene_id"
    else:
        bind_obj.set_error("????????????????????????{}????????????????????????".format(str(exp_matrix.columns)))
    exp_matrix = exp_matrix.set_index(good_id)
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_exp_matrix_whole_transcriptome(data, option_name, dir_path, bind_obj=None):
    if "__format__task_id" in str(data):
        db = Config().get_mongo_client(mtype='whole_transcriptome')[Config().get_mongo_dbname('whole_transcriptome')]
        main_conn = db['exp']
        main_conn2 = db['specimen_group']
        try:
            samples = list(chain(*list(main_conn2.find_one({"task_id": data.split('__format__')[0]})['specimen_names'])))
            data = str(main_conn.find_one({"task_id": data.split('__format__')[0]})['main_id'])
        except:
            raise Exception("???????????????????????????task_id:{}?????????main_id(???exp_id)".format(data))
    else:
        group_dict = bind_obj.sheet.option('rna_group_dict')
        group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
        samples = list()
        for each in group_dict:
            samples += group_dict[each]

    ##############################################################################################################
    project_type = 'whole_transcriptome'
    # add by xuxi 20210527
    db = selecet_db(project_type, "exp_detail","exp_id", ObjectId(data), bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    conn = db['exp_detail']
    #
    exp_maininfo = db['exp'].find_one({"main_id": ObjectId(data)})
    if exp_maininfo['level'] == "G":
        show_id = "gene_id"
    else:
        show_id = "transcript_id"
    if exp_maininfo['is_rmbe'] == "false":
        is_rmbe_suffix = ""
    else:
        is_rmbe_suffix = "_batch"

    # target_cols = OrderedDict(transcript_id=1, _id=0)
    target_cols = OrderedDict()
    target_cols[show_id] = 1
    target_cols['_id'] = 0
    for each in samples:
        target_cols[each+is_rmbe_suffix] = 1
    exp_records = conn.find({"exp_id": ObjectId(data)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.rename({show_id: 'seq_id'}, axis=1)
    exp_matrix = exp_matrix.set_index('seq_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_exp_matrix_prok(data, option_name, dir_path, bind_obj=None):
    project_type = 'prok_rna'
    # add by fwy 20210128
    db = selecet_db(project_type,"sg_exp_detail","exp_id",data,bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    conn = db['sg_exp_detail']
    group_dict = bind_obj.sheet.option('rna_group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(seq_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    exp_records = conn.find({"exp_id": ObjectId(data)}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('seq_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_relation_tab(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    main_conn = db['sg_p2g_relationship']
    conn = db['sg_p2g_relationship_prvenn']
    output = os.path.join(dir_path, option_name)
    target_cols = dict(
        _id=0,
        related=1,
        failed_proteins=1,
        failed_transcripts=1
    )
    try:
        main_id = main_conn.find_one({"task_id": data})['main_id']
    except:
        bind_obj.set_error("?????????????????????{}????????????????????????".format(data))
    result = conn.find_one({"rela_id": main_id}, target_cols)
    # print(result)
    if not result or 'related' not in result:
        bind_obj.set_error("?????????????????????{}????????????????????????".format(data))
    with open(output, 'w') as rel_w:
        rel_w.write('related' + '\t' + result['related'] + '\n')
        rel_w.write('failed_transcripts' + '\t' + result['failed_transcripts'] + '\n')
        rel_w.write('failed_proteins' + '\t' + result['failed_proteins'] + '\n')
    return output

def export_diff_labelfree(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    # diff_id = bind_obj.sheet.option('protein_diff_id')
    diff_id = data
    compare_group = bind_obj.sheet.option('protein_compare_group')
    target_cols = OrderedDict(accession_id=1, log2fc=1, significant=1, regulate=1,  _id=0)

    bind_obj.logger.debug("?????????????????? {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id),"compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    diff_exp_matrix = diff_exp_matrix.set_index('accession_id')
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export protein diff matrix')
    return output

def export_diff_rna(data, option_name, dir_path, bind_obj=None):
    project_type = data
    # add by fwy 20210128
    table_name = 'sg_diff_detail'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
    db = selecet_db(project_type, table_name, "diff_id",ObjectId(bind_obj.sheet.option('gene_diff_id')), bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    diff_id = bind_obj.sheet.option('gene_diff_id')
    compare_group = bind_obj.sheet.option('gene_compare_group')
    target_cols = OrderedDict(seq_id=1, log2fc=1, significant=1, regulate=1,  _id=0)
    bind_obj.logger.debug("?????????????????? {}".format(target_cols))
    conn = db[table_name]
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id),"compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    diff_exp_matrix = diff_exp_matrix.set_index('seq_id')
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export protein diff matrix')
    return output

def export_fc_info(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    protein_id, rna_id, rna_type = data.split(',')
    db_rna = Config().get_mongo_client(mtype=rna_type)[Config().get_mongo_dbname(rna_type)]
    conn_protein = db['sg_diff']
    conn_rna = db['sg_diff']
    diff_protein = conn_protein.find_one({"_id": ObjectId(protein_id)})
    diff_rna = conn_rna.find_one({"_id": ObjectId(rna_id)})
    if diff_protein and 'params' in diff_protein:
        diff_info = json.loads(diff_protein['params'])
        try:
            protein_up = str(diff_info['fc_up'])
            protein_down = str(diff_info['fc_down'])
        except:
            protein_up = '1.2'
            protein_down = '0.83'
    else:
        protein_up = '1.2'
        protein_down = '0.83'
    if diff_rna and 'params' in diff_rna:
        diff_info = json.loads(diff_rna['params'])
        try:
            rna_up = str(diff_info['fc'])
            rna_down = str(1/float(rna_up))
        except:
            rna_up = '2'
            rna_down = '0.5'
    else:
        rna_up = '2'
        rna_down = '0.5'
    output = os.path.join(dir_path, option_name)
    with open(output, 'w') as fc_w:
        fc_w.write('protein_up' + '\t' + protein_up + '\n')
        fc_w.write('protein_down' + '\t' + protein_down + '\n')
        fc_w.write('rna_up' + '\t' + rna_up + '\n')
        fc_w.write('rna_down' + '\t' + rna_down + '\n')
    return output

def export_relaset_list(data, option_name, dir_path, bind_obj=None):
    relaset_list_path = os.path.join(dir_path, "relaset.list")
    bind_obj.logger.debug("????????????????????????????????????")
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    collection = db['sg_relaset_detail']
    main_collection = db['sg_relaset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("????????????????????????????????????????????????".format(ObjectId(data)))
    result = collection.find_one({"relaset_id": ObjectId(data)}, {'_id': -1, "seq_list": 1})
    with open(relaset_list_path, "wb") as f:
        if 'seq_list' in result and result['seq_list']:
            f.write("\n".join(result['seq_list']))
        else:
            Exception("?????????????????????????????????")
    return relaset_list_path

def export_proteinset_list(data, option_name, dir_path, bind_obj=None):
    proteinset_list_path = os.path.join(dir_path, "proteinset.list")
    bind_obj.logger.debug("????????????????????????????????????")
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    my_result = main_collection.find_one({'main_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("????????????????????????????????????????????????".format(ObjectId(data)))
    result = collection.find_one({"proteinset_id": ObjectId(data)}, {'_id': -1, "seq_list": 1})
    with open(proteinset_list_path, "wb") as f:
        if 'seq_list' in result and result['seq_list']:
            f.write("\n".join(result['seq_list']))
        else:
            Exception("?????????????????????????????????")
    return proteinset_list_path

def export_geneset_list(data, option_name, dir_path, bind_obj=None):
    geneset_list_path = os.path.join(dir_path, "geneset.list")
    bind_obj.logger.debug("????????????????????????????????????")
    main_id, project_type = data.split(',')
    table_name = 'sg_geneset'
    table_name2 = 'sg_geneset_detail'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    # add by fwy 20210128
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    db = selecet_db(project_type, table_name, 'main_id',ObjectId(main_id), bind_obj=bind_obj)
    collection = db[table_name2]
    main_collection = db[table_name]
    my_result = main_collection.find_one({'main_id': ObjectId(main_id)})
    if not my_result:
        bind_obj.set_error("????????????????????????????????????????????????".format(ObjectId(main_id)))
    result = collection.find_one({"geneset_id": ObjectId(main_id)}, {'_id': -1, "seq_list": 1})
    with open(geneset_list_path, "wb") as f:
        if 'seq_list' in result and result['seq_list']:
            f.write("\n".join(result['seq_list']))
        else:
            Exception("?????????????????????????????????")
    return geneset_list_path

def export_rna_protein_group(data, option_name, dir_path, bind_obj=None):
    protein_group_dict = bind_obj.sheet.option('protein_group_dict')
    rna_group_dict = bind_obj.sheet.option('rna_group_dict')
    protein_group_dict = json.loads(protein_group_dict, object_pairs_hook=OrderedDict)
    rna_group_dict = json.loads(rna_group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in protein_group_dict:
            for each in protein_group_dict[key]:
                f.write('{}_protein\t{}_protein\n'.format(each, key))
        for key in rna_group_dict:
            for each in rna_group_dict[key]:
                f.write('{}_rna\t{}_rna\n'.format(each, key))
    return group_out

def export_protein_group(data, option_name, dir_path, bind_obj=None):
    protein_group_dict = bind_obj.sheet.option('protein_group_dict')
    protein_group_dict = json.loads(protein_group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in protein_group_dict:
            for each in protein_group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    return group_out

def export_rna_group(data, option_name, dir_path, bind_obj=None):
    rna_group_dict = bind_obj.sheet.option('protein_group_dict')
    rna_group_dict = json.loads(rna_group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in rna_group_dict:
            for each in rna_group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    return group_out

def get_proteinsets_ids(data, bind_obj):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    proteinset_collection = db["sg_proteinset"]
    proteins = list()
    for proteinset_id in data.split(","):
        proteinset_result = proteinset_collection.find_one({"main_id": ObjectId(proteinset_id)})
        if not proteinset_result:
            bind_obj.set_error("????????????:??????????????????_id???{}??????????????????".format(proteinset_id))
        collection = db['sg_proteinset_detail']
        results = collection.find_one({"proteinset_id": ObjectId(proteinset_id)})
        seqs = results["seq_list"]
        proteins += seqs
    proteins = list(set(proteins))
    return proteins

def get_genesets_ids(data, bind_obj):
    project_type = data.split(',')[0]
    # add by fwy 20210128
    table_name = 'sg_geneset'
    table_name2 = 'sg_geneset_detail'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    db = selecet_db(project_type, table_name, "main_id",ObjectId(data.split(",")[1]), bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    geneset_collection = db[table_name]
    genes = list()
    for geneset_id in data.split(",")[1:]:
        geneset_result = geneset_collection.find_one({"main_id": ObjectId(geneset_id)})
        if not geneset_result:
            bind_obj.set_error("????????????:??????????????????_id???{}??????????????????".format(geneset_id))
        collection = db[table_name2]
        results = collection.find_one({"geneset_id": ObjectId(geneset_id)})
        seqs = results["seq_list"]
        genes += seqs
    genes = list(set(genes))
    return genes

def export_go_protein(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    go_path = os.path.join(dir_path, option_name)
    bind_obj.logger.debug("????????????{}".format(go_path))
    proteins = get_proteinsets_ids(data, bind_obj)
    # go_collection = db["sg_annotation_go"]
    go_level_collection = db["sg_annotation_go_detail"]
    # go_id = go_collection.find_one({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option('type')})["main_id"]
    go_id = bind_obj.sheet.option("protein_go_id")
    one_record = go_level_collection.find_one({'go_id': ObjectId(go_id)})
    if not one_record:
        bind_obj.set_error("????????????:?????????go_id???{}?????????go????????????".format(go_id))
    with open(go_path, "wb") as w:
        w.write("Term type\tTerm\tGO\tLevel\t" + "protein_num\tprotein_percent\tprotein_list" + "\n")
        go_results = go_level_collection.find({'go_id': ObjectId(go_id)})
        for gr in go_results:
            seq_list = gr["seq_list"].split(";")
            total_num = len(proteins)
            in_list = list(set(proteins) & set(seq_list))
            percent = str(float(len(in_list)) / total_num)
            level = str(gr['level'])
            if in_list:
                if level == '2':
                    w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"], level) +
                        "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '3':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_3"], gr["goid_3"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '4':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_4"], gr["goid_4"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
        # term_list = ["biological_process", "cellular_component", "molecular_function"]
    return go_path

def export_go_rna(data, option_name, dir_path, bind_obj=None):
    project_type = data.split(',')[0]
    # add by fwy 20210128
    if project_type in ['denovo_rna_v2',"ref_rna_v2","ref_rna","medical_transcriptome"]:
        select_table = "sg_annotation_go_detail"
    elif project_type == "prok_rna":
        select_table = "sg_annotation_go_graph"
    elif project_type == "whole_transcriptome":
        select_table = "annotation_go_detail"
    else:
        select_table = "sg_annotation_go_detail"
    db = selecet_db(project_type, select_table,'go_id', ObjectId(bind_obj.sheet.option("gene_go_id")), bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    go_path = os.path.join(dir_path, option_name)
    bind_obj.logger.debug("????????????{}".format(go_path))
    genes = get_genesets_ids(data, bind_obj)
    go_id = bind_obj.sheet.option("gene_go_id")
    # go_collection = db["sg_annotation_go"]
    w = open(go_path, "wb")
    w.write("Term type\tTerm\tGO\tLevel\t" + "gene_num\tgene_percent\tgene_list" + "\n")
    if project_type == 'denovo_rna_v2':
        go_level_collection = db["sg_annotation_go_detail"]
        one_record = go_level_collection.find_one({'go_id': ObjectId(go_id), "anno_type": "G"})
        if not one_record:
            bind_obj.set_error("????????????:?????????go_id???{}?????????go????????????".format(go_id))
        go_results = go_level_collection.find({'go_id': ObjectId(go_id), "anno_type" : "G", "level" : 2})
        for gr in go_results:
            seq_list = gr["seq_list"].split(";")
            total_num = len(genes)
            in_list = list(set(genes) & set(seq_list))
            percent = str(float(len(in_list)) / total_num)
            level = str(gr['level'])
            if in_list:
                if level == '2':
                    w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"], level) +
                            "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '3':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_3"], gr["goid_3"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '4':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_4"], gr["goid_4"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
    if project_type == 'ref_rna':
        go_level_collection = db["sg_annotation_go_detail"]
        one_record = go_level_collection.find_one({'go_id': ObjectId(go_id), "anno_type": "gene"})
        if not one_record:
            bind_obj.set_error("????????????:?????????go_id???{}?????????go????????????".format(go_id))
        go_results = go_level_collection.find({'go_id': ObjectId(go_id), "anno_type": "gene", "seq_type": "ref", "level" : 2})
        for gr in go_results:
            seq_list = gr["seq_list"].split(";")
            total_num = len(genes)
            in_list = list(set(genes) & set(seq_list))
            percent = str(float(len(in_list)) / total_num)
            level = str(gr['level'])
            if in_list:
                if level == '2':
                    w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"], level) +
                            "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '3':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_3"], gr["goid_3"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '4':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_4"], gr["goid_4"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
    if project_type == 'ref_rna_v2' or project_type == 'medical_transcriptome':
        go_level_collection = db["sg_annotation_go_detail"]
        one_record = go_level_collection.find_one({'go_id': ObjectId(go_id), "anno_type": "G"})
        if not one_record:
            bind_obj.set_error("????????????:?????????go_id???{}?????????go????????????".format(go_id))
        go_results = go_level_collection.find({'go_id': ObjectId(go_id), "anno_type": "G", "seq_type": "ref", "level" : 2})
        for gr in go_results:
            # print(gr)
            seq_list = gr["seq_list"].split(";")
            total_num = len(genes)
            in_list = list(set(genes) & set(seq_list))
            percent = str(float(len(in_list)) / total_num)
            level = str(gr['level'])
            if in_list:
                if level == '2':
                    w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"], level) +
                            "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '3':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_3"], gr["goid_3"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '4':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_4"], gr["goid_4"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
    if project_type == 'whole_transcriptome':
        go_level_collection = db["annotation_go_detail"]
        one_record = go_level_collection.find_one({'go_id': ObjectId(go_id), "anno_type": "G"})
        if not one_record:
            raise Exception("????????????:?????????go_id???{}?????????go????????????".format(go_id))
        go_results = go_level_collection.find({'go_id': ObjectId(go_id), "anno_type": "G", "seq_type": "ref", "level" : 2})
        for gr in go_results:
            # print(gr)
            seq_list = gr["seq_list"].split(";")
            total_num = len(genes)
            in_list = list(set(genes) & set(seq_list))
            percent = str(float(len(in_list)) / total_num)
            level = str(gr['level'])
            if in_list:
                if level == '2':
                    w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"], level) +
                            "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '3':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_3"], gr["goid_3"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '4':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_4"], gr["goid_4"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
    if project_type == 'prok_rna':
        go_level_collection = db["sg_annotation_go_graph"]
        one_record = go_level_collection.find_one({'go_id': ObjectId(go_id)})
        if not one_record:
            bind_obj.set_error("????????????:?????????go_id???{}?????????go????????????".format(go_id))
        go_results = go_level_collection.find({'go_id': ObjectId(go_id)})
        for gr in go_results:
            # seq_list = gr["seq_list"].split(";")
            seq_list = gr["seq_list"]
            total_num = len(genes)
            in_list = list(set(genes) & set(seq_list))
            percent = str(float(len(in_list)) / total_num)
            try:
                level = str(gr['level'])
            except Exception as e:
                level = "2" # new def 'add_annotation_go_graph_v2'(prokrna_annotation.py) of prokrna ,level always == '2'
            if in_list:
                if level == '2':
                    w.write("{}\t{}\t{}\t{}\t".format(gr["term_type"], gr["go_term"], gr["go_id2"], level) +
                        "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '3':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["term_type"], gr["go_term3"], gr["go_id3"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
                # if level == '4':
                #     w.write("{}\t{}\t{}\t{}\t".format(gr["term_type"], gr["go_term4"], gr["go_id4"], level) +
                #             "{}\t{}\t{}\t".format(str(len(in_list)), percent, ';'.join(in_list)) + '\n')
    w.close()
    return go_path

def export_add_info_protein(data, option_name,dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("????????????add_info??????")
    col = db["sg_annotation_kegg"]
    result = col.find_one({"_id":ObjectId(data)})
    insert_id = result["main_id"]
    col = db["sg_annotation_kegg_level"]
    results = col.find({"kegg_id":insert_id})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info

def export_add_info_rna(data, option_name,dir_path, bind_obj=None):
    project_type, data = data.split(',')
    table_name = 'sg_annotation_kegg'
    table_name2 = 'sg_annotation_kegg_level'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    #modify by fwy 20210128
    db = selecet_db(project_type, table_name, "_id",ObjectId(data), bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("????????????add_info??????")
    col = db[table_name]
    result = col.find_one({"_id":ObjectId(data)})
    insert_id = result["main_id"]
    col = db[table_name2]
    results = col.find({"kegg_id":insert_id, "anno_type":"G"})
    if project_type == 'ref_rna':
        results = col.find({"kegg_id": insert_id, "seq_type": "all", "anno_type": "gene"})
    if project_type == 'prok_rna':
        results = col.find({"kegg_id": insert_id})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info

def export_kegg_table_protein(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    kegg_path = os.path.join(dir_path, 'protein_kegg_table.xls')
    bind_obj.logger.debug("??????????????????%s???kegg_table???????????????:%s" % (option_name, kegg_path))
    my_result = db["sg_annotation_kegg"].find_one({"_id": ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????annotation_kegg_id:{}???sg_annotation_kegg???????????????".format(data))
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Protein id)\tKO_name(Protein name)\tHyperlink\tPaths\tKEGG_gene_id\n')
        kegg_id = my_result["main_id"]
        bind_obj.logger.debug(kegg_id)
        results = db['sg_annotation_kegg_table'].find({'kegg_id': kegg_id})
        one_record = db['sg_annotation_kegg_table'].find_one({'kegg_id': kegg_id})
        if not one_record:
            bind_obj.set_error("??????kegg_table?????????kegg_id:{}???sg_annotation_kegg_table???????????????".format(ObjectId(kegg_id)))
        for result in results:
            if 'hyperlink' not in result:
                bind_obj.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                result['hyperlink'] = 'None'
            w.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths'], result['kegg_gene_id']))
    return kegg_path

def export_kegg_table_rna(data, option_name, dir_path, bind_obj=None):
    project_type, data = data.split(',')
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    #modify by fwy 20210128
    table_name = 'sg_annotation_kegg'
    table_name2 = 'sg_annotation_kegg_table'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    db = selecet_db(project_type, table_name,"_id", ObjectId(data), bind_obj=bind_obj)
    kegg_path = os.path.join(dir_path, 'gene_kegg_table.xls')
    bind_obj.logger.debug("??????????????????%s???kegg_table???????????????:%s" % (option_name, kegg_path))
    my_result = db[table_name].find_one({"_id": ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????annotation_kegg_id:{}???sg_annotation_kegg???????????????".format(data))
    with open(kegg_path, 'wb') as w:
        w.write('#Query\tKO_ID(Protein id)\tKO_name(Protein name)\tHyperlink\tPaths\tKEGG_gene_id\n')
        kegg_id = my_result["main_id"]
        bind_obj.logger.debug(kegg_id)
        results = db[table_name2].find({'kegg_id': kegg_id, "anno_type": "G"})
        if project_type == 'ref_rna':
            results = db[table_name2].find({"kegg_id": kegg_id, "anno_type": "gene"})
        if project_type == 'prok_rna':
            results = db[table_name2].find({"kegg_id": kegg_id})
        one_record = db[table_name2].find_one({'kegg_id': kegg_id})
        if not one_record:
            bind_obj.set_error("??????kegg_table?????????kegg_id:{}???sg_annotation_kegg_table???????????????".format(ObjectId(kegg_id)))
        for result in results:
            if 'hyperlink' not in result:
                bind_obj.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                result['hyperlink'] = 'None'
            w.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'], result['hyperlink'], result['paths'], ""))
    return kegg_path

def export_kegg_level_table_protein(data, option_name, dir_path, bind_obj=None):
    kegg__level_path = os.path.join(dir_path, 'protein_kegg_level_table.xls')
    bind_obj.logger.debug("??????????????????%s???kegg_table???????????????:%s" % (option_name, kegg__level_path))
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    my_result = db["sg_annotation_kegg"].find_one({"_id": ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????annotation_kegg_id:{}???sg_annotation_kegg???????????????".format(data))
    with open(kegg__level_path, 'wb') as w:
        w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        kegg_id = my_result["main_id"]
        bind_obj.logger.debug(kegg_id)
        if not kegg_id:
            bind_obj.set_error("???????????????annotation_kegg_id:{}???sg_annotation_kegg???????????????".format(kegg_id))
        results = db["sg_annotation_kegg_level"].find({"kegg_id": kegg_id})
        one_record = db['sg_annotation_kegg_level'].find_one({'kegg_id': kegg_id})
        if not one_record:
            bind_obj.set_error("??????kegg_table?????????kegg_id:{}???sg_annotation_kegg_level???????????????".format(ObjectId(kegg_id)))
        for result in results:
            if 'hyperlink' not in result:
                bind_obj.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                result['hyperlink'] = 'None'
            w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['seq_num'], result['pathway_definition'], result['first_category'], '', result['hyperlink'], result['seq_list'], '', result['second_category']))
    return kegg__level_path

def export_kegg_level_table_rna(data, option_name, dir_path, bind_obj=None):
    kegg__level_path = os.path.join(dir_path, 'gene_kegg_level_table.xls')
    bind_obj.logger.debug("??????????????????%s???kegg_table???????????????:%s" % (option_name, kegg__level_path))
    project_type, data = data.split(',')
    # add by fwy 20210128
    table_name = 'sg_annotation_kegg'
    table_name2 = 'sg_annotation_kegg_level'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    db = selecet_db(project_type, table_name, "_id",ObjectId(data), bind_obj=bind_obj)
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    my_result = db[table_name].find_one({"_id": ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????annotation_kegg_id:{}???sg_annotation_kegg???????????????".format(data))
    with open(kegg__level_path, 'wb') as w:
        w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
        kegg_id = my_result["main_id"]
        bind_obj.logger.debug(kegg_id)
        if not kegg_id:
            bind_obj.set_error("???????????????annotation_kegg_id:{}???sg_annotation_kegg???????????????".format(kegg_id))
        results = db[table_name2].find({"kegg_id": kegg_id, "seq_type": "all", "anno_type": "G"})
        if project_type == 'ref_rna':
            results = db[table_name2].find({"kegg_id": kegg_id, "seq_type": "all", "anno_type": "gene"})
        if project_type == 'prok_rna':
            results = db[table_name2].find({"kegg_id": kegg_id})
        if project_type == 'denovo_rna_v2':
            results = db[table_name2].find({"kegg_id": kegg_id, "anno_type": "G"})
        one_record = db[table_name2].find_one({'kegg_id': kegg_id})
        if not one_record:
            bind_obj.set_error("??????kegg_table?????????kegg_id:{}???sg_annotation_kegg_level???????????????".format(ObjectId(kegg_id)))
        for result in results:
            if 'hyperlink' not in result:
                bind_obj.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                result['hyperlink'] = 'None'
            if project_type == 'prok_rna':
                # w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['seq_num'], result['pathway_definition'], result['first_category'], '', result['hyperlink'], result['seq_list'], result['anno_type'], result['second_category']))
                w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['seq_num'], result['pathway_definition'], result['first_category'], '', result['hyperlink'], result['seq_list'], '', result['second_category']))
            else:
                w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '', result['number_of_seqs'],
                                                                          result['pathway_definition'],
                                                                          result['first_category'], '',
                                                                          result['hyperlink'], result['seq_list'],
                                                                          result['anno_type'],
                                                                          result['second_category']))
    return kegg__level_path

def export_genepro_list(data, option_name, dir_path, bind_obj=None):
    list_path = os.path.join(dir_path, "gene_protein.list")
    project_type, genesets, proteinsets = data.split('___')
    bind_obj.logger.debug("??????????????????????????????????????????")
    db = Config().get_mongo_client(mtype='labelfree')[Config().get_mongo_dbname('labelfree')]
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    proteins = list()
    for proteinset in proteinsets.split(','):
        p_result = main_collection.find_one({'main_id': ObjectId(proteinset)})
        if not p_result:
            bind_obj.set_error("????????????????????????????????????????????????".format(ObjectId(proteinset)))
        result = collection.find_one({"proteinset_id": ObjectId(proteinset)}, {'_id': -1, "seq_list": 1})
        proteins += result['seq_list']

    genes = list()
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    #modify by fwy 2021028
    table_name = 'sg_geneset'
    table_name2 = 'sg_geneset_detail'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    rna_db = selecet_db(project_type, table_name,'main_id' ,ObjectId(genesets.split(',')[0]), bind_obj=bind_obj)
    collection = rna_db[table_name2]
    main_collection = rna_db[table_name]
    for geneset in genesets.split(','):
        p_result = main_collection.find_one({'main_id': ObjectId(geneset)})
        if not p_result:
            bind_obj.set_error("????????????????????????????????????????????????".format(ObjectId(geneset)))
        result = collection.find_one({"geneset_id": ObjectId(geneset)}, {'_id': -1, "seq_list": 1})
        genes += result['seq_list']

    with open(list_path, "wb") as f:
        f.write('protein\t' + ','.join(set(proteins)) + '\n')
        f.write('gene\t' + ','.join(set(genes)) + '\n')
    return list_path

def export_proteinsets_list(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    out_paths = list()
    for proteinset in data.split(','):
        collection = db['sg_proteinset_detail']
        main_collection = db['sg_proteinset']
        my_result = main_collection.find_one({'main_id': ObjectId(proteinset)})
        if not my_result:
            bind_obj.set_error("???????????????proteinset_id:{}???sg_proteinset???????????????".format(ObjectId(data)))
        name = my_result['name']
        bind_obj.logger.debug("?????????????????????-->%s"%name)
        results = collection.find_one({"proteinset_id": ObjectId(proteinset)})
        protein_list_path = os.path.join(dir_path, "%s_protein.list" % name)
        with open(protein_list_path, "wb") as f:
            protein_list = results["seq_list"]
            for protein_id in protein_list:
                f.write(protein_id + "\n")
        out_paths.append(protein_list_path)
    return ','.join(out_paths)

def export_genesets_list(data, option_name, dir_path, bind_obj=None):
    project_type, data = data.split(';')
    table_name = 'sg_geneset'
    table_name2 = 'sg_geneset_detail'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    #modify by fwy 20210128
    db = selecet_db(project_type, table_name,'main_id', ObjectId(data.split(',')[0]), bind_obj=bind_obj)
    out_paths = list()
    for geneset in data.split(','):
        collection = db[table_name2]
        main_collection = db[table_name]
        my_result = main_collection.find_one({'main_id': ObjectId(geneset)})
        if not my_result:
            bind_obj.set_error("???????????????geneset_id:{}???sg_geneset???????????????".format(ObjectId(data)))
        name = my_result['name']
        bind_obj.logger.debug("?????????????????????-->%s"%name)
        results = collection.find_one({"geneset_id": ObjectId(geneset)})
        gene_list_path = os.path.join(dir_path, "%s_gene.list" % name)
        with open(gene_list_path, "wb") as f:
            gene_list = results["seq_list"]
            for gene_id in gene_list:
                f.write(gene_id + "\n")
        out_paths.append(gene_list_path)
    return ','.join(out_paths)

def export_go_list_rna(data, option_name, dir_path, bind_obj=None):
    project_type, data = data.split(';')
    table_name = 'sg_annotation_go'
    table_name2 = 'sg_annotation_go_list'
    if project_type == "whole_transcriptome":
        table_name = table_name[3:]
        table_name2 = table_name2[3:]
    # db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    db = selecet_db(project_type, table_name,"_id", ObjectId(data), bind_obj=bind_obj)
    go_list_path = os.path.join(dir_path, "GO_rna.list")
    bind_obj.logger.debug("????????????%s??????:%s" % (option_name, go_list_path))
    my_result = db[table_name].find_one({"_id": ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????annotation_go_id:{}???sg_annotation_go???????????????".format(data))
    collection = db[table_name2]
    results = collection.find({"go_id": ObjectId(data), "anno_type": "G"})
    if project_type == 'ref_rna':
        results = collection.find({"go_id": ObjectId(data), "anno_type": "gene"})
    if project_type == 'prok_rna':
        results = collection.find({"go_id": ObjectId(data)})
    one_record = collection.find_one({"go_id": ObjectId(data)})
    if not one_record:
        bind_obj.set_error("??????gos_list?????????annotation_id:{}???sg_annotation_gos_list???????????????".format(ObjectId(data)))
    with open(go_list_path, "wb") as w:
        for result in results:
            protein_id = result["gene_id"]
            go_list = result["gos_list"]
            w.write(protein_id + "\t" + go_list + "\n")
    return go_list_path

def export_go_list_protein(data, option_name, dir_path, bind_obj=None):
    project_type = 'labelfree'
    db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
    go_list_path = os.path.join(dir_path, "GO_protein.list")
    bind_obj.logger.debug("????????????%s??????:%s" % (option_name, go_list_path))
    my_result = db["sg_annotation_go"].find_one({"_id": ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????annotation_go_id:{}???sg_annotation_go???????????????".format(data))
    collection = db["sg_annotation_go_list"]
    results = collection.find({"go_id": ObjectId(data)})
    one_record = collection.find_one({"go_id": ObjectId(data)})
    if not one_record:
        bind_obj.set_error("??????gos_list?????????annotation_id:{}???sg_annotation_gos_list???????????????".format(ObjectId(data)))
    with open(go_list_path, "wb") as w:
        for result in results:
            protein_id = result["gene_id"]
            go_list = result["gos_list"]
            w.write(protein_id + "\t" + go_list + "\n")
    return go_list_path

# ???????????????id???sg_proteinset?????????task_id???type???????????????????????????sg_annotation_kegg???????????????????????????????????????????????????????????????2???proteinset_id?????????
# ??????????????????????????????????????????proteinset_id?????????????????????????????????????????????????????????
# ?????????one_record?????????????????????find???????????????????????????????????????????????????




def export_kegg_pdf(data, option_name, dir_path, bind_obj=None):
    fs = gridfs.GridFS(db)
    annotation_collection = db["sg_annotation_kegg"]
    task_id = data.split("\t")[0]
    anno_type = data.split("\t")[1]
    main_id = annotation_collection.find_one({"task_id":task_id, "type": bind_obj.sheet.option("type")})["main_id"]
    kegg_level_collection = db["sg_annotation_kegg_level"]
    results = kegg_level_collection.find({"kegg_id":main_id, "anno_type":anno_type})
    anno_path = dir_path + "/png"
    if not os.path.exists(anno_path):
        os.mkdir(anno_path)
    for result in results:
        graph_id = result["graph_png_id"]
        pathway_id = result["pathway_id"]
        bind_obj.logger.info("????????????{}???png??????".format(pathway_id))
        with open(anno_path + "/" + pathway_id + ".png", "w") as fw:
            content = fs.get(graph_id).read()
            fw.write(content)
    return anno_path


def export_all_list(data, option_name, dir_path, bind_obj=None):
    all_list = os.path.join(dir_path, "all_protein.list")
    bind_obj.logger.debug("??????????????????????????????{}".format(all_list))
    collection = db['sg_express_detail']
    exp_collection = db['sg_express']
    main_collection = db['sg_proteinset']
    my_result = main_collection.find_one({"main_id": ObjectId(data)})
    task_id = my_result["task_id"]
    # proteinset_type = my_result["type"]
    bind_obj.logger.debug(task_id)
    exp_result = exp_collection.find_one({'task_id': bind_obj.sheet.option("task_id")})
    if not exp_result:
        bind_obj.set_error("???????????????task_id:{}??????????????????sg_proteinset???????????????".format(data))
    exp_id = exp_result["main_id"]
    results = collection.find({"express_id": ObjectId(exp_id)})
    with open(all_list, "wb") as f:
        for result in results:
            protein_id = result['seq_id']
            f.write(protein_id + "\n")
    return all_list


def export_cog_class(data, option_name, dir_path, bind_obj=None):
    cog_path = os.path.join(dir_path, 'cog_class_table.xls')
    bind_obj.logger.debug("????????????{}".format(cog_path))
    proteinsets, table_title, task_id, seq_type = get_proteinset_detail(data, bind_obj)
    cog_collection = db["sg_annotation_cog"]
    cog_detail_collection = db["sg_annotation_cog_detail"]
    cog_id = cog_collection.find_one({"task_id": bind_obj.sheet.option("task_id"), 'type': bind_obj.sheet.option('type')})["main_id"]
    print("cog_id:", cog_id)
    cog_results = cog_detail_collection.find({'cog_id': cog_id})
    new_table_title = []
    for tt in table_title:
        # new_tt = [tt + "_COG", tt + "_NOG", tt + "_KOG", tt + "_COG_list", tt + "_NOG_list", tt + "_KOG_list"]
        new_tt = [tt + "_COG", tt + "_NOG", tt + "_COG_list", tt + "_NOG_list"]
        new_table_title = new_table_title + new_tt
    bind_obj.logger.debug(table_title)
    with open(cog_path, "wb") as w:
        w.write("Type\tFunctional Categoris\t" + "\t".join(new_table_title) + "\n")
        for cr in cog_results:
            # kog_list = set([])
            nog_list = set(cr["nog_list"].split(";") if cr["nog_list"] else [])
            cog_list = set(cr["cog_list"].split(";") if cr["cog_list"] else [])
            write_line = {}
            for gt in proteinsets:
                # kog_count = list(kog_list & proteinsets[gt][1])
                nog_count = list(nog_list & proteinsets[gt][1])
                cog_count = list(cog_list & proteinsets[gt][1])
                # if not len(kog_count) + len(nog_count) + len(cog_count) == 0:
                #     write_line[gt] = [str(len(cog_count)), str(len(nog_count)), str(len(kog_count)), ";".join(cog_count), ";".join(nog_count), ";".join(kog_count)]
                if not len(nog_count) + len(cog_count) == 0:
                    write_line[gt] = [str(len(cog_count)), str(len(nog_count)), ";".join(cog_count), ";".join(nog_count)]
            print('here:',  nog_list)
            print('here:',  proteinsets[gt][1])

            if len(write_line) > 0:
                w.write("{}\t{}\t".format(cr["type"], cr["function_categories"]))
                for tt in table_title:
                    w.write("\t".join(write_line[tt]) + "\t") if tt in write_line else w.write("0\t0\tnone\tnone\t")
                w.write("\n")
    return cog_path


def get_proteinset_detail(data, bind_obj):
    proteinset_collection = db["sg_proteinset"]
    proteinsets = {}
    names = []
    task_id = ""
    proteinset_type = "P"
    for proteinset_id in data.split(","):
        proteinset_result = proteinset_collection.find_one({"main_id": ObjectId(proteinset_id)})
        if not proteinset_result:
            bind_obj.set_error("????????????:??????????????????_id???{}??????????????????".format(proteinset_id))
        task_id = proteinset_result["task_id"]
        #proteinset_type = proteinset_result["type"]
        proteinset_name = proteinset_result["name"]
        names.append(proteinset_name)
        proteinsets[proteinset_name] = [proteinset_type]
        collection = db['sg_proteinset_detail']
        results = collection.find_one({"proteinset_id": ObjectId(proteinset_id)})
        proteinset_names = set(results["seq_list"])
        proteinsets[proteinset_name].append(proteinset_names)
    # print(proteinsets)
    return proteinsets, names, task_id, proteinset_type


def export_go_class(data, option_name, dir_path, bind_obj=None):
    go_path = os.path.join(dir_path, 'go_class_table.xls')
    bind_obj.logger.debug("????????????{}".format(go_path))
    proteinsets, names, task_id, seq_type = get_proteinset_detail(data, bind_obj)
    bind_obj.logger.debug("### ????????????{} {} {} {}".format(proteinsets, names, task_id, seq_type))
    go_collection = db["sg_annotation_go"]
    go_level_collection = db["sg_annotation_go_detail"]
    go_id = go_collection.find_one({"task_id": bind_obj.sheet.option("task_id"), "type": bind_obj.sheet.option('type')})["main_id"]
    one_record = go_level_collection.find_one({'go_id': go_id, "level": 2})
    if not one_record:
        bind_obj.set_error("????????????:?????????go_id???{}??????????????????".format(go_id))
    new_table_title = []
    for gt in proteinsets:
        new_table_title.append(gt + " num")
        new_table_title.append(gt + " percent")
        new_table_title.append(gt + " list")
    bind_obj.logger.debug(new_table_title)
    with open(go_path, "wb") as w:
        w.write("Term type\tTerm\tGO\t" + "\t".join(new_table_title) + "\n")
        term_list = ["biological_process", "cellular_component", "molecular_function"]
        for item in term_list:
            go_results = go_level_collection.find({'go_id': go_id, "level": 2})
            for gr in go_results:
                if gr["goterm"] == item:
                    seq_list = set(gr["seq_list"].split(";"))
                    write_line = {}
                    for gt in proteinsets:
                        total_protein_num = len(proteinsets[gt][1])
                        go_count = list(seq_list & proteinsets[gt][1])
                        if not len(go_count) == 0:
                            write_line[gt] = str(len(go_count)) + "\t" + str(len(go_count)/total_protein_num) + \
                                             "(" + str(len(go_count)) + "/" + str(total_protein_num) + ")" + "\t" + ";".join(go_count)
                    if len(write_line):
                        w.write("{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"]))
                        for tt in proteinsets:
                            w.write(write_line[tt] + "\t") if tt in write_line else w.write("0\t0\tnone\t")
                        w.write("\n")
    return go_path


# ############???????????????

def export_group_table_by_detail(data, option_name, dir_path, bind_obj=None):
    """
    ??????????????????????????????group???
    ?????????????????????workflow???option???group_detal????????????
    """
    file_path = os.path.join(dir_path, "%s_input.group.xls" % option_name)
    bind_obj.logger.debug("??????????????????%s???GROUP????????????????????????:%s" % (option_name, file_path))
    if data in ["all", "All", "ALL"]:
        with open(file_path, "wb") as f:
            f.write("#sample\t" + "##empty_group##" + "\n")
        return file_path
    data = _get_objectid(data, bind_obj=bind_obj)
    group_detail = bind_obj.sheet.option('group_detail')  #????????????
    group_table = db['sg_specimen_group']
    if not isinstance(group_detail, dict):
        try:
            table_dict = json.loads(group_detail)
        except Exception:
            bind_obj.set_error("??????group?????????????????????{}???????????????????????????????????????????????????".format(option_name))
    if not isinstance(table_dict, dict):
        bind_obj.set_error("??????group?????????????????????{}???????????????????????????????????????????????????".format(option_name))
    group_schema = group_table.find_one({"main_id": ObjectId(data)})
    if not group_schema:
        bind_obj.set_error("?????????????????????group_id:{}???sg_specimen_group_compare???????????????????????????".format(data))
    schema_name = re.sub("\s", "_", group_schema["group_name"])  # ?????????????????????????????????????????????
    with open(file_path, "wb") as f:
        f.write("#sample\t" + "group" + "\n")
    sample_table_name = 'sg_specimen'
    sample_table = db[sample_table_name]
    with open(file_path, "ab") as f:
        for k in table_dict:
            for sp_id in table_dict[k]:
                sp = sample_table.find_one({"main_id": ObjectId(sp_id)})
                if not sp:
                    bind_obj.set_error("group_detal????????????_id:{}????????????{}????????????".format(sp_id, sample_table_name))
                else:
                    sp_name = sp["specimen_name"]
                f.write("{}\t{}\n".format(sp_name, k))
    return file_path


def _get_objectid(data, bind_obj=None):
    if not isinstance(data, ObjectId):
        if not isinstance(data, StringTypes):
            bind_obj.set_error("{}??????ObjectId?????????????????????????????????".format(data))
        else:
            try:
                data = ObjectId(data)
            except:
                bind_obj.set_error("{}??????ObjectId?????????????????????????????????".format(data))
    return data


def export_control_file(data, option_name, dir_path, bind_obj=None):  #??????????????? ????????????
    file_path = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("????????????????????????:%s" % file_path)
    collection = db['sg_specimen_group_compare']
    result = collection.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("???????????????control_id:{}???sg_specimen_group_compare???????????????".format(ObjectId(data)))
    group_id = result['specimen_group_id']
    if group_id not in ['all', 'All', 'ALL']:
        """??????group_id?????????"""
        if isinstance(group_id, types.StringTypes):
            group_id = ObjectId(group_id)
        group_coll = db['sg_specimen_group']
        g_result = group_coll.find_one({'main_id': group_id})
        if not g_result:
            bind_obj.set_error("???????????????control_file???group_id:{}???sg_specimen_group???????????????".format(group_id))
    control_detail = json.loads(result['compare_names'])
    with open(file_path, 'wb') as w:
        w.write('#control\t{}\n'.format('group'))
        for i in control_detail:    #??????????????????, ??????????????????
            # w.write('{}\t{}\n'.format(i.keys()[0], i.values()[0]))
            control_other = i.split("|")
            w.write('{}\t{}\n'.format(control_other[0], control_other[1]))
    return file_path


def _get_protein_id(proteinset, proteinset_detail, _id):
    try:
        results = proteinset_detail.find_one({"proteinset_id":ObjectId(_id)})
        seq_id = results['seq_list']
    except Exception:
        bind_obj.set_error("{}???sg_proteinset_detail??????????????????!")
    try:
        my_result = proteinset.find_one({"main_id":ObjectId(_id)})
        _name = my_result['name']
    except Exception:
        bind_obj.set_error("{}???sg_proteinset??????????????????!")
    return seq_id, _name

# ??????"add_info":proteinset_info['task_id'] + "\t" + data.proteinset_type?????????sg_proteinset?????????task_id????????????????????????annotation??????annotation1
# ????????????origin??????latest
def export_add_info(data, option_name,dir_path, bind_obj=None):
    task_id = data.split("\t")[0]
    # anno_type = data.split("\t")[1]
    add_info = os.path.join(dir_path, '{}.txt'.format(option_name))
    bind_obj.logger.debug("????????????add_info??????")
    col = db["sg_annotation_kegg"]
    result = col.find_one({"task_id":task_id, "type": bind_obj.sheet.option('type')})
    insert_id = result["main_id"]
    col = db["sg_annotation_kegg_level"]
    results = col.find({"kegg_id":insert_id})
    with open(add_info, "w") as fw:
        fw.write("pathway\thyperlink\n")
        for result in results:
            fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
    return add_info


def export_multi_protein_list(data, option_name, dir_path, bind_obj=None):
    proteinset_id = data.split(",")
    multi_proteinset_path = dir_path + "/multi_proteinset_list"
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    f = open(multi_proteinset_path, "wb")
    for n, gi in enumerate(proteinset_id):
        my_result = main_collection.find_one({'main_id': ObjectId(gi)})
        if not my_result:
            bind_obj.set_error("???????????????proteinset_id:{}???sg_proteinset???????????????".format(ObjectId(gi)))
        f.write(my_result["name"] + "\t")
        results = collection.find_one({"proteinset_id": ObjectId(gi)})
        f.write(",".join(results["seq_list"]) + "\n")
    return multi_proteinset_path


# ---------------------???????????????gdq----------------------------------


def export_exp_matrix_scaled(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_express_detail']
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1
    con_main = db['sg_express']
    express_id = con_main.find_one({"task_id": data, "type": "scaled"})["main_id"]
    exp_records = conn.find({"express_id": express_id}, target_cols)
    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('accession_id')
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output


def export_group(data, option_name, dir_path, bind_obj=None):
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    group_out = os.path.join(dir_path, option_name)
    with open(group_out, 'w') as f:
        f.write('#sample\tgroup\n')
        for key in group_dict:
            for each in group_dict[key]:
                f.write('{}\t{}\n'.format(each, key))
    return group_out


def export_compare(data, option_name, dir_path, bind_obj=None):
    conn = db['sg_specimen_group_compare']
    result = conn.find_one({'main_id': ObjectId(data)})
    if not result:
        bind_obj.set_error("control_id: {} is not found in sg_specimen_group_compare".format(ObjectId(data)))
    cmp_info = json.loads(result['compare_names'])
    cmp_out = os.path.join(dir_path, option_name)
    with open(cmp_out, 'w') as f:
        f.write('#control\tother\n')
        for each in cmp_info:
            f.write(each.replace('|', '\t') + '\n')
    return cmp_out
# ---------------------???????????????----------------------------------


# ---------------------??????????????????gdq----------------------------------
def export_proteinset_exp_matrix(data, option_name, dir_path, bind_obj=None):
    exp_id, proteinset_id = data.split(";")
    group_dict = bind_obj.sheet.option('group_dict')
    group_dict = json.loads(group_dict, object_pairs_hook=OrderedDict)
    samples = list()
    for each in group_dict:
        samples += group_dict[each]
    target_cols = OrderedDict(accession_id=1, _id=0)
    for each in samples:
        target_cols[each] = 1

    bind_obj.logger.debug("?????????????????? {}".format(target_cols))
    # get proteinset
    conn = db['sg_proteinset_detail']
    proteinset_records = conn.find_one({"proteinset_id": ObjectId(proteinset_id)})
    proteinset = proteinset_records['seq_list']
    # get exp matrix
    conn = db['sg_express_detail']
    exp_records = conn.find({"express_id": ObjectId(exp_id)}, target_cols)

    bind_obj.logger.debug("????????????express_id??? {}".format(exp_id))

    exp_matrix = pd.DataFrame(list(exp_records))
    exp_matrix = exp_matrix.set_index('accession_id')
    exp_matrix = exp_matrix.loc[proteinset, :]
    output = os.path.join(dir_path, option_name)
    exp_matrix.to_csv(output, sep='\t', header=True, index=True)
    print('success to export expression matrix')
    return output

def export_compare_exp_fc(data, option_name, dir_path, bind_obj=None):
    '''
    ??????????????????accession_id ??? fc
    '''

    diff_id = bind_obj.sheet.option('diff_id')
    compare_group = bind_obj.sheet.option('compare_group')
    target_cols = OrderedDict(accession_id=1, log2fc=1, _id=0)

    bind_obj.logger.debug("?????????????????? {}".format(target_cols))
    conn = db['sg_diff_detail']
    diff_exp_records = conn.find({"diff_id": ObjectId(diff_id),"compare": compare_group}, target_cols)
    diff_exp_matrix = pd.DataFrame(list(diff_exp_records))
    output = os.path.join(dir_path, option_name)
    diff_exp_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export expression matrix')
    return output


def export_kegg_enrich_matrix(data, option_name, dir_path, bind_obj=None):
    '''
    ??????kegg????????????
    '''
    kegg_enrich_id = bind_obj.sheet.option('enrich_id')
    target_cols = OrderedDict(_id=0, id=1, term=1, pvalue=1, corrected_pvalue=1, seq_list=1, kegg_type=1)
    bind_obj.logger.debug("??????KEGG?????????")
    # get proteinset
    conn = db['sg_proteinset_kegg_enrich_detail']
    bind_obj.logger.debug("???????????????{}".format(target_cols))
    bind_obj.logger.debug("???????????????{}".format(kegg_enrich_id))
    kegg_enrich_records = conn.find({"kegg_enrich_id": ObjectId(kegg_enrich_id)}, target_cols)
    kegg_enrich_matrix = pd.DataFrame(list(kegg_enrich_records))
    # exp_matrix = exp_matrix.loc[proteinset, :]
    output = os.path.join(dir_path, option_name)
    kegg_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export kegg_enrich matrix')
    return output

def export_go_enrich_matrix(data, option_name, dir_path, bind_obj=None):

    go_enrich_id = bind_obj.sheet.option('enrich_id')
    go_type = bind_obj.sheet.option('go_type')

    target_cols = OrderedDict(go_id=1, go_type=1, discription=1, p_corrected=1, p_uncorrected=1, seq_list=1, depth=1, _id=0)
    bind_obj.logger.debug("??????GO?????????")
    # get proteinset
    conn = db['sg_proteinset_go_enrich_detail']
    if go_type == 'ALL' or go_type == 'All':
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id)}, target_cols)
    else:
        go_enrich_records = conn.find({"go_enrich_id": ObjectId(go_enrich_id), "go_type":go_type}, target_cols)
    go_enrich_matrix = pd.DataFrame(list(go_enrich_records))
    # exp_matrix = exp_matrix.loc[proteinset, :]
    output = os.path.join(dir_path, option_name)
    go_enrich_matrix.to_csv(output, sep='\t', header=True, index=False)
    print('success to export go_enrich matrix')
    return output

def export_protein_list_ppi(data, option_name, dir_path, bind_obj=None):
    gene_list_path = os.path.join(dir_path, "%s.txt" % option_name)
    bind_obj.logger.debug("?????????????????????")
    collection = db['sg_proteinset_detail']
    main_collection = db['sg_proteinset']
    my_result = main_collection.find_one({'_id': ObjectId(data)})
    if not my_result:
        bind_obj.set_error("???????????????proteinset_id:{}???sg_proteinset???????????????".format(ObjectId(data)))
    results = collection.find_one({"proteinset_id": ObjectId(data)})["seq_list"]
    with open(gene_list_path, "wb") as f:
        f.write("accession_id" + "\n")
        for result in results:
            f.write(result + "\n")
    bind_obj.logger.debug("????????????????????????")
    return gene_list_path


def selecet_db(project_type,table_name,id_name,table_id,bind_obj=None):
    if isinstance(table_id, types.StringTypes):
        main_id = ObjectId(table_id)
    elif isinstance(table_id, ObjectId):
        main_id = table_id
    else:
        bind_obj.set_error("table_id??????????????????????????????ObjectId??????!")
    bind_obj.logger.info("project_type={}\n".format(project_type))
    bind_obj.logger.info("table_name={}\n".format(table_name))
    bind_obj.logger.info("id_name={}\n".format(id_name))
    bind_obj.logger.info("table_id={}\n".format(main_id))
    # bind_obj.logger.info("project_type={}\n".format(project_type))
    raw_workflow_db_version = bind_obj.config.DBVersion
    bind_obj.logger.info("?????????db_version={}\n".format(bind_obj.config.DBVersion))
    try:
        # bind_obj.config.DBVersion = 0
        old_config = Config()
        old_config.DBVersion = 0
        db =  Config().get_mongo_client(mtype=project_type, db_version=0)[ Config().get_mongo_dbname(project_type, db_version=0)]
        target_collection =db[table_name].find_one({id_name: main_id})
        target_collection["_id"]
        bind_obj.config.DBVersion = raw_workflow_db_version
        return db
    except:
        # bind_obj.config.DBVersion = 1
        new_config = Config()
        new_config.DBVersion = 1
        db = Config().get_mongo_client(mtype=project_type, db_version=1)[Config().get_mongo_dbname(project_type, db_version=1)]
        target_collection = db[table_name].find_one({id_name: main_id})
        target_collection["_id"]
        bind_obj.config.DBVersion = raw_workflow_db_version
        return db

    if not target_collection:
        bind_obj.set_error("????????????{}????????????id???{}??????".format(table_name,table_id))


def main():
    pass
