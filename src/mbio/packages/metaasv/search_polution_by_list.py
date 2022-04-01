# -*- coding:utf-8 -*-
##__author__ = 'guanqing.zou'
## modify by qingchen.zhang

import pandas as pd
import numpy as np
import pymongo
import os
from biocluster.config import Config
import datetime


def get_polution(ori_file,search_col,polution_list,threshold_dic,outdir=None):
    data = pd.read_table(ori_file,sep='\t',header=0)
    sample_cols = data.columns.tolist()
    sample_cols.remove(search_col)
    sample_num = len(sample_cols)
    sub_data = data[data[search_col].apply(lambda x: x in polution_list)]
    if outdir:
        sub_data.to_csv(outdir+'/sub_abund.xls',sep='\t',index=False)

    ret_dic = {}
    for p in polution_list:
        ret_dic[p] = {
            'max_abund':0,
            'max_sample':'',
            'over_ave': 0,
            'over_samples' : [],
            'over_num' : 0,
            'pol_num' : 0,
            'pol_samples' : [],
            'is_pol' : 'no'
        }

    for i in range(len(sub_data)):
        line_data = sub_data.iloc[i,]
        pol_name = line_data[search_col]
        #print pol_name
        threshold_value = threshold_dic[pol_name]  #不同的菌，不同的阈值
        max_abund = 0
        max_name = ''
        pol_samples = []
        over_vaule = []
        over_samples = []
        for sample in sample_cols:
            tmp_value = float(line_data[sample])
            if tmp_value > max_abund:
                max_abund = tmp_value
                max_name = sample

            if tmp_value > threshold_value:
                over_samples.append(sample)
                over_vaule.append(tmp_value)

            if tmp_value > 0:
                pol_samples.append(sample)


        ret_dic[pol_name]['max_sample'] = max_name
        ret_dic[pol_name]['max_abund'] = max_abund
        ret_dic[pol_name]['pol_num'] = len(pol_samples)
        ret_dic[pol_name]['pol_samples'] = pol_samples
        ret_dic[pol_name]['over_samples'] = over_samples
        ret_dic[pol_name]['over_num'] = len(over_samples)
        if over_samples:
            ret_dic[pol_name]['over_ave'] = np.mean(over_vaule)
        if  max_abund >  threshold_value:
            ret_dic[pol_name]['is_pol'] = 'yes'

    ##输出ret_dic 到文件
    if outdir:
        with open(outdir+'/pol.xls','w') as fw:
            head = 'pollution\tall_num\tpol_num\tpol_samples\tover_num\tover_ave\tover_samples\tmax_abund\tmax_sample\tis_pol'
            fw.write(head +'\n')
            for pk in sorted(ret_dic.keys()):
                out_list = [
                    pk,
                    str(sample_num),
                    str(ret_dic[pk]['pol_num']),
                    ','.join(ret_dic[pk]['pol_samples']),
                    str(ret_dic[pk]['over_num']),
                    str(ret_dic[pk]['over_ave']),
                    ','.join(ret_dic[pk]['over_samples']),
                    str(ret_dic[pk]['max_abund']),
                    ret_dic[pk]['max_sample'],
                    str(ret_dic[pk]['is_pol'])
                ]

                fw.write('\t'.join(out_list)+'\n')

    return ret_dic

# database : mongo 数据库
def import_mongo(infile,task_id,project_type,threshold_dic,database=None,samples_abunds_file=None,database_detail=None):
    pollution_list = ['pseud','cyano','mitoc','cyan_ch']
    after = ['_max','_spe','_dect_num','_dect_list','_over_num','_over_list','_ave']
    pollution_name_map = {
        'g__Pseudomonas' :'pseud',
        'p__Cyanobacteria' :'cyano',
        'f__Mitochondria': 'mitoc',
        'p__Cyanobacteria_Chloroplast':'cyan_ch'
    }

    desc_map_mongo_file = {}
    #mongo_file_map_desc = {}
    pn = 1
    for p in pollution_list:
        an = 1
        for a in after:
            k = 'field'+str(pn)+str(an)
            desc_map_mongo_file[p+a] = k
            #mongo_file_map_desc[k] = p+a
            an+=1
        pn+=1

    data = pd.read_table(infile,sep='\t',header=0)
    data = data.fillna('-')
    cols = data.columns.tolist()
    map = {
        'pol_num':'_dect_num',
        'pol_samples':'_dect_list',
        'over_num' : '_over_num',
        'over_ave' : '_ave',
        'over_samples' : '_over_list',
        'max_abund' :'_max',
        'max_sample' : '_spe',

    }


    change_name_threshold_dic = {}
    for k in   threshold_dic:
        new_k =  pollution_name_map[k]
        change_name_threshold_dic[new_k] = threshold_dic[k]


    insert_data = {
        'task_id': task_id,
        'project_type':project_type,
        'targets': change_name_threshold_dic,
        'all_num': data['all_num'][0],
        'key_map' : desc_map_mongo_file  #mongo_file_map_desc
    }

    is_pollution = 0
    for i in range(len(data)):
        pollution_name = data['pollution'][i]
        monog_p_pre = pollution_name_map[pollution_name]

        for head in map.keys():
            mongo_k_desc = monog_p_pre+ map[head]
            mongo_k = desc_map_mongo_file[mongo_k_desc]
            insert_data[mongo_k] = data[head][i]
        if  data['is_pol'][i] == 'yes':
            is_pollution = 1
    insert_data['pollute'] =  is_pollution

    if database:
        main_id = database.insert_one(insert_data).inserted_id
        if samples_abunds_file:
            detail_data = import_mongo_detail(samples_abunds_file,main_id,pollution_name_map)
            if database_detail:
                database_detail.insert_many(detail_data)
                # return 'Done'
            else:
                pass
                # return  detail_data
        else:
            # return main_id
            pass

    else:
        #return insert_data
        pass
    return is_pollution


def import_mongo_detail(samples_abunds_file,main_table_id,pollution_name_map):
    samples_info = {}
    abund_data = pd.read_table(samples_abunds_file,sep='\t',header=0)
    samples = abund_data.columns.tolist()[1:]
    pollution_head = abund_data.columns.tolist()[0]
    for s in samples:
        samples_info[s] = {'sample':s,'tas_id':main_table_id}

    for i in range(len(abund_data)):
        pollution_name = abund_data[pollution_head][i]
        k =  pollution_name_map[pollution_name]
        for s in samples_info:
            samples_info[s][k] = abund_data[s][i]

    return samples_info.values()


def check_pollution_pip(summary_dir,task_id,out_dir,is_sanger='tsanger',sanger_mongo=None):

    if sanger_mongo:
        database=pymongo.MongoClient(sanger_mongo)
        db = database['project']
    else:
        db_name = 'project'
        db = Config().get_mongo_client(mtype=db_name)[Config().get_mongo_dbname(db_name)]

    main_collection = db['pollute']
    detail_collection = db['pollute_detail']
    os.system('cat {0}/asv_taxon_Genus.percent.full.xls {0}/asv_taxon_Family.percent.full.xls {0}/asv_taxon_Phylum.percent.full.xls > {1}/g_f_p.xls'.format(summary_dir,out_dir))
    os.system('sed -i "s/.*;//" {0}/g_f_p.xls'.format(out_dir))
    os.system('sed -i "s/f__mitochondria/f__Mitochondria/" {0}/g_f_p.xls'.format(out_dir))  #解决greengene的f__mitochondria是小写。统一成大写
    summary_file = '{0}/g_f_p.xls'.format(out_dir)
    pollution_list = ['g__Pseudomonas','p__Cyanobacteria','f__Mitochondria','p__Cyanobacteria_Chloroplast']
    threshold_dic = {'g__Pseudomonas':0.2,'p__Cyanobacteria':0.3,'f__Mitochondria':0.2,'p__Cyanobacteria_Chloroplast':0.2}
    out_file = '{0}/pol.xls'.format(out_dir)
    pol_res = get_polution(summary_file,'#OTU ID',pollution_list,threshold_dic,outdir=out_dir)
    #task_id = sys.argv[2]
    project_type = 'metaasv.meta_asv'
    sub_abund = out_dir+'/sub_abund.xls'
    is_pollu = import_mongo(out_file,task_id,project_type,threshold_dic,database=main_collection,samples_abunds_file=sub_abund,database_detail=detail_collection)
    return is_pollu


def get_pollution_from_mongo(task_ids,pol_list=None,sanger_mongo=None):
    if not pol_list:
        pol_list = ['g__Pseudomonas','p__Cyanobacteria','f__Mitochondria','p__Cyanobacteria_Chloroplast']
        
    query_map  = {
        'g__Pseudomonas' : {'key':'g__','value':'g__Pseudomonas'},
        'p__Cyanobacteria' : {'key': 'p__', 'value':'p__Cyanobacteria'},
        'f__Mitochondria' : {'key': 'f__', 'value':'f__Mitochondria'},
        'p__Cyanobacteria_Chloroplast' : {'key':'p__', "value":"p__Cyanobacteria_Chloroplast"}
    }

    if sanger_mongo:
        database = pymongo.MongoClient(sanger_mongo)
        db = database['sanger']
    else:
        db_name = 'metaasv'
        db = Config().get_mongo_client(mtype=db_name)[Config().get_mongo_dbname(db_name)]
    otu_main_table = db['asv']
    otu_detail_table = db['asv_detail']
    otu_seq_table = db['asv_specimen']
    task_ids_result_files = {}
    for id in task_ids:
        print id
        try:
            ret_main = otu_main_table.find_one({'task_id': id,'name':'ASV_Origin'})
            main_id = ret_main['_id']

            samples_info = otu_seq_table.find({"asv_id":main_id})
            samples_reads = {}
            if samples_info:
                for sample_i in  samples_info:
                    specimen_name = sample_i['specimen_name']
                    samples_reads[specimen_name] = sample_i['read_number']

            sample_pol_abund = {}
            for  p in pol_list:
                q = query_map[p]['key']
                q_v = query_map[p]['value']

                detail_rets = otu_detail_table.find({"otu_id":main_id,q:q_v})
                if p == 'f__Mitochondria':
                    if not detail_rets:  # 有的库的名字是f__Mitochondria，少数库是f__mitochondria
                        q_v = 'f__mitochondria'
                        detail_rets = otu_detail_table.find({"otu_id":main_id,q:q_v})

                if detail_rets:
                    sample_pol_abund[p] = {}
                    for sample in  samples_reads.keys():
                        sample_pol_abund[p][sample] = 0


                    for detail_ret in detail_rets:
                        for sample in samples_reads.keys():
                            sample_pol_abund[p][sample] += detail_ret[sample]


                    for sample in samples_reads.keys():

                        sample_pol_abund[p][sample] = float(sample_pol_abund[p][sample])/samples_reads[sample]

            task_ids_result_files[id] = str(id)+'_pollution.xls'
            with open(str(id)+'_pollution.xls','w') as fw:
                sample_list = sorted(samples_reads.keys())
                fw.write('#ASV ID\t'+'\t'.join(sample_list)+'\n')
                for p in  sample_pol_abund.keys():
                    fw.write(p)
                    for sample  in sample_list:
                        fw.write('\t'+str(sample_pol_abund[p][sample]))
                    fw.write('\n')
        except Exception  as e:
            print 'NOT DONE %s' % id
            print e

    return task_ids_result_files


def get_task_id(start_time=None,end_time=None,sanger_mongo=None):
    #return set(['i-sanger_201722','sanger_175247','i-sanger_178128'])
    #return ['i-sanger_202560']
    if sanger_mongo:
        database = pymongo.MongoClient(sanger_mongo)
        db = database['sanger']
    else:
        db_name = 'meta'
        db = Config().get_mongo_client(mtype=db_name)[Config().get_mongo_dbname(db_name)]
    task_table = db['sg_task']
    all_task = task_table.find()
    task_ids = []

    if start_time:
        print start_time
        start_time = datetime.datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
    if end_time:
        print end_time
        end_time = datetime.datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")

    for task in all_task:
        if 'task_id' not in task.keys():
            continue

        if  not start_time and not end_time:
            task_ids.append(task['task_id'])

        else:
            if 'created_ts' not in task.keys():
                continue
            ts = datetime.datetime.strptime(task['created_ts'], "%Y-%m-%d %H:%M:%S")
            if start_time and not end_time:
                if ts > start_time:
                    task_ids.append(task['task_id'])
            elif end_time and not start_time:
                if ts < end_time:
                    task_ids.append(task['task_id'])
            elif start_time and end_time:
                if ts > start_time and  ts < end_time:
                    task_ids.append(task['task_id'])

    return set(task_ids)

def update_old_project_mongo_pip(out_dir=None,start_time=None,end_time=None,is_sanger='tsanger',sanger_mongo=None,task_mongo=None):
    if not out_dir:
        out_dir = './tmp_mongo_dir'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    if sanger_mongo:
        database=pymongo.MongoClient(sanger_mongo)
        db = database['project']
    else:
        db_name = 'project'
        db = Config().get_mongo_client(mtype=db_name)[Config().get_mongo_dbname(db_name)]
    main_collection = db['pollute']
    detail_collection = db['pollute_detail']
    pollution_list = ['g__Pseudomonas','p__Cyanobacteria','f__Mitochondria','p__Cyanobacteria_Chloroplast']
    threshold_dic = {'g__Pseudomonas':0.2,'p__Cyanobacteria':0.3,'f__Mitochondria':0.2,'p__Cyanobacteria_Chloroplast':0.2}
    project_type = 'meta.meta_base'

    task_ids = get_task_id(start_time,end_time,sanger_mongo=task_mongo)
    #print task_ids
    print len(task_ids)


    pol_files = get_pollution_from_mongo(task_ids,sanger_mongo=task_mongo)

    for id in pol_files:
        summary_file = pol_files[id]
        out_old_file = '{0}/pol.xls'.format(out_dir)
        out_file = '{0}/{1}_pol.xls'.format(out_dir,id)
        pol_res = get_polution(summary_file,'#OTU ID',pollution_list,threshold_dic,outdir=out_dir)
        os.system('mv %s %s' %(out_old_file,out_file))
        sub_abund = summary_file
        import_mongo(out_file,id,project_type,threshold_dic,database=main_collection,samples_abunds_file=sub_abund,database_detail=detail_collection)
