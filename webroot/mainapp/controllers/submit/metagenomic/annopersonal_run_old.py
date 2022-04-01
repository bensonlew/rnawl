# -*- coding: utf-8 -*-
#!/usr/bin/python
#__author__ = 'guanqing.zou'


import os
import json
import web
import json
import datetime

from mainapp.controllers.project.metagenomic_controller import MetagenomicController

class AnnopersonalRunAction(MetagenomicController):
    def __init__(self):
        super(AnnopersonalRunAction,self).__init__(instant=False)

    def POST(self):
        data = web.input()
        print(data)
        default_argu = ['task_id']
        for argu in default_argu:
            if not hasattr(data, argu):
                info = {'success':False, "info" : "parameters %s is missing"}
                return json.dumps(info)
        task_name = "metagenomic.report.personal_anno"
        geneset_table = self.metagenomic.find_origin_geneset_info(data.task_id)
        geneset_id = str(geneset_table['_id'])
        gene_length = geneset_table['gene_list_length']
        reads_num_file =geneset_table['reads_num']
        faa_file = reads_num_file.split('gene_profile')[0] + 'uniGeneset/gene.uniGeneset.faa'
        pre_summary = reads_num_file.split('gene_profile')[0] + 'anno_overview.xls'
        find_table = self.metagenomic.common_find_one('anno_overview',{'task_id': data.task_id})
        if find_table:
            if 'new_sum' in find_table.keys():
                pre_summary = find_table['new_sum']

        project_sn = geneset_table["project_sn"]
        database_list = []
        if hasattr(data, 'database'):
            database_list = data.database.split(',')   #go,phi,mvirdb,tcdb,qs,pfam,secretory_protein,t3ss,probio,p450
        nr_list = []
        if hasattr(data,'nr'):
            nr_list = data.nr.split(',')

        self.run_list = []
        self.do_api_list = []
        self.main_id_list = []
        self.down_files_path = {}
        self.done_list = []
        for i in database_list:
            table_suffix = i
            anno_type = i
            if i == 'mvirdb': table_suffix = 'mvir'
            if i == 'p450':
                table_suffix = 'cyps'
                anno_type = 'cyps'
            if i == 'secretory_protein' :
                table_suffix = 'sec'
                anno_type = 'sec'
                i ='sec'
            if i == 't3ss' : table_suffix = 'ttss'
            result = self.metagenomic.common_find_one('anno_personal_run',{'task_id':data.task_id, 'anno_type': anno_type})
            if result:
                if result['status'] == 'start' or result['status'] == 'end' :
                    continue

            if i == 'probio':
                for t in  nr_list + ['best_hit']:
                    find_table = ''
                    data_name = 'probio'
                    if t == 'lca':
                        data_name = 'probio_lca'
                        find_table = self.metagenomic.common_find_one('anno_probio',{'task_id': data.task_id,'is_origin':1, 'name':'Probio_Origin_LCA'})
                    elif t == 'de_unclassified':
                        data_name = 'probio_de_unclassified'
                        find_table = self.metagenomic.common_find_one('anno_probio',{'task_id': data.task_id,'is_origin':1, 'name':'Probio_Origin_Deunclassified'})
                    else:
                        find_table = self.metagenomic.common_find_one('anno_probio',{'task_id': data.task_id,'is_origin':1, 'name':'Probio_Origin'})
                    if find_table:
                        if find_table['status'] == 'hide':
                            self.main_id_list.append(str(find_table['_id']))
                            self.do_api_list.append(data_name)
                            anno_path = find_table['anno_file'].split('/')
                            anno_path.pop()
                            anno_dir = '/'.join(anno_path) + '/'
                            #self.down_files_path.append(anno_dir)
                            if 'personal_anno' in anno_dir:   ###临时
                                anno_dir = anno_dir.replace('personal_anno', 'anno_personal')
                            self.down_files_path[data_name +'_file'] = anno_dir

                        elif find_table['status'] in ['failed','Failed']:
                            self.run_list.append(data_name)
                        elif find_table['status'] in ['end','End','END'] :
                            self.done_list.append(data_name)
                    else:
                        self.run_list.append(data_name)
            elif i == 't3ss':
                find_table = self.metagenomic.common_find_one('anno_ttss',{'task_id': data.task_id,'is_origin':1, 'name':'TTSS_Origin'})
                if find_table:
                    main_id = find_table['_id']
                    if find_table['status'] == 'hide':
                        self.main_id_list.append(str(find_table['_id']))
                        self.do_api_list.append('t3ss')
                        anno_path = find_table['anno_file'].split('/')
                        anno_path.pop()
                        anno_dir = '/'.join(anno_path) + '/'
                        self.down_files_path['t3ss_file'] = anno_dir
                    else:
                        for t in  nr_list + ['best_hit']:
                            find_table2 = self.metagenomic.common_find_one('anno_ttss_stat',{'ttss_id':main_id,"type":t})
                            if not find_table2:
                                if i not in self.run_list:
                                    self.run_list.append(i)
                else:
                    self.run_list.append(i)

            else:
                find_table = self.metagenomic.find_anno_main_table('anno_'+ table_suffix, data.task_id)
                if find_table:
                    #if True:
                    if find_table['status'] == 'hide':
                        self.main_id_list.append(str(find_table['_id']))
                        self.do_api_list.append(i)
                        anno_path = find_table['anno_file'].split('/')
                        anno_path.pop()
                        anno_dir = '/'.join(anno_path) + '/'
                        #self.down_files_path.append(anno_dir)
                        if 'personal_anno' in anno_dir:   ###临时
                            anno_dir = anno_dir.replace('personal_anno', 'anno_personal')
                        self.down_files_path[i+'_file'] = anno_dir
                    elif find_table['status'] in ['failed','Failed']:
                        self.run_list.append(i)
                    elif find_table['status'] in ['end','End','END'] :
                        self.done_list.append(i)
                    else:
                        self.run_list.append(i)

                else:
                    self.run_list.append(i)

        geneset = self.metagenomic.common_find_one('geneset', {'task_id': data.task_id,'name':'GENESET_Origin', "type" : 1})
        if geneset:
            specimen = geneset['specimen'].split(',')
        else:
            info = {'success':False, "info" : "specimen info  not Found"}
            return json.dumps(info)

        options = {
            "task_id" : data.task_id,
            "geneset_id" : geneset_id ,
            "group_detail" : json.dumps({'all':specimen}),
            "gene_length_table" : gene_length ,
            "specimen_group" : "all",
            "pre_summary" : pre_summary,
            #"summary_path" : pre_summary
        }
        do_nr = []
        nr_only_api = []
        nr_main_id_list = []
        if hasattr(data, 'nr'):
            if 'lca' in nr_list:        #de_unclassified,lca
                r =self.metagenomic.common_find_one('anno_nr', {'task_id':data.task_id,'name':'NR_Origin_LCA'})
                if not r:
                    personal_run_find = self.metagenomic.common_find_one('anno_personal_run',{'task_id':data.task_id, 'anno_type': 'lca'})
                    if personal_run_find:
                        if personal_run_find['status'] =='start' or personal_run_find['status'] =='end':
                             pass
                        else:
                            do_nr.append('lca')
                    else:
                         do_nr.append('lca')
                else:
                    #if True:
                    if r['status'] == 'hide':
                        nr_only_api.append('lca')
                        nr_main_id_list.append(str(r['_id']))
                        options['nr_lca_dir'] = '/'.join(r['anno_file'].split('/')[:-1]) + '/'
                    options['nr_gene_anno_lca'] = r['anno_file']
            if 'de_unclassified' in nr_list :
                r =self.metagenomic.common_find_one('anno_nr', {'task_id':data.task_id,'name':'NR_Origin_Deunclassified'})
                if not r:
                    personal_run_find = self.metagenomic.common_find_one('anno_personal_run',{'task_id':data.task_id, 'anno_type': 'de_unclassified'})
                    if personal_run_find:
                        if personal_run_find['status'] =='start' or personal_run_find['status'] =='end':
                            pass
                        else:
                            do_nr.append('de_unclassified')
                    else:
                        do_nr.append('de_unclassified')
                else:
                    #if True:
                    if r['status'] == 'hide':
                        nr_only_api.append('de_unclassified')
                        nr_main_id_list.append(str(r['_id']))
                        options['nr_de_dir'] = '/'.join(r['anno_file'].split('/')[:-1]) + '/'
                    options['nr_gene_anno_de'] = r['anno_file']

            options['nr_method'] = ','.join(do_nr)
            options['nr_method_ori'] = data.nr
            options['nr_only_api'] = ','.join(nr_only_api)
            options['nr_main_id_list'] = ','.join(nr_main_id_list)


        to_file = []
        if len(self.run_list)!=0 or len(do_nr) != 0 :
            options['run_list'] = ";".join(self.run_list)
            options["reads_profile_table"] = reads_num_file
            options["query"] = faa_file
            if 'qs' in self.run_list:
                options['group_table'] = geneset_id
                to_file.append('metagenomic.export_group_table_by_detail_2(group_table)')
            if 'go' in self.run_list:
                nr_table = self.metagenomic.common_find_one('anno_nr', {'task_id':data.task_id,'name':'NR_Origin'})
                options['blastout'] = '/'.join(nr_table['anno_file'].split('/')[:-1]) + '/nr_align_table.xls'



        if len(self.do_api_list) != 0:
            options['do_api_list'] = ";".join(self.do_api_list)
            for k in self.down_files_path.keys():
                options[k] = self.down_files_path[k]
            options["main_id_list"]  = ";".join(self.main_id_list)
        if 'probio' in self.run_list or 'sec' in self.run_list or 't3ss' in self.run_list:
            nr_table = self.metagenomic.common_find_one('anno_nr', {'task_id':data.task_id,'name':'NR_Origin'})
            if nr_table:
                options['nr_gene_anno'] = nr_table['anno_file']
            if 'lca' in nr_list:
                nr_lca_table = self.metagenomic.common_find_one('anno_nr', {'task_id':data.task_id,'name':'NR_Origin_LCA'})
                if nr_lca_table:
                    options['nr_gene_anno_lca'] = nr_lca_table['anno_file']
            if 'de_unclassified' in nr_list:
                nr_de_table = self.metagenomic.common_find_one('anno_nr', {'task_id':data.task_id,'name':'NR_Origin_Deunclassified'})
                if nr_de_table:
                    options['nr_gene_anno_de'] = nr_de_table['anno_file']


        for i in self.run_list:
            if i == 'mvirdb': i = 'mvir'
            if i == 'p450': i = 'cyps'
            submit = i
            anno_type = i
            name = i.capitalize()+'_Origin'
            if i == 't3ss':
                submit = 'ttss'
                name = 'T3ss_Origin'
            elif i == 'probio_lca':
                submit = 'probio'
                name = 'Probio_Origin_LCA'
                anno_type = 'probio'
            elif i == 'probio_de_unclassified':
                submit = 'probio'
                name = 'Probio_Origin_Deunclassified'
                anno_type = 'probio'
            result = self.metagenomic.common_find_one('anno_personal_run',{'task_id':data.task_id, 'name': name})
            if result:
                self.metagenomic.common_update_one('anno_personal_run',result['_id'],{'status':'start','desc':'任务在运行'})
                continue
            mongo_data = [('anno_type',anno_type),('desc','任务在运行'),
                          ('status','start'),('submit_location','anno'+submit),
                          ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                          ('name',name),('task_id',data.task_id)]
            self.metagenomic.insert_main_table("anno_personal_run", mongo_data)

        for i in self.do_api_list:
            if i == 'mvirdb': i = 'mvir'
            if i == 'p450': i = 'cyps'
            submit = i
            anno_type = i
            name = i.capitalize()+'_Origin'
            if i == 't3ss':
                submit = 'ttss'
                name = 'T3ss_Origin'
            elif i == 'probio_lca':
                submit = 'probio'
                name = 'Probio_Origin_LCA'
                anno_type = 'probio'
            elif i == 'probio_de_unclassified':
                submit = 'probio'
                name = 'Probio_Origin_Deunclassified'
                anno_type = 'probio'
            result = self.metagenomic.common_find_one('anno_personal_run',{'task_id':data.task_id, 'name': name})
            if result:
                self.metagenomic.common_update_one('anno_personal_run',result['_id'],{'status':'start','desc':'任务在运行'})
                continue
            mongo_data = [('anno_type',anno_type),('desc','任务在运行'),
                          ('status','start'),('submit_location','anno'+submit),
                          ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                          ('name',name),('task_id',data.task_id)]
            self.metagenomic.insert_main_table("anno_personal_run", mongo_data)

        for i in self.done_list:
            if i == 'mvirdb': i = 'mvir'
            if i == 'p450': i = 'cyps'
            submit = i
            anno_type = i
            name = i.capitalize()+'_Origin'
            if i == 't3ss':
                submit = 'ttss'
                name = 'T3ss_Origin'
            elif i == 'probio_lca':
                submit = 'probio'
                name = 'Probio_Origin_LCA'
                anno_type = 'probio'
            elif i == 'probio_de_unclassified':
                submit = 'probio'
                name = 'Probio_Origin_Deunclassified'
                anno_type = 'probio'
            result = self.metagenomic.common_find_one('anno_personal_run',{'task_id':data.task_id, 'name': name})
            if result:
                self.metagenomic.common_update_one('anno_personal_run', result['_id'], {'status':'end','desc':'任务结束'})
                continue
            mongo_data = [('anno_type',anno_type),('desc','任务结束'),
                          ('status','end'),('submit_location','anno'+submit),
                          ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                          ('name',name),('task_id',data.task_id)]
            self.metagenomic.insert_main_table("anno_personal_run", mongo_data)

        for i in nr_list:
            result = self.metagenomic.common_find_one('anno_personal_run',{'task_id':data.task_id, 'anno_type': i})
            if result:
                if i not in do_nr:
                    self.metagenomic.common_update_one('anno_personal_run', result['_id'], {'status':'end','desc':'任务结束'})
                continue
            mongo_data = [('anno_type',i),('desc','任务在运行'),
                          ('status','start'),('submit_location','annonr'),
                          ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                          ('name',i.capitalize()+'_Origin'),('task_id',data.task_id)]
            self.metagenomic.insert_main_table("anno_personal_run", mongo_data)

        if len(self.run_list) == 0 and len(self.do_api_list) == 0 and len(do_nr) == 0  and len(nr_only_api) == 0 :
            info = {'success':False, "info" : "RUN END!"}
            return json.dumps(info)



        self.set_sheet_data(name=task_name, options=options, main_table_name= 'Personal_Anno',
                            module_type='workflow', project_sn=project_sn, task_id=data.task_id, to_file=to_file)
        task_info = super(AnnopersonalRunAction, self).POST()
        #if task_info['success']:
        #    task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

















