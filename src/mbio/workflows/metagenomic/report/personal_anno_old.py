# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifiy =  2081116

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mainapp.models.mongo.metagenomic import Metagenomic
import os
import re
import shutil
import json
import types
from biocluster.api.file.lib.transfer import MultiFileTransfer
from mainapp.libs.param_pack import group_detail_sort
import copy
from biocluster.core.exceptions import OptionError

class PersonalAnnoWorkflow(Workflow):
    """
    宏基因组个性化注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PersonalAnnoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name":"do_api_list","type":"string","default":""},   # 需导表的分析项，不包括先运行注释再导表的分析
            {"name":"main_id_list","type":"string"},
            #{"name":"down_files_path","type":"string"},
            {"name": "run_list","type":"string"},
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "geneset_id", "type": "string", "default":""},
            {"name": "group_detail","type": "string"},
            {"name": "nr_gene_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "nr_gene_anno_lca", "type": "infile", "format": "sequence.profile_table"},
            {"name": "nr_gene_anno_de", "type": "infile", "format": "sequence.profile_table"},
            {"name": "nr_method", "type" : "string", "default" : ""} ,  # deunclassied  or lca
            {"name": "nr_method_ori", "type" : "string", "default" : ""} ,
            {"name": "specimen_group", "type" : "string", "default" : ""},
            {"name": "gene_length_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "blastout", "type": "infile", "format": "sequence.profile_table"},
            {"name": "task_id", "type":"string"},
            {"name": "pre_summary", "type": "infile", "format": "sequence.profile_table"},
            #{"name": "summary_path", "type": "string" },
            # {"name": "tcdb_file","type": "infile", "format": "sequence.profile_table"},
            # {"name": "mvirdb_file","type": "infile", "format": "sequence.profile_table"},
            # {"name": "qs_file","type": "infile", "format": "sequence.profile_table"},
            # {"name": "probio_file","type": "infile", "format": "sequence.profile_table"},
            # {"name": "phi_file","type": "infile", "format": "sequence.profile_table"},
            # {"name": "pfam_file","type": "infile", "format": "sequence.profile_table"},
            # {"name": "go_file","type": "infile", "format": "sequence.profile_table"},
            {"name": "tcdb_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "mvirdb_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "qs_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "probio_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "probio_lca_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "probio_de_unclassified_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "phi_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "pfam_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "go_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "sec_file","type": "infile", "format": "meta_genomic.sec_dir"},
            {"name": "p450_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "t3ss_file","type": "infile", "format": "meta_genomic.simple_dir"},
            {"name": "nr_only_api","type": "string", 'default':""},
            {"name": "nr_main_id_list", "type":"string", "default":""},
            {"name": "nr_lca_dir", "type": "infile", "format": "meta_genomic.nr_dir"},
            {"name": "nr_de_dir", "type": "infile", "format": "meta_genomic.nr_dir"},

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())


        self.remote_dir = self._sheet.output
        self.geneset_id = self.option('geneset_id')

        self.group_detail = eval(self.option('group_detail'))
        self.specimen_group = self.option('specimen_group')


        self.need_run_list= []
        self.do_api_list = []
        self.nr_run = []
        self.nr_only_api = []
        self.nr_main_id = []
        if self.option('run_list'):
            self.need_run_list = self.option('run_list').split(';')

        if self.option('do_api_list'):
            self.do_api_list = self.option('do_api_list').split(';')
            self.main_id_list = self.option('main_id_list').split(';')
            #self.down_files_list =  self.option('down_files_path').split(';')
        if self.option('nr_method'):
            self.nr_run = self.option('nr_method').split(',')
        if self.option('nr_only_api'):
            self.nr_only_api = self.option('nr_only_api').split(',')
            self.nr_main_id = self.option('nr_main_id_list').split(',')

        self.anno_nr_lca= self.add_module('annotation.mg_common_anno_stat')
        self.anno_nr_deunclass= self.add_module('annotation.mg_common_anno_stat')
        self.anno_module = self.add_module('annotation.personal_anno')
        self.diamond_nr = self.add_module('align.meta_diamond')
        self.overview = self.add_tool('annotation.database_merge')
        self.metagenomic = Metagenomic()
        self.metagenomic._config = Config()

        self.need_run_result = {}
        self.anno_dic = {'mvirdb':'gene_mvirdb_anno.xls','tcdb':'gene_tcdb_anno.xls',
                    'probio' : 'gene_probio_anno.xls','pfam':'gene_pfam_anno.xls',
                    'p450':'gene_cyps_anno.xls','t3ss':'ttss_predict.txt',
                    'phi':'gene_phi_anno.xls','go':'all.go.annotation.xls',
                    'qs': 'gene_qs_anno.xls','probio_lca':'gene_probio_anno.xls',
                    'probio_de_unclassified':'gene_probio_anno.xls'
                    }

        if len(self.need_run_list) != 0 or self.option('nr_method') != '' :
            self.reads_profile_table = self.option('reads_profile_table').prop['path']
            self.query = self.option('query').prop['path']

        for k in self.need_run_list:
            target_out = self.output_dir + '/' + k.capitalize()
            if k == 't3ss':
                target_out = self.output_dir + '/' + 'ttss'.capitalize()
            elif k == 'probio_lca':
                target_out = self.output_dir + '/Probio_LCA/'
            elif k == 'probio_de_unclassified':
                target_out = self.output_dir + '/Probio_Deunclassifed/'
            if k not in self.need_run_result.keys():
                self.need_run_result[k] = {}
            self.need_run_result[k]['profile_dir'] = target_out


        if 'probio' in self.need_run_list or 't3ss' in self.need_run_list:
            if not self.option("nr_gene_anno").is_set and not self.option("nr_gene_anno_lca") and not self.option("nr_gene_anno_de"):
                raise OptionError("annotation must input nr_gene_anno file!", code="12803801")
        if 'qs' in self.need_run_list and not self.option("group_table").is_set:
            raise OptionError("must input group_table", code="12803802")
        if 'go' in self.need_run_list and not self.option("blastout").is_set:
            raise OptionError("go annotation must input blastout file!", code="12803803")

        self.database_version = {
            "mvir" : "mvirDB_v20110724",
            "cyps" : "CYPED_v6.0",
            # "probio" : "probio_v2016.06.20",
            "phi": "phi_v4.9",
            "tcdb": "tcdb_v20200917",
            "sec": "v4.1",
            "t3ss": "v1.0.1",
            "pfam": "pfam_v33.1",
            "go": "go_v1.2",
            "lca": "nr_v20200604",
            "de_unclassified": "nr_v20200604",
            #"qs" : "-" # 没有版本号
        }

    def run(self):

        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("Start Run Workflow")

        if self._sheet.id in ['i-sanger_183577_44543_128356']:
            # 导表以完成，直接上传文件
            import gevent
            self.IMPORT_REPORT_DATA = False
            gevent.spawn_later(5, self.end)
            super(PersonalAnnoWorkflow, self).run()
            return

        self.on_list = []
        self.first_module = [self.run_dimond_nr]
        if self.option('nr_method') == 'de_unclassified':
            self.diamond_nr.on('end',self.run_anno_nr_deunclass)
            self.on_list.append(self.anno_nr_deunclass)
        elif self.option('nr_method') == 'lca':
            self.diamond_nr.on('end',self.run_anno_nr_lca)
            self.on_list.append(self.anno_nr_lca)
        elif self.option('nr_method') == 'lca,de_unclassified' or self.option('nr_method') == 'de_unclassified,lca':
            self.diamond_nr.on('end',self.run_anno_nr_lca)
            self.on_list.append(self.anno_nr_lca)
            self.diamond_nr.on('end',self.run_anno_nr_deunclass)
            self.on_list.append(self.anno_nr_deunclass)
        else:
            self.first_module = []

        if len(self.on_list) > 1:
            if len(self.need_run_list) != 0:
                self.on_rely(self.on_list, self.run_anno_module)
                self.anno_module.on('end', self.run_overview)
                self.overview.on('end', self.set_db)
            else:
                self.on_rely(self.on_list, self.run_overview)
                self.overview.on('end', self.set_db)

        elif len(self.on_list) == 1:
            if len(self.need_run_list) != 0:
                self.on_list[0].on('end', self.run_anno_module)
                self.anno_module.on('end', self.run_overview)
                self.overview.on('end', self.set_db)
            else:
                self.on_list[0].on('end', self.run_overview)
                self.overview.on('end', self.set_db)
        elif len(self.on_list) == 0:
            if len(self.need_run_list) != 0:
                self.first_module = [self.run_anno_module]
                self.anno_module.on('end', self.run_overview)
                self.overview.on('end', self.set_db)
            else:
                self.first_module = [self.run_overview]
                self.overview.on('end',self.set_db)


        if len(self.do_api_list) != 0 :
            self.run_down_api_file(self.first_module)

        else:
            for func in self.first_module:
                func()

        super(PersonalAnnoWorkflow, self).run()



    def run_down_api_file(self, next_functions):
        self.api_main_dic = {}
        self.api_files_dic = {}
        for k in self.do_api_list:
            i = self.do_api_list.index(k)
            tar_file = self.option(k + '_file').prop['path']
            # tar_file = self.down_files_list[i]
            # if re.search(r'://', tar_file):
            #     transfer = MultiFileTransfer()
            #     if tar_file[-1] != '/':
            #         tar_file += '/'
            #     transfer.add_download(tar_file ,self.work_dir + '/Down2Api/' + k + '/')
            #     transfer.perform()
            #     tar_file = self.work_dir + '/Down2Api/' + k + '/'
            main_id = self.main_id_list[i]
            self.api_main_dic[k] = main_id
            if k not in self.api_files_dic.keys():
                self.api_files_dic[k] = {}
            #tar = '/'.join(tar_file.split('/')[:-1]) + '/'
            self.api_files_dic[k]['profile_dir'] = tar_file

        for func in next_functions:
            func()

    def run_dimond_nr(self):
        opts = {
            'query': self.option('query'),
            'query_type': "prot",
            'database': 'nr',
            'lines': '50000',
            "target_num": 5,
        }
        self.diamond_nr.set_options(opts)
        self.diamond_nr.run()

    def run_anno_nr_lca(self):
        opts = {
            'reads_profile_table':  self.option("reads_profile_table"),
            "nr_method": "lca",
            "nr_xml_dir": self.diamond_nr.option('outxml_dir'),
            "out_type" :1   #zouguanqing 20190319
        }

        self.anno_nr_lca.set_options(opts)
        self.anno_nr_lca.run()

    def run_anno_nr_deunclass(self):
        opts = {
            'reads_profile_table':  self.option("reads_profile_table"),
            "nr_method": "deunclassied",
            "nr_xml_dir": self.diamond_nr.option('outxml_dir'),
            "out_type" :1  #zouguangqing 20190319
        }

        self.anno_nr_deunclass.set_options(opts)
        self.anno_nr_deunclass.run()


    def run_anno_module(self):
        part = copy.deepcopy(self.need_run_list)
        mk_rm = 0
        for d in ['probio_lca', 'probio_de_unclassified']:
            if d in self.need_run_list:
                mk_rm =1
                part.remove(d)
        if mk_rm != 0:
            part.append('probio')
        database = ';'.join(part)
        set_dic = {
            "database": database,
            "query" : self.option('query'),
            "reads_profile_table" : self.option("reads_profile_table"),
        }
        if self.option('nr_gene_anno').is_set:
            set_dic['nr_gene_anno'] = self.option('nr_gene_anno')
        if 'de_unclassified' in self.option('nr_method_ori').split(','):
            if 'de_unclassified' in self.option('nr_method').split(','):
                set_dic['nr_gene_anno_de'] = self.anno_nr_deunclass.option('tax_level_dir').prop['path'] +'/gene_nr_anno.xls'
            else:
                set_dic['nr_gene_anno_de'] = self.option('nr_gene_anno_de')
        if 'lca' in  self.option('nr_method_ori').split(','):
            if 'lca' in  self.option('nr_method').split(','):
                set_dic['nr_gene_anno_lca'] = self.anno_nr_lca.option('tax_level_dir').prop['path'] + '/gene_nr_anno.xls'
            else:
                set_dic['nr_gene_anno_lca'] = self.option('nr_gene_anno_lca')
        if 'qs' in self.option('run_list'):
            set_dic['group_table'] = self.option('group_table')
        if 'go' in self.need_run_list:
            set_dic['blastout'] = self.option('blastout')


        self.anno_module.set_options(set_dic)
        self.anno_module.run()

    def run_overview(self):
        opts = {
            "gene_length_table" : self.option('gene_length_table')
        }
        types = []
        files = []
        for anno in self.need_run_list:
            if anno == 'sec':
                for j in ['sec_Gram_neg','sec_Gram_pos', 'sec_Euk']:
                    types.append(j)
                for j in ['signalp_Gram-_SignalP.txt','signalp_Gram+_SignalP.txt','signalp_Euk_SignalP.txt']:
                    files.append(self.anno_module.output_dir + '/'+ anno.capitalize() + '/' + j)
            elif anno in ['probio_lca', 'probio_de_unclassified']:
                continue
            else:
                types.append(anno)
                if anno == 't3ss':
                    files.append(self.anno_module.output_dir + '/Ttss/' + self.anno_dic[anno])
                else:
                    files.append(self.anno_module.output_dir + '/'+ anno.capitalize() + '/' + self.anno_dic[anno])


        for anno in self.do_api_list:
            if anno == 'sec':
                for j in ['sec_Gram_neg','sec_Gram_pos', 'sec_Euk']:
                    types.append(j)
                for j in ['signalp_Gram-_SignalP.txt','signalp_Gram+_SignalP.txt','signalp_Euk_SignalP.txt']:
                    files.append(self.api_files_dic[anno]['profile_dir'] + '/' + j)
            elif anno == 'probio_lca' or anno == 'probio_de_unclassified':
                continue
            else:
                types.append(anno)
                files.append(self.api_files_dic[anno]['profile_dir']  + '/' + self.anno_dic[anno])

        for anno in self.nr_only_api:
            if anno == 'lca':
                types.append('nr_lca')
                files.append(self.option('nr_gene_anno_lca').prop['path'])
            elif anno == 'de_unclassified':
                types.append('nr_de_unclassified')
                files.append(self.option('nr_gene_anno_de').prop['path'])


        for anno in self.nr_run:
            if anno == 'lca':
                types.append('nr_lca')
                files.append(self.anno_nr_lca.option('tax_level_dir').prop['path'] +'/gene_nr_anno.xls')
            elif anno == 'de_unclassified':
                types.append('nr_de_unclassified')
                files.append(self.anno_nr_deunclass.option('tax_level_dir').prop['path'] +'/gene_nr_anno.xls')


        opts['other'] = ','.join(types)
        opts['otherf'] = ','.join(files)
        if  opts['other'] == '':
            opts['other'] = 'None'
            opts['otherf'] = 'None'
        if self.option('pre_summary').is_set:
            opts['pre_sum'] = self.option('pre_summary')
        self.overview.set_options(opts)
        self.overview.run()



    def set_db(self):
        """
        保存结果output，导mongo数据库
        """

        for k in self.need_run_result.keys():
            target_out = self.need_run_result[k]['profile_dir']
            if os.path.exists(target_out):
                shutil.rmtree(target_out)
            os.mkdir(target_out)
            dir_name = k.capitalize()
            if k == 'probio_lca':
                dir_name = 'Probio_LCA'
            elif k == 'probio_de_unclassified':
                dir_name = 'Probio_Deunclassifed'
            elif k == 't3ss':
                dir_name = 'Ttss'
            for i in os.listdir(self.anno_module.output_dir + '/'+ dir_name):
                basename = os.path.basename(i)
                ori_path = self.anno_module.output_dir + '/'+ dir_name+ '/' + i
                self.logger.info(ori_path)
                self.logger.info(target_out+'/'+ basename)
                os.link(ori_path, target_out+'/'+ basename)


        for k in self.nr_run:
            if k == 'lca':
                target_out = self.output_dir + '/nr_lca'
                if not os.path.exists(target_out):
                    os.mkdir(target_out)
                tax_level_dir = self.anno_nr_lca.option('tax_level_dir').prop['path']
                for f in os.listdir(tax_level_dir):
                    ori_path = os.path.join(tax_level_dir, f)
                    if os.path.exists(target_out + '/' + f):
                        os.remove(target_out + '/' + f)
                    os.link(ori_path, target_out + '/' + f)
            elif k == 'de_unclassified':
                target_out = self.output_dir + '/nr_deunclassied'
                if not os.path.exists(target_out):
                    os.mkdir(target_out)
                tax_level_dir = self.anno_nr_deunclass.option('tax_level_dir').prop['path']
                for f in os.listdir(tax_level_dir):
                    ori_path = os.path.join(tax_level_dir, f)
                    if os.path.exists(target_out + '/' + f):
                        os.remove(target_out + '/' + f)
                    os.link(ori_path, target_out + '/' + f)


        for i in ['gene_overview_anno.xls']:
            overview_file = self.overview.output_dir + '/' + i
            if os.path.exists(overview_file):
                target = self.output_dir + '/' + i
                if os.path.exists(target):
                    os.remove(target)
                os.link(overview_file, target)


        self.logger.info("正在写入mongo数据库")

        common_param = {}
        self.personal_api = {}
        common_param["geneset_id"] = self.geneset_id
        common_param["group_detail"] = self.group_detail
        common_param["anno_dir"] = self.remote_dir #os.path.join(self.remote_dir, "anno_personal")
        self.personal_api_anno = self.api.api("metagenomic.personal_anno")
        self.personal_api["go"] = self.api.api("metagenomic.anno_go")
        self.personal_api["qs"] = self.api.api("metagenomic.qs_anno")
        self.personal_api["probio"] = self.api.api("metagenomic.probiotics")
        self.personal_api["pfam"] = self.api.api("metagenomic.mg_anno_pfam")
        self.personal_api["p450"] = self.api.api("metagenomic.mg_anno_cyps")
        self.personal_api["tcdb"] = self.api.api("metagenomic.tcdb")
        self.personal_api["mvirdb"] = self.api.api("metagenomic.mvirdb")
        self.personal_api["sec"] = self.api.api("metagenomic.mg_anno_sec")
        self.personal_api["t3ss"]  = self.api.api("metagenomic.mg_anno_ttss")
        self.personal_api["phi"] = self.api.api("metagenomic.mg_anno_phi")
        self.nr_api = self.api.api("metagenomic.mg_anno_nr")
        self.overview_api = self.api.api('metagenomic.mg_anno_overview')
        self.personal_params = {
            "geneset_id": str(self.geneset_id),
            "group_detail": group_detail_sort(self.group_detail),
            "group_id": str(self.specimen_group),
            "database": "",
            "submit_location": "",
            "task_type": 2
        }
        if len(self.need_run_list) != 0:
            part = copy.deepcopy(self.need_run_list)
            method_dic = {'probio_lca':'lca','probio_de_unclassified':"de_unclassied"}
            for d in ['probio_lca', 'probio_de_unclassified']:
                if d in self.need_run_list:
                    part.remove(d)
                    probio_dir = {'probio':self.need_run_result[d]}
                    self.personal_api_anno.add_main('probio', self.personal_api, self.personal_params, common_param,
                                                    task_id=self.option('task_id'),nr_method=method_dic[d])
                    self.logger.info(d +' anno_probio _id: ' + str(self.personal_api_anno.main_dic))
                    self.personal_api_anno.add_all_detail('probio', self.personal_api, self.personal_api_anno.main_dic, probio_dir)
            if part != []:
                database = ';'.join(part)
                self.personal_api_anno.add_main(database, self.personal_api, self.personal_params, common_param,task_id=self.option('task_id'))
                self.personal_api_anno.add_all_detail(database, self.personal_api, self.personal_api_anno.main_dic, self.need_run_result)

        if self.option('do_api_list') != '':
            part = copy.deepcopy(self.do_api_list)
            for d in ['probio_lca', 'probio_de_unclassified']:
                if d in self.do_api_list:
                    part.remove(d)
                    probio_main = {'probio':self.api_main_dic[d]}
                    probio_file = {'probio': self.api_files_dic[d]}
                    self.personal_api_anno.add_all_detail('probio', self.personal_api, probio_main, probio_file)
            if part !=[]:
                part_database = ';'.join(part)
                self.personal_api_anno.add_all_detail(part_database, self.personal_api, self.api_main_dic, self.api_files_dic)

        if 'lca' in self.option('nr_method').split(',') :  # deunclassied  or lca
            nr_path = os.path.join(self.remote_dir, 'nr_lca/gene_nr_anno.xls')
            return_id = self.nr_api.add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group,
                                                   group_detail=self.group_detail, name = "NR_Origin_LCA", task_id=self.option('task_id'))
            self.nr_api.add_anno_nr_detail(return_id, self.output_dir + "/nr_lca")
        elif 'lca' in self.nr_only_api:
            i = self.nr_only_api.index('lca')
            self.nr_api.add_anno_nr_detail(self.nr_main_id[i] ,self.option('nr_lca_dir').prop['path'])

        if 'de_unclassified' in self.option('nr_method').split(',') :
            nr_path = os.path.join(self.remote_dir, 'nr_deunclassied/gene_nr_anno.xls')
            return_id = self.nr_api.add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group,
                                           group_detail=self.group_detail, name = "NR_Origin_Deunclassified", task_id=self.option('task_id'))
            self.nr_api.add_anno_nr_detail(return_id, self.output_dir + "/nr_deunclassied")
        elif 'de_unclassied' in self.nr_only_api:
            i = self.nr_only_api.index('de_unclassied')
            self.nr_api.add_anno_nr_detail(self.nr_main_id[i] ,self.option('nr_de_dir').prop['path'])

        self.logger.info('开始导总览表')
        if os.path.exists(self.output_dir + '/gene_overview_anno.xls'):
            find_pre_sum = self.metagenomic.common_find_one('anno_overview',{'task_id':self.option('task_id')})
            self.anno_overview_id = self.overview_api.add_anno_overview(self.geneset_id,task_id=self.option('task_id'))
            anno_path = os.path.join(self.remote_dir, 'gene_overview_anno.xls')
            self.overview_api.add_anno_overview_detail(self.anno_overview_id, self.output_dir + '/gene_overview_anno.xls')
            #self.metagenomic.common_update_one('anno_overview', anno_overview_id, {'new_sum':anno_path})
            self.logger.info('总览表导完')
            if find_pre_sum:
                change_task_id = find_pre_sum['task_id'] + '_1'
                self.metagenomic.common_update_one('anno_overview',find_pre_sum['_id'],{'task_id':change_task_id})
                self.logger.info("将原主表更新完成")

                delete_task = self.metagenomic.common_find_one('anno_overview',{'task_id':change_task_id})
                if delete_task:
                    task_id = delete_task['task_id']
                    main_id = delete_task['_id']
                    self.logger.info("开始删除原详情表！")
                    self.logger.info("开始删除原详情表overview_detail表数据！")
                    self.metagenomic.rm_main_detail(task_id, 'anno_overview_detail', main_id, 'anno_overview_id')
                    self.logger.info("删除原详情表overview_detail表数据完成！")
                    self.logger.info("开始删除原主表")
                    self.metagenomic.rm_main(task_id, 'anno_overview', main_id)
                    self.logger.info("删除原来主表完成！")
        else:
            self.logger.info('没有路径：' + self.output_dir + '/gene_overview_anno.xls')
        # else:
        #     self.overview_api.add_anno_overview_detail(anno_overview_id, self.output_dir + '/gene_overview_anno.xls' )


        tmp_list = []
        tmp_list.extend(self.need_run_list)
        tmp_list.extend(self.do_api_list)
        tmp_list.extend(self.option('nr_method').split(','))
        for i in tmp_list:
            if i =='': continue
            if i == 'mvirdb': i = 'mvir'
            if i == 'p450': i = 'cyps'
            if i == 'probio_lca':
                result = self.metagenomic.common_find_one('anno_personal_run', {'task_id':self.option('task_id'), 'name':'Probio_Origin_LCA'})
            elif i == 'probio_de_unclassified':
                result = self.metagenomic.common_find_one('anno_personal_run', {'task_id':self.option('task_id'), 'name': 'Probio_Origin_Deunclassified'})
            else:
                result = self.metagenomic.common_find_one('anno_personal_run', {'task_id':self.option('task_id'), 'anno_type': i})
            if result:
                self.metagenomic.common_update_one('anno_personal_run', result['_id'], {'status':'end','desc':'任务结束'})
                if i in self.database_version:## 根据任务判断是新老项目，根据有没有用到nr注释结果来更新数据库版本
                    result2 = self.metagenomic.common_find_one("sg_task", {'task_id':self.option('task_id')})
                    if result2:
                        if not "database" in result2:
                            if i not in ['lca', 'de_unclassified', 'go', 'probio', 'sec', 't3ss']:
                                self.metagenomic.common_update_one('anno_personal_run', result['_id'], { "version":self.database_version[i]})
                        else:
                            self.metagenomic.common_update_one('anno_personal_run', result['_id'], { "version":self.database_version[i]})
                    else:
                        if i not in ['lca', 'de_unclassified', 'go', 'probio', 'sec', 't3ss']:
                            self.metagenomic.common_update_one('anno_personal_run', result['_id'], { "version":self.database_version[i]})
        self.end()


    def end(self):

        # if os.path.exists(self.output_dir + '/Sec'):
        #     for i in os.listdir(self.output_dir + '/Sec'):
        #         if 'SignalP.txt' in i:      #不能删，否则总览表报错
        #             os.remove(self.output_dir + '/Sec/'+i)

        # if os.path.exists(self.output_dir + '/Ttss'):
        #     for i in os.listdir(self.output_dir + '/Ttss'):
        #         if 'predict.txt' in i:            #不能删，否则总览表报错
        #             os.remove(self.output_dir + '/Ttss/'+i)

        #if self.option('summary_path')!= '':
        #    self.upload_to_s3(self.output_dir + '/gene_overview_anno.xls', self.option('summary_path')) ##
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            ["Tcdb", "", "TCDB功能注释结果目录", 0, "120058"],
            ["Tcdb/class_abund_out.xls", "xls", "各样品Class丰度表", 0, "120060"],
            ["Tcdb/subclass_abund_out.xls", "xls", "各样品SubClass丰度表",0,"120297"],
            ["Tcdb/family_abund_out.xls", "xls", "各样品family丰度表",0,"120298"],
            ["Tcdb/tcdb_abund_out.xls", "xls", "各样品tcdb丰度表",0,"120299"],
            ["Tcdb/Class_gene_stat.xls", "xls", "各样品Class基因列表", 0, "120060"],
            ["Tcdb/TCDB/SubClass_gene_stat.xls", "xls", "各样品SubClass基因列表",0,"120300"],
            ["Tcdb/Family_gene_stat.xls", "xls", "各样品family基因列表",0,"120301"],
            ["Tcdb/TCDB_gene_stat.xls", "xls", "各样品tcdb基因列表",0,"120302"],
            ["Tcdb/gene_tcdb_anno.xls", "xls", "每条基因的TCDB功能注释表",0,"120303"],
            ["Mvirdb", "", "MvirDB功能注释结果目录", 0, "120058"],
            ["Mvirdb/Factor_abund_out.xls", "xls", "各样品Factor丰度表", 0, "120060"],
            ["Mvirdb/Type_abund_out.xls", "xls", "各样品Type丰度表",0,"120304"],
            ["Mvirdb/Factor_gene_stat.xls", "xls", "各样品Factor基因列表",0,"120305"],
            ["Mvirdb/Type_gene_stat.xls", "xls", "各样品Type基因列表",0,"120306"],
            ["Mvirdb/gene_mvirdb_anno.xls", "xls", "每条基因的MvirDB功能注释表",0,"120307"],
            ["Probio", "", "益生菌功能注释结果目录",0,"120308"],
            ["Probio/gene_probio_anno.xls", "xls", "每条基因的益生菌注释表",0,"120309"],
            ["Go", "", "go功能注释结果目录",0,"120310"],
            ["Go/all.go1.function.xls", "xls", "level1水平各样品丰度表",0,"120311"],
            ["Go/all.go12.function.xls", "xls", "level2水平各样品丰度表",0,"120312"],
            ["Go/all.go123.function.xls", "xls", "level3水平各样品丰度表",0,"120313"],
            ["Go/all.go11234.function.xls", "xls", "level4水平各样品丰度表",0,"120314"],
            ["Go/all.go.annotation.xls", "xls", "各样品go注释gene对应注释表",0,"120315"],
            ["P450", "", "P450蛋白功能注释结果目录",0,"120316"],
            ["P450/gene_cyps_anno.xls", "xls", "每条蛋白的CYPS功能注释表",0,"120317"],
            ["P450/cyps_homo_profile.xls", "xls", "各样品CYPS Homologous_family丰度表",0,"120318"],
            ["P450/cyps_super_profile.xls", "xls", "各样品CYPS Superfamily丰度表",0,"120319"],
            ["P450/cyps_sid_profile.xls", "xls", "样品注释的最低层级表",0,"120320"],
            ["P450/cyps_anno_stat.xls", "xls", "pfam注释信息统计表",0,"120321"],
            ["Nr", "", "NR功能注释结果目录", 0, "120035"],
            ["Nr/gene_nr_anno.xls", "xls", "每条基因的物种注释表", 0, "120036"],
            ["Nr/nr_align_table.xls", "", "物种序列比对结果", 0, "120037"],
            ["Nr/tax_d.xls", "xls", "域注释丰度表", 0, "120038"],
            ["Nr/tax_k.xls", "xls", "界注释丰度表", 0, "120039"],
            ["Nr/tax_p.xls", "xls", "门注释丰度表", 0, "120040"],
            ["Nr/tax_c.xls", "xls", "纲注释丰度表", 0, "120041"],
            ["Nr/tax_o.xls", "xls", "目丰注释度表", 0, "120042"],
            ["Nr/tax_f.xls", "xls", "科注释丰度表", 0, "120043"],
            ["Nr/tax_g.xls", "xls", "属注释丰度表", 0, "120044"],
            ["Nr/tax_s.xls", "xls", "种注释丰度表", 0, "120045"],
            ["Pfam", "", "Pfam结构域注释结果目录",0,"120322"],
            ["Pfam/gene_pfam_anno.xls", "xls", "每条基因的Pfam功能注释表",0,"120323"],
            ["Pfam/pfam_type_profile.xls", "xls", "各样品的Pfam Type丰度表",0,"120324"],
            ["Pfam/pfam_acc_profile.xls", "xls", "各样品的Pfam Pfam丰度表/最低层级表",0,"120325"],
            ["Pfam/pfam_clan_profile.xls", "xls", "各样品的Pfam CLAN丰度表",0,"120326"],
            ["Pfam/gene_pfam_anno_stat.xls", "xls", "Pfam 注释基因统计",0,"120327"],
            ["Phi", "", "Phi结果目录",0,"120328"],
            ["Phi/gene_phi_anno.xls", "xls", "PHI功能注释结果信息表", 0, "120329"],
            ["Phi/phi_host_profile.xls", "xls", "各样品Host丰度表", 0, "120330"],
            ["Phi/phi_pathogen_profile.xls", "xls", "各样品Pathogen丰度表", 0, "120331"],
            ["Phi/phi_phenotype_profile.xls", "xls", "各样品Phenotype丰度表", 0, "120332"],
            ["Phi/phi_protein_profile.xls", "xls", "各样品Protein丰度表", 0, "120333"],
            ["Qs", "", "QS功能注释结果目录", 0, "120058"],
            ["Qs/qs_class_profile.xls", "xls", "各样品QS的class水平丰度表",0,"120334"],
            ["Qs/gene_qs_anno.xls", "xls", "各样品gene对应注释表",0,"120335"],
            ["Qs/qs_lowest_profile.xls", "xls", "各样品QS的最低水平丰度表",0,"120336"],
            ["Qs/anno_qs_graph.xls", "xls", "各样品QS画图的数据表",0,"120337"],
            ["Sec", "", "分泌蛋白预测结果目录", 0, "120338"],
            ["T3ss", "", "分泌蛋白预测结果目录",0,"120345"],
            ['gene_overview_anno.xls', 'xls', '个性化注释总览表',0,"120346"]
            #['gene_overview_anno_all.xls', 'xls', '注释总览表']
        ])
        result_dir.add_regexp_rules([
            [r"Sec/.*fisher\.txt", "txt", "分泌蛋白物种注释结果", 0, "120339"],
            [r"Sec/.*summary\.txt", "txt", "分泌蛋白个数统计", 0, "120340"],
            [r"Ttss/.*fisher\.txt", "txt", "Ttss注释结果",0,"120347"],
            [r"Ttss/.*summary\.txt", "txt", "Ttss个数统计",0,"120348"],
            [r"Ttss/.*predict\.txt", "txt", "Ttss基因列表",0,"120349"]
        ])
        if os.path.exists(self.output_dir + '/gene_overview_anno.xls'):  # update mongo table after results are uploaded to s3
            if self._sheet.id == 'i-sanger_183577_44543_128356':
                self.anno_overview_id = "5f3bd8a5a4e1af33c5b582e7"
            anno_path = os.path.join(self.remote_dir, 'gene_overview_anno.xls')
            self.metagenomic.common_update_one('anno_overview', self.anno_overview_id, {'new_sum': anno_path})

        super(PersonalAnnoWorkflow, self).end()

