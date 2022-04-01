# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""微生物基因组分析工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
import os,re
import json
import shutil
from biocluster.config import Config
import time
import datetime
import gevent
import functools
import shutil
from Bio import SeqIO
import types
from collections import defaultdict
import pandas as pd
from bson.objectid import ObjectId
from mainapp.models.mongo.bacgenome import Bacgenome
import glob
from mbio.packages.bacgenome.change_result_seq_name import ChangeResultSeqName
from mbio.packages.meta.delete_mongo import DeleteDemoMongo

def time_count(func):  # 统计导表时间
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run %s at %s' % (func_name, start_time))
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，该阶段导表已进行{}s".format(func_name, end - start))

    return wrapper

def tryforgood(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

class BacgenomeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        微生物基因组workflow option参数设置
        """
        self._sheet = wsheet_object
        super(BacgenomeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "raw_dir", "type": "infile", "format": "bacgenome.raw_dir"},  ###rawdata的文件目录
            {"name": "asse_dir", "type": "infile", "format": "bacgenome.asse_dir"},  ###组装的文件目录
            {"name": "analysis", "type": "string", "default": "uncomplete"}, ###流程分析模式complete，uncomplete
            {"name": "gff_dir", "type":"infile", "format":"bacgenome.gff_dir"},  ##   #bacgenome.gff_dir
            {"name": "fna_dir","type": "infile", "format": "bacgenome.fna_dir"},  # fna文件夹。注释流程
            {"name": "data_type", "type": "string"},
            {"name": "trans_code", "type" : "string","default":"11"},
            #{"name": "ncbi_id", "type": "string", "default":""},  # 注释流程，完成图选择NCBI的assemble id进行分析
            {"name": "nr_evalue","type":"string","default":"1e-5"},
            {"name": "swissprot_evalue","type":"string","default":"1e-5"},
            {"name": "cog_evalue","type":"string","default":"1e-5"},
            {"name": "kegg_evalue","type":"string","default":"1e-5"},
            {"name": "go_evalule","type": "string", "default":"1e-5"},
            {"name": "software_list", "type":"string", "default":"glimmer"} , #基因预测使用的软件，逗号分割
            {"name": "assemble_id","type":"string","default":""},  # ,分割， 注释流程的完成图的从NCBI中选择 的流程使用 #test ： GCA_000007165.1_ASM716v1
            {"name": "test", "type": "string","default":"T"},
            {"name": "data_from", "type": "string"},  # 注释流程 数据来源：根据ncib的编号或者上传文件
            {"name": "p_trans_code","type": "string","default":"11"}, #
            {"name": "p_software_list", "type": "string","default":"genemark"},  #质粒的基因预测使用的软件，逗号分割
            {'name': 'qc_tool', 'type': 'string', 'default': 'fastp'}, #质控软件选择
            {'name': 'assemble_tool', 'type': 'string'},  # 组装软件选择
            {'name': 'kmer', 'type': 'string', 'default': '21-47'},  # 质控软件选择

        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.samples ={}
        self.modules = []
        self.trees = []
        self.hgene_tree = self.add_tool('graph.phy_tree')
        self.tree = self.add_tool('graph.phy_tree')
        self.list =[self.hgene_tree,self.tree]
        self.ssu_align = self.add_tool("bacgenome.ssu_align")
        self.logger.info(self._sheet.output)
        self.step.add_steps('hgene_tree', 'tree', 'ssu_align')
        self.remote_dir = self._sheet.output
        self.samples_seq_map = {}
        self.samples_gene_map = {}

        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'bacgenome')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        """
        检查参数
        """

        if not self.option("analysis"):
            raise OptionError("请提供流程分析模式！", code="11400101")
        if self.option("gff_dir").is_set or self.option("assemble_id") !="" :
            self.is_anno_pip = True
        else:
            self.is_anno_pip = False
        if not self.is_anno_pip:
            if not self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                raise OptionError("必须输入原始序列文件夹或组装序列文件夹其中一个！", code="11400102")

        self.option('software_list', self.option('software_list').replace('genemarks','genemark'))
        self.option('p_software_list',self.option('p_software_list').replace('genemarks','genemark'))
        self.option('trans_code',self.option('trans_code').replace("genetic code ",""))
        self.option('p_trans_code',self.option('p_trans_code').replace("genetic code ",""))


    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def set_run(self, opts, module, event, step, start=True):
        module.set_options(opts)
        module.on('start', self.set_step, {'start': step})
        module.on('end', self.set_step, {'end': step})
        module.on('end', self.set_output, event)
        if start:
            module.run()

    def run_bacgenome(self):
        self.get_list()
        if len(self.samples.keys()) > 3:
            for sample in self.samples.keys():
                asse_path = self.work_dir + '/' + sample + '/assemble'
                raw_path = self.work_dir + '/' + sample + '/data'
                self.bacgenome = self.add_module('bacgenome.bacgenome')
                self.step.add_steps('bacgenome{}'.format(' analysis, sample: ' + sample))
                opts = ''
                if self.option("analysis") in ["uncomplete"]:
                    if os.path.exists(raw_path) and not os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'analysis': self.option('analysis'),
                        }
                    elif not os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                    elif os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                elif self.option("analysis") in ["complete"]:
                    if not os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                    elif os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                    elif os.path.exists(raw_path) and not os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'analysis': self.option('analysis'),
                        }
                opts['trans_code'] = self.option('trans_code')
                opts['nr_evalue'] = self.option('nr_evalue')
                opts['go_evalue'] = self.option('go_evalule')
                opts['cog_evalue'] = self.option('cog_evalue')
                opts['kegg_evalue'] = self.option('kegg_evalue')
                opts['swissprot_evalue'] = self.option('swissprot_evalue')
                opts['software_list'] = self.option('software_list')
                opts['p_software_list'] = self.option('p_software_list')
                opts['p_trans_code'] = self.option('p_trans_code')
                opts['qc_tool'] = self.option('qc_tool')
                self.bacgenome.set_options(opts)
                step = getattr(self.step, 'bacgenome{}'.format(' analysis, sample: ' + sample))
                step.start()
                self.step.update()
                self.bacgenome.on('end', self.set_output, 'bacgenome{}'.format(' analysis, sample: ' + sample))
                self.bacgenome.on('end', self.finish_update, 'bacgenome{}'.format(' analysis, sample: ' + sample))
                self.modules.append(self.bacgenome)
            self.logger.info(self.modules)
            if self.option("analysis") in ["uncomplete"]:
                self.on_rely(self.modules, self.run_hgene_tree)
            elif self.option("analysis") in ["complete"]:
                self.on_rely(self.modules, self.run_tree)
            for module in self.modules:
                module.run()
                gevent.sleep(0)
        elif len(self.samples.keys()) <= 3:
            for sample in self.samples.keys():
                asse_path = self.work_dir + '/' + sample + '/assemble'
                raw_path = self.work_dir + '/' + sample + '/data'
                self.bacgenome = self.add_module('bacgenome.bacgenome')
                self.step.add_steps('bacgenome{}'.format(' analysis, sample: ' + sample))
                opts = ''
                if self.option("analysis") in ["uncomplete"]:
                    if os.path.exists(raw_path) and not os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'analysis': self.option('analysis'),
                        }
                    elif not os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                    elif os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                elif self.option("analysis") in ["complete"]:
                    if not os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                    elif os.path.exists(raw_path) and os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'asse_dir': asse_path,
                            'analysis': self.option('analysis'),
                        }
                    elif os.path.exists(raw_path) and not os.path.exists(asse_path):
                        opts = {
                            'sample_name': sample,
                            'raw_dir': raw_path,
                            'analysis': self.option('analysis'),
                        }
                opts['trans_code'] = self.option('trans_code') #
                opts['nr_evalue'] = self.option('nr_evalue')
                opts['go_evalue'] = self.option('go_evalule')
                opts['cog_evalue'] = self.option('cog_evalue')
                opts['kegg_evalue'] = self.option('kegg_evalue')
                opts['swissprot_evalue'] = self.option('swissprot_evalue')
                opts['software_list'] = self.option('software_list')
                opts['p_software_list'] = self.option('p_software_list')
                opts['p_trans_code'] = self.option('p_trans_code')
                self.bacgenome.set_options(opts)
                step = getattr(self.step, 'bacgenome{}'.format(' analysis, sample: ' + sample))
                step.start()
                self.step.update()
                self.bacgenome.on('end', self.finish_update, 'bacgenome{}'.format(' analysis, sample: ' + sample))
                self.modules.append(self.bacgenome)
            self.logger.info(self.modules)
            if 1 < len(self.samples.keys()) <= 3:
                self.on_rely(self.modules, self.set_output)
            elif len(self.samples.keys()) == 1:
                self.modules[0].on('end', self.set_output)
            for module in self.modules:
                module.run()

    def run_bacgenome_2(self):
        if self.option('assemble_id'):
            self.pre_anno_pipline_from_ncbi()
        else:
            self.pre_anno_pipeline()
        self.write_map_file()  #生成 序列对应的新序列名称 map file
        for sample in self.samples.keys():
            self.bacgenome = self.add_module('bacgenome.bacgenome_anno')
            opts = {
                "seq_dir" :  self.work_dir +'/'+sample + '/seq_dir',
                "gff_file" : self.work_dir + '/'+sample + '/' + sample+'.gff',
                "all_fasta" : self.work_dir + '/'+sample + '/all.fasta',
                "sample_name" : sample,
                "analysis" : self.option("analysis"),
            }
            if self.option("analysis") in ['complete']:
                opts['txt_info'] = str(self.com_chr_pls_num[sample])
            opts['nr_evalue'] = self.option('nr_evalue')
            opts['go_evalue'] = self.option('go_evalule')
            opts['cog_evalue'] = self.option('cog_evalue')
            opts['kegg_evalue'] = self.option('kegg_evalue')
            opts['swissprot_evalue'] = self.option('swissprot_evalue')
            self.bacgenome.set_options(opts)
            if len(self.samples.keys()) >3:
                self.bacgenome.on('end', self.set_output, 'bacgenome{}'.format(' analysis, sample: ' + sample))
            self.modules.append(self.bacgenome)

    def run_tree(self):
        self.ssu_align.on("end",self.run_16s_tree)
        self.run_hgene_tree()
        self.run_ssu_align()

    def run_hgene_tree(self):
        if os.path.exists(self.work_dir + '/all.cor_gene.fasta'):
            os.remove(self.work_dir + '/all.cor_gene.fasta')
        files = glob.glob(self.work_dir+'/Bacgenome*/output/*/tree/hgene/*.corgene.fa')
        for file in files:
            os.system('cat %s >> %s' % (file, self.work_dir + '/all.cor_gene.fasta'))
        opts = {
            "fasta": self.work_dir + '/all.cor_gene.fasta',
            "sequence_type": 'amino_acid',
        }
        self.set_run(opts, self.hgene_tree, 'hgene_tree', self.step.hgene_tree)

    def run_16s_tree(self):
        opts = {
            "fasta": self.ssu_align.option('out'),
        }
        self.set_run(opts, self.tree, 'tree', self.step.tree)

    def run_ssu_align(self):
        if os.path.exists(self.work_dir + '/all.16s.fasta'):
            os.remove(self.work_dir + '/all.16s.fasta')
        files = glob.glob(self.work_dir + '/Bacgenome*/output/*/tree/16s/*.16s.fa')
        if len(files) >=4:
            for file in files:
                os.system('cat %s >> %s' % (file, self.work_dir + '/all.16s.fasta'))
            opts = {
                "16s_fa": self.work_dir + '/all.16s.fasta'
            }
            self.set_run(opts, self.ssu_align, "ssu_align", self.step.ssu_align)

    def run(self):
        """
        运行:genome_workflow
        :return:
        """
        if self.sheet.id in ["tsg_248925"]:
            self.IMPORT_REPORT_DATA = False
            self.get_list()
            gevent.spawn_later(5, self.end)
            super(BacgenomeWorkflow, self).run()
            return
        else:
            self.logger.info("a1")
            if not self.is_anno_pip:
                self.logger.info("a2")
                self.run_split_dir()
            self.logger.info("a3")
            time.sleep(2)
            task_info = self.api.api('task_info.bacg_task_info')
            self.logger.info("a4")
            task_info.add_task_info()
            self.logger.info("a5")
            if not self.is_anno_pip:
                self.logger.info("aaa"+str(len(self.samples.keys())))
                if len(self.samples.keys()) > 3:
                    if self.option("analysis") in ["uncomplete"]:
                        self.hgene_tree.on('end', self.end)
                        self.run_bacgenome()
                    elif self.option("analysis") in ["complete"]:
                        self.on_rely(self.list, self.end)
                        self.run_bacgenome()
                elif 1 < len(self.samples.keys()) <= 3:
                    self.run_bacgenome()
                elif len(self.samples.keys()) == 1:
                    self.logger.info("bbb" + str(len(self.samples.keys())))
                    self.run_bacgenome()
            else:
                self.run_bacgenome_2()
                sample_num = len(self.modules)
                if sample_num > 3 :
                    if self.option("analysis") in ["uncomplete"]:
                        self.on_rely(self.modules, self.run_hgene_tree)
                        self.hgene_tree.on('end', self.end)
                    elif self.option("analysis") in ["complete"]:
                        self.on_rely(self.modules, self.run_tree)
                        self.on_rely(self.list, self.end)
                elif 1<sample_num <= 3:
                    self.on_rely(self.modules, self.set_output)
                elif sample_num == 1:
                    self.modules[0].on('end',self.set_output)
                else:
                    raise Exception('sample num is 0')

                for module in self.modules:
                    module.run()
            super(BacgenomeWorkflow, self).run()

    def set_output(self, event):
        """
        将各个模块的结果输出至output
        """
        self.logger.info('execute set_output')
        self.logger.info('sample num {}'.format(len(self.samples.keys())))
        if len(self.samples.keys()) <= 3:
            self.logger.info('开始移动文件到workflow的output')
            for module in self.modules:
                self.move_dir(module.output_dir,self.output_dir)
            self.end()
        elif len(self.samples.keys()) > 3:
            obj = event["bind_object"]
            self.logger.info(obj)
            self.logger.info(event["data"])
            if re.match(r"bacgenome", event["data"]):
                self.move_dir(obj.output_dir, self.output_dir)
            if self.option("analysis") in ["uncomplete"]:
                if event['data'] == 'hgene_tree':
                    if os.path.exists(self.output_dir + '/all.cor_gene.phylo_tree.nwk'):
                        os.remove(self.output_dir + '/all.cor_gene.phylo_tree.nwk')
                    os.link(self.hgene_tree.work_dir + '/phylo_tree.nwk',
                            self.output_dir + '/all.cor_gene.phylo_tree.nwk')
            elif self.option("analysis") in ["complete"]:
                if event['data'] == 'hgene_tree':
                    if os.path.exists(self.output_dir + '/all.cor_gene.phylo_tree.nwk'):
                        os.remove(self.output_dir + '/all.cor_gene.phylo_tree.nwk')
                    os.link(self.hgene_tree.work_dir + '/phylo_tree.nwk',
                            self.output_dir + '/all.cor_gene.phylo_tree.nwk')
                if event['data'] == 'tree':
                    if os.path.exists(self.output_dir + '/all.16s.phylo_tree.nwk'):
                        os.remove(self.output_dir + '/all.16s.phylo_tree.nwk')
                    os.link(self.tree.work_dir + '/phylo_tree.nwk', self.output_dir + '/all.16s.phylo_tree.nwk')

    def wait_file(self, path, wait_times=1):
        '''
        增加等待文件结果方法
        :param path: 结果文件路径
        :param wait_times: 等待次数
        :return: 文件路径
        :time: 20180425
        '''
        self.logger.info(">>>in wait_file")
        while wait_times < 11:
            if not os.path.exists(path):
                time.sleep(10)
                wait_times += 1
                self.wait_file(path, wait_times=wait_times)
            else:
                self.logger.info(">>> esists path:%s" % path)
            return path
        self.logger.info("超过文件等待次数，需检查文件%s" % path)
        return

    def end(self):
        if self.sheet.id in ["tsg_248925"]:
            self.send_files()
            super(BacgenomeWorkflow, self).end()
            return
        else:
            self.run_api()
            self.send_files()
            super(BacgenomeWorkflow, self).end()

    def run_api(self, test=False):
        self.stop_timeout_check()
        self.change_back_result_files_name()  ##注释流程 返回序列名称
        task_id =self._sheet.id
        project_sn =self._sheet.project_sn
        self.gbk = self.api.api("bacgenome.gbk")
        self.srna = self.api.api("bacgenome.srna")
        self.plasmid = self.api.api("bacgenome.plasmid")
        self.house_keep = self.api.api("bacgenome.house_keep")
        self.circos = self.api.api("bacgenome.circos")
        self.crispr = self.api.api("bacgenome.crispr")
        self.island = self.api.api("bacgenome.island")
        self.prephage = self.api.api("bacgenome.prephage")
        self.datastat = self.api.api("bacgenome.genome_qc")
        self.genome_size = self.api.api("bacgenome.genome_size")
        self.assess_gc = self.api.api("bacgenome.assess_gc")
        self.assemble = self.api.api("bacgenome.assemble")
        self.assess_kmer = self.api.api("bacgenome.assess_kmer")
        self.anno_antismash = self.api.api("bacgenome.anno_antismash")
        self.anno_card = self.api.api("bacgenome.anno_card")
        self.anno_cazy = self.api.api("bacgenome.anno_cazy")
        self.anno_cog = self.api.api("bacgenome.anno_cog")
        self.anno_go = self.api.api("bacgenome.anno_go")
        self.anno_kegg = self.api.api("bacgenome.anno_kegg")
        self.anno_nr = self.api.api("bacgenome.anno_nr")
        self.anno_pfam = self.api.api("bacgenome.anno_pfam")
        self.anno_phi = self.api.api("bacgenome.anno_phi")
        self.anno_ref = self.api.api("bacgenome.anno_ref")
        self.anno_summary = self.api.api("bacgenome.anno_summary")
        self.anno_swissprot = self.api.api("bacgenome.anno_swissprot")
        self.anno_tcdb = self.api.api("bacgenome.anno_tcdb")
        self.anno_tmhmm = self.api.api("bacgenome.anno_tmhmm")
        self.anno_vfdb = self.api.api("bacgenome.anno_vfdb")
        self.gene_graph = self.api.api("bacgenome.gene_graph")
        self.gene_predict = self.api.api("bacgenome.gene_predict")
        self.promote = self.api.api("bacgenome.promote")
        self.repeat_predict = self.api.api("bacgenome.repeat_predict")
        self.rrna_predict = self.api.api("bacgenome.rrna_predict")
        self.secretory = self.api.api("bacgenome.secretory")
        self.summary_map = self.api.api("bacgenome.summary_map")
        self.trna_predict = self.api.api("bacgenome.trna_predict")
        self.two_component_api= self.api.api("bacgenome.common_api")
        self.methylation_api = self.api.api("bacgenome.methylation")
        self.software = self.api.api("bacgenome.software")
        self.api_integron = self.api.api('bacgenome.anno_integron')
        self.api_is = self.api.api('bacgenome.anno_is')
        self.up_repeatmasker = self.api.api('bacgenome.repeatmasker')
        self.api_resfinder = self.api.api('bacgenome.resfinder')

        self.bacgenome_model = Bacgenome()
        self.bacgenome_model._config = Config()
        self.bacgenome_model.update_mongo('sg_task',{"task_id":task_id}, {"version":"3.1"})
        database_type = json.dumps({
            "nr": "v20200604",
            "cazy": "V8",
            "card": "v3.0.9",
            "kegg": "v94.2",
            "vfdb": "v20200703",
            'pfam': "v33.1",
            "swissprot": "v20200617"
        },sort_keys=True, separators=(',', ':'))
        self.bacgenome_model.update_mongo('sg_task',{"task_id":task_id}, {"database":database_type})

        ### 将软件与任务关联
        soft_params = {"version": "update_database_202011"}
        new_main_id = self.software.add_software(task_id=task_id,params=soft_params)
        self.software.add_software_detail(new_main_id)
        self.logger.info("software finish")

        if self.option('raw_dir').is_set:
            seq_type = self.get_assemble_type(self.option('raw_dir').prop['path'] + '/list.txt')

        if self.option('qc_tool') == "fastp":
            qc_software = "Fastp"
        elif self.option('qc_tool') == "old_mode":
            qc_software = "Trimmomatic, SeqPrep, Sickle, FastqTotalHighQualityBase.jar"
        else:
            qc_software = ""
        if not self.is_anno_pip:
            ############ 质控
            if self.option("analysis") in ["uncomplete"]:
                if self.option("raw_dir").is_set:
                    datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                             qc_software, "30", "20")
                elif self.option("asse_dir").is_set:
                    datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', "", "assemble", "", "",
                                                             "")
            elif self.option("analysis") in ["complete"]:
                if self.option("raw_dir").is_set:
                    if re.search(r'PE', seq_type):
                        datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                                 qc_software, "30", "20")
                    else:
                        datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', seq_type, "rawdata",
                                                                 "smrtanalysis", "", "")
                    if re.search(r'Pacbio', seq_type) or re.search(r'pacbio', seq_type):
                        methylation_id = self.methylation_api.add_methy()

                elif self.option("asse_dir").is_set:
                    datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', "", "assemble", "", "","")

            ############ 组装
            if self.option("analysis") in ["uncomplete"]:
                if self.option("raw_dir").is_set:
                    assemble_id = self.assemble.add_assemble(self.option('analysis'), seq_type, 'rawdata', '组装评估',
                                                             'kmer: {}'.format(self.option('kmer')), 'SOAPdenovo v2.04,GapCloser v1.12')
                elif self.option("asse_dir").is_set:
                    assemble_id = self.assemble.add_assemble(self.option('analysis'), '', 'assemble', '组装评估', '', '')
            elif self.option("analysis") in ["complete"]:
                if self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                    assemble_id = self.assemble.add_assemble(self.option('analysis'), seq_type, 'rawdata', '组装评估', '', '')
                elif not self.option("raw_dir").is_set and self.option("asse_dir").is_set:
                    assemble_id = self.assemble.add_assemble(self.option('analysis'), '', 'assemble', '组装评估', '', '')
                if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                    assemble_id = self.assemble.add_assemble(self.option('analysis'), seq_type, 'rawdata', '组装评估', '', 'Unicycler')
            ########### 预测
            if self.option("analysis") in ["uncomplete"]:
                gene_id = self.gene_predict.add_gene_predict(self.remote_dir , '/assembly_predict/predict/CDS_predict/',
                                                             '_CDS.', "Prodigal")
            elif self.option("analysis") in ["complete"]:
                gene_id = self.gene_predict.add_gene_predict(self.remote_dir , '/assembly_predict/predict/CDS_predict/',
                                                             '_whole_genome_CDS.', "Prodigal and GeneMarkS")
            trna_id = self.trna_predict.add_trna_predict(self.remote_dir, '/assembly_predict/predict/tRNA/', '_tRNA.',
                                                         params='tRNAscan-SE')
            rrna_id = self.rrna_predict.add_rrna_predict(self.remote_dir, '/assembly_predict/predict/rRNA/', '_rRNA.',
                                                         params='barrnap V0.9')
            repeat_id = self.repeat_predict.add_repeat_predict(params='TRF')
            if (self.option("asse_dir").is_set) and (not self.option("raw_dir").is_set):## add by qingchen.zhang@20201214原因是只上传组装序列的时候不跑gc和基因组大小
                pass
            else:
                gc_id = self.assess_gc.add_assess_gc("GC_Depth",link_path =self.remote_dir)
                size_id = self.genome_size.add_assess_size("genome size")
        ##注释流程
        else:
            datastat_id = self.datastat.add_datastat(self.option('analysis'), 'datastat主表', "", "assemble", "", "","")
            gene_id = self.gene_predict.add_gene_predict(self.remote_dir , '/gene/', '.cds.', 'none', params='{"soft": ""}')  #注释流程，加类似基因预测的结果。后续交互分析需要
            assemble_id = self.assemble.add_assemble_anno_pipline(analysis_type=self.option('analysis'))

        ## 注释
        nr_id = self.anno_nr.add_anno_nr(params={"NR":"Diamond","evalue":self.option('nr_evalue')})
        cog_id = self.anno_cog.add_anno_cog(params={"COG":"Diamond","evalue":self.option('cog_evalue')})
        gbk_id = self.gbk.add_gbk(project_sn=project_sn, task_id=task_id)
        go_id = self.anno_go.add_anno_go(params={"GO":"blast2go","evalue":self.option('go_evalule')})
        antismash_id =self.anno_antismash.add_anno_antismash(params={"antismash":"antismash"})
        cazy_id = self.anno_cazy.add_anno_cazy(params={"CAZy":"hmmscan"})
        kegg_id = self.anno_kegg.add_anno_kegg(self.remote_dir, "/annotation/KEGG/", "_kegg_pathway_img", params={"KEGG":"Diamond","evalue":self.option("kegg_evalue")})
        pfam_id = self.anno_pfam.add_anno_pfam(params={"Pfam":"hmmer3"})
        swissport_id = self.anno_swissprot.add_anno_swissprot(params={"Swiss-prot":"blast","evalue":self.option("swissprot_evalue")})
        summary_id = self.anno_summary.add_anno_summary(params="{NR,KEGG,COG,GO,Swissprot,Pfam}")
        island_id = self.island.add_island("isand")
        prephage_id = self.prephage.add_prephage("prephage","Phigaro")
        crispr_id = self.crispr.add_crispr("crispr")
        two_component_id = self.two_component_api.add_main('anno_regulator', name='two_component')
        integron_id = self.api_integron.add_anno_integron(project_sn=project_sn, task_id=task_id,params={"submit_location":"integron","task_id":task_id,"task_type":2}, specimen_id=",".join(self.samples.keys()))
        is_id = self.api_is.add_anno_is(project_sn=project_sn, task_id=task_id,params={"submit_location":"is","task_id":task_id,"task_type":2}, specimen_id=",".join(self.samples.keys()))
        interpersed_id = self.up_repeatmasker.add_repeat_predict(params ={"submit_location":"interpersedrepeat","task_id":task_id,"task_type":2})
        res_id = self.api_resfinder.add_resfinder(params={"specimen_id": ",".join(self.samples.keys()),"identity": float(80),"coverage": float(60),"species_name": "all","task_id": task_id, "task_type": 2,"submit_location": "resfinder"})
        srna_id = self.srna.add_srna("rRNA_predict", "Infernal", "Rfam")
        plasmid_id = self.plasmid.add_plasmid("plasmid_predict", "PlasFlow、BLAST", "PLSDB")
        house_keep_id = self.house_keep.add_house("house_keep", "BLAST", "House-keeping Gene")

        ### 导详情表
        n = 1
        for sample in self.samples.keys():
            self.logger.info(">>>start wait_file,file route: %s" % (self.output_dir + '/' + sample))
            self.wait_file(self.output_dir + '/' + sample)
            self.logger.info(">>>wait_file end")
            ##  原流程和注释流程共用
            #gbk
            if self.option("analysis") in ["uncomplete"]:
                gbk_fi = os.path.join(self.remote_dir,sample + "/project_overview/Analysis_files/")
                self.gbk.add_gbk_detail(path = self.output_dir + '/' + sample + '/gbk/gbk',main_id =gbk_id,sample_name= sample,gbk_path= gbk_fi,task_id=task_id)
            elif self.option("analysis") in ["complete"]:
                gbk_fi = os.path.join(self.remote_dir,sample + "/project_overview/Analysis_files/")
                self.gbk.add_gbk_detail(path=self.output_dir + '/' + sample + '/gbk/seq_gbk', main_id=gbk_id,sample_name=sample, gbk_path=gbk_fi, task_id=task_id)
            #project.summary
            if self.option("analysis") in ["uncomplete"]:
                self.datastat.add_datastat_uncomplete_summary(datastat_id,
                                                              self.output_dir + '/' + sample + '/' + sample + '.project.summary')
            elif self.option("analysis") in ["complete"]:
                self.datastat.add_datastat_complete_summary(datastat_id,
                                                            self.output_dir + '/' + sample + '/' + sample + '.project.summary')
            ## tree
            if self.option("analysis") in ["uncomplete"]:
                for type in ['hgene']:
                    if os.path.exists(self.output_dir + '/' + sample + '/tree/' + type + '/' + sample + '.phylo_tree.nwk'):
                        self.datastat.add_datastat_tree_detail(datastat_id,
                                    self.output_dir + '/' + sample + '/tree/' + type + '/' + sample + '.phylo_tree.nwk', sample, type, 'single')

            elif self.option("analysis") in ["complete"]:
                for type in ['16s', 'hgene']:
                    if os.path.exists(self.output_dir + '/' + sample + '/tree/' + type + '/' + sample + '.phylo_tree.nwk'):
                        self.datastat.add_datastat_tree_detail(datastat_id,
                                        self.output_dir + '/' + sample + '/tree/' + type + '/' + sample + '.phylo_tree.nwk',sample, type, 'single')

            ##核心基因进化树比对的结果
            if os.path.exists(self.output_dir + '/' + sample + '/tree/hgene/'+ sample + '.house_blast.xls'):
                self.datastat.add_datastat_blast(self.output_dir + '/' + sample + '/tree/hgene/'+ sample + '.house_blast.xls',datastat_id,"hgene") #zouguanqing
            if os.path.exists(self.output_dir + '/' + sample + '/tree/16s/'+ sample + '.16s_blast.xls'):
                self.datastat.add_datastat_blast(self.output_dir + '/' + sample + '/tree/16s/'+ sample + '.16s_blast.xls',datastat_id,"16s") #gaohao

            if self.is_anno_pip:
                if self.option("analysis") in ["complete"]:
                    self.datastat.add_gene(datastat_id, sample, gene_pre=self.samples_gene_map, map=self.samples_seq_map[sample]["map"])  ##20190709
                else:
                    self.datastat.add_gene(datastat_id, sample, gene_pre=self.samples_gene_map)
                self.datastat.add_sample(datastat_id,sample)

                ## 基因结果
                self.gene_predict.add_gene_predict_detail(gene_id, sample, self.output_dir + '/' + sample + '/gene/' + sample + ".gene.gff", anno_pipeline=True)

                # self.gene_predict.add_gene_predict_specimen(gene_id, sample,
                #                                             self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS_statistics.xls")
                self.gene_predict.add_gene_predict_seq(gene_id, sample,
                                                       self.output_dir + '/' + sample + '/gene/' + sample + ".cds.ffn",
                                                       self.output_dir + '/' + sample + '/gene/' + sample + ".cds.faa")
                # self.gene_predict.add_gene_predict_bar(gene_id, sample,
                #                                        self.output_dir + '/' + sample + '/assembly_predict/predict/' + "length_distribute.txt")
                self.assemble.add_assemble_anno_pipline(main_id=assemble_id, specimen_id=sample,
                                                        len_file=self.output_dir + '/' + sample + '/gene/'+sample+".sequence_len.xls") #20191202
            if not self.is_anno_pip:
                ## 增加上传的样本层级上的gff文件和genome文件
                genome_path = self.remote_dir  + sample + '/' + sample + '.all.fna'
                gff_path = self.remote_dir + sample + '/' + sample + '.all.gff'
                if self.option("analysis") in ["uncomplete"]:
                    self.datastat.add_datastat_specimen(datastat_id,
                                                        self.output_dir + '/' + sample + '/' + 'specimen.data', genome_path=genome_path, gff_path=gff_path)
                    self.datastat.add_datastat_uncomplete_gene(datastat_id,
                                                               self.output_dir + '/' + sample + '/' + 'gene.data')
                elif self.option("analysis") in ["complete"]:
                    self.datastat.add_datastat_specimen(datastat_id,
                                                        self.output_dir + '/' + sample + '/' + 'specimen.data', genome_path=genome_path, gff_path=gff_path)
                    self.datastat.add_datastat_complete_gene(datastat_id,
                                                             self.output_dir + '/' + sample + '/' + 'gene.data')
                ####质粒注释
                if os.path.exists(self.output_dir + '/'+sample+'/plasmid_predict/all.stat.xls'):
                    seq_dict = self.get_seq_dict(self.output_dir + '/'+sample+'/plasmid_predict/all.stat.xls')
                self.plasmid.add_plasmid_stat(plasmid_id, self.output_dir + '/'+sample+'/plasmid_predict/'+ sample+'.stat.xls')
                if os.path.exists(self.output_dir + '/'+sample+'/plasmid_predict/'+ sample+'.plasmid_anno.xls'):
                    self.plasmid.add_plasmid_anno(plasmid_id, self.output_dir + '/'+sample+'/plasmid_predict/'+ sample+'.plasmid_anno.xls', sample)
                #####sRNA
                if os.path.exists(self.output_dir + '/' + sample + '/srna_predict/'+sample+".stat.xls"):
                    self.srna.add_srna_stat(srna_id, self.output_dir + '/' + sample + '/srna_predict/'+sample+".stat.xls", sample)
                    self.srna.add_srna_detail(srna_id, self.output_dir + '/' + sample + '/srna_predict/' + sample + ".sRNA.xls")
                ####house_keep
                if os.path.exists(self.output_dir + '/' + sample + '/tree/hgene/' + sample + ".cor_blast.xls"):
                    self.house_keep.add_house_detail(house_keep_id, self.output_dir + '/' + sample + '/tree/hgene/' + sample + ".cor_blast.xls", sample)
                ##甲基化
                if self.option("analysis") in ['complete'] and self.option('raw_dir').is_set:
                    motif_csv = self.output_dir + '/'+sample+'/methylation/'+ sample+'.motifs.csv'
                    motif_detail = self.output_dir + '/'+sample+'/methylation/'+ sample+'.motif_detail.xls'
                    if os.path.exists(motif_csv):
                        self.methylation_api.add_methy_stat(methylation_id, sample, motif_csv)
                        self.methylation_api.add_methy_detail(methylation_id,sample, motif_detail)
                        self.bacgenome_model.update_mongo('sg_task',{"task_id":task_id}, {"methylation":1})
                    else:
                        self.bacgenome_model.update_mongo('sg_task',{"task_id":task_id}, {"methylation":0})
                    motif_gff = self.output_dir + '/'+sample+'/methylation/'+ sample+'.motif.gff'
                    if os.path.exists(motif_gff):# 新的甲基化分析多了一个文件motifmaker ## 没有找到版本
                        self.bacgenome_model.update_mongo('methylation',{"task_id":task_id},
                                                          {"software":"MotifMaker_v0.1"})

                ### 质控
                if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                    my_seq = self.get_lib_type(self.option('raw_dir').prop['path'] + '/list.txt')
                    self.logger.info(my_seq)
                    if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls'):
                        self.datastat.add_qc_stat_uncomplete(datastat_id,
                                                             self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls',self.option('qc_tool'))
                    if sample in my_seq.keys():
                        self.logger.info(my_seq[sample].keys())
                        for type in my_seq[sample].keys():
                            self.logger.info(type)
                            for s in ['raw', 'clean']:
                                self.datastat.add_datastat_graphic(datastat_id,
                                                                   self.output_dir + '/' + sample + '/fastx/' + type + '_l.' + s + '_fastxstat',
                                                                   sample, s, "left", type)
                                self.datastat.add_datastat_graphic(datastat_id,
                                                                   self.output_dir + '/' + sample + '/fastx/' + type + '_r.' + s + '_fastxstat',
                                                                   sample, s, "right", type)

                elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                    my_seq = self.get_lib_type(self.option('raw_dir').prop['path'] + '/list.txt')
                    if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls'):
                        self.datastat.add_qc_stat_uncomplete(datastat_id,
                                                                  self.output_dir + '/' + sample + '/data_QC/' + sample + '_Illumina_statistics.xls',self.option('qc_tool'))
                        if sample in my_seq.keys():
                            for type in my_seq[sample].keys():
                                for s in ['raw', 'clean']:
                                    self.datastat.add_datastat_graphic(datastat_id,
                                                                       self.output_dir + '/' + sample + '/fastx/' + type + '_l.' + s + '_fastxstat',
                                                                       sample, s, "left", type)
                                    self.datastat.add_datastat_graphic(datastat_id,
                                                                       self.output_dir + '/' + sample + '/fastx/' + type + '_r.' + s + '_fastxstat',
                                                                       sample, s, "right", type)
                    if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '.PacBio_statistics.xls'):
                        self.datastat.add_qc_stat_complete(datastat_id,
                                                                self.output_dir + '/' + sample + '/data_QC/' + sample + '.PacBio_statistics.xls',
                                                                'pacbio',sample)
                        for s in ['raw', 'clean']:
                            if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '.' + s + '.len.xls'):
                                self.datastat.add_datastat_pacbio_graphic(datastat_id,self.output_dir + '/' + sample + '/data_QC/' +sample+ '.' + s + '.len.xls',sample, s, "len")
                            if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + 'pacbio.' + s + '.qv.xls'):
                                self.datastat.add_datastat_pacbio_graphic(datastat_id, self.output_dir + '/' + sample + '/data_QC/' + 'pacbio.' + s + '.qv.xls',sample, s, "qv")
                    elif os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '.Nanopore_statistics.xls'):
                        self.datastat.add_qc_stat_complete(datastat_id,
                                                                self.output_dir + '/' + sample + '/data_QC/' + sample + '.Nanopore_statistics.xls',
                                                                'nanopore',sample)
                        for s in ['raw', 'clean']:
                            if os.path.exists(self.output_dir + '/' + sample + '/data_QC/' + sample + '.' + s + '.len.xls'):
                                self.datastat.add_datastat_pacbio_graphic(datastat_id,self.output_dir + '/' + sample + '/data_QC/' +sample+ '.' + s + '.len.xls',sample, s, "len")
                if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "1k",self.remote_dir + sample + "/genomic_assessment/depth_gc_1000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "3k",self.remote_dir + sample + "/genomic_assessment/depth_gc_3000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "5k",self.remote_dir + sample + "/genomic_assessment/depth_gc_5000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "8k",self.remote_dir + sample + "/genomic_assessment/depth_gc_8000/")
                    self.assess_gc.add_assess_gc_detail(gc_id, sample, "10k",self.remote_dir + sample + "/genomic_assessment/depth_gc_10000/")
                    kmer_id = self.assess_kmer.add_assess_kmer(sample, "kmer")
                    self.assess_kmer.add_assess_kmer_detail(kmer_id,
                                                            self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + ".frequency.xls")
                    self.genome_size.add_assess_size_detail(size_id,
                                                            self.output_dir + '/' + sample + "/genomic_assessment/genome_size/" + sample + ".summary.xls",
                                                            sample)
                elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                    if re.search(r'PE', seq_type):
                        self.assess_gc.add_assess_gc_detail(gc_id, sample, "1k",self.remote_dir + sample + "/genomic_assessment/depth_gc_1000/")
                        self.assess_gc.add_assess_gc_detail(gc_id, sample, "3k",self.remote_dir + sample + "/genomic_assessment/depth_gc_3000/")
                        self.assess_gc.add_assess_gc_detail(gc_id, sample, "5k",self.remote_dir + sample + "/genomic_assessment/depth_gc_5000/")
                        self.assess_gc.add_assess_gc_detail(gc_id, sample, "8k",self.remote_dir + sample + "/genomic_assessment/depth_gc_8000/")
                        self.assess_gc.add_assess_gc_detail(gc_id, sample, "10k",self.remote_dir + sample + "/genomic_assessment/depth_gc_10000/")
                        kmer_id = self.assess_kmer.add_assess_kmer(sample, "kmer")
                        self.assess_kmer.add_assess_kmer_detail(kmer_id,
                                                                self.output_dir + '/' + sample + "/genomic_assessment/kmer_frequency/" + sample + ".frequency.xls")
                        self.genome_size.add_assess_size_detail(size_id,
                                                                self.output_dir + '/' + sample + "/genomic_assessment/genome_size/" + sample + ".summary.xls",
                                                                sample)
                ## 组装
                if self.option("analysis") in ["uncomplete"]:
                    if self.option('raw_dir').is_set:
                        self.assemble.add_assemble_stat_uncomplete(assemble_id,
                                                                   self.output_dir + '/' + sample + '/assembly_predict/assembly/assembly/' + sample + '_assembly_summary.xls',
                                                                   sample, file2 =self.output_dir + '/' + sample + '/assembly_predict/assembly/busco/' + sample + '_busco.xls')
                    else:
                        self.assemble.add_assemble_stat_uncomplete(assemble_id,
                                                               self.output_dir + '/' + sample + '/assembly_predict/assembly/assembly/' + sample + '_assembly_summary.xls',
                                                               sample)
                    dict12 = self.get_seq_dict2(self.output_dir + '/' + sample + '/plasmid_predict/all.stat.xls')
                    for type in ['scaffold', 'contig']:
                        if type == "scaffold":
                            self.assemble.add_assemble_seq(assemble_id,
                                                       self.output_dir + '/' + sample + '/assembly_predict/assembly/assembly/' + sample + '_assembly_' + type + '_details.xls',
                                                       sample, type,self.remote_dir, dict12)
                        else:
                            self.assemble.add_assemble_seq(assemble_id,
                                                           self.output_dir + '/' + sample + '/assembly_predict/assembly/assembly/' + sample + '_assembly_' + type + '_details.xls', sample, type, self.remote_dir)
                        for window in ['1k', '2k', '5k']:
                            if window == '1k':
                                self.assemble.add_assemble_graphic(assemble_id,
                                                                   self.output_dir + '/' + sample + '/assembly_predict/assembly/len/' + sample + '.1000.' + type + 's.len.xls',
                                                                   sample, type, '1k')
                            elif window == '2k':
                                self.assemble.add_assemble_graphic(assemble_id,
                                                                   self.output_dir + '/' + sample + '/assembly_predict/assembly/len/' + sample + '.2000.' + type + 's.len.xls',
                                                                   sample, type, '2k')
                            elif window == '5k':
                                self.assemble.add_assemble_graphic(assemble_id,
                                                                   self.output_dir + '/' + sample + '/assembly_predict/assembly/len/' + sample + '.5000.' + type + 's.len.xls',
                                                                   sample, type, '5k')
                elif self.option("analysis") in ["complete"]:
                    dict12 = self.get_seq_dict2(self.output_dir + '/' + sample + '/plasmid_predict/all.stat.xls')
                    if self.option('raw_dir').is_set:
                        self.assemble.add_assemble_stat_complete(assemble_id,
                                                                 self.output_dir + '/' + sample + '/assembly_predict/assembly/' + sample + '_assembly_summary.xls',
                                                                 sample,file2 =self.output_dir + '/' + sample + '/assembly_predict/assembly/busco/' + sample + '_busco.xls')
                        self.assemble.add_assemble_complete_seq(assemble_id,
                                                                self.output_dir + '/' + sample + '/assembly_predict/assembly/' + sample + '_assembly_details.xls',
                                                                sample, self.remote_dir, dict12)
                    else:
                        self.assemble.add_assemble_stat_complete(assemble_id,
                                                                 self.output_dir + '/' + sample + '/assembly_predict/assembly/' + sample + '_assembly_summary.xls',
                                                                 sample)
                        self.assemble.add_assemble_complete_seq(assemble_id,
                                                                self.output_dir + '/' + sample + '/assembly_predict/assembly/' + sample + '_assembly_details.xls',
                                                                sample, self.remote_dir, dict12)

                ## 预测
                if self.option("analysis") in ["uncomplete"]:
                    self.gene_predict.add_gene_predict_detail(gene_id, sample,
                                                              self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS.gff")
                    self.gene_predict.add_gene_predict_specimen(gene_id, sample,
                                                                self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS_statistics.xls")
                    self.gene_predict.add_gene_predict_seq(gene_id, sample,
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS.fnn",
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_CDS.faa")
                    self.gene_predict.add_gene_predict_bar(gene_id, sample,
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/' + "length_distribute.txt")
                elif self.option("analysis") in ["complete"]:
                    self.gene_predict.add_gene_predict_detail(gene_id, sample,
                                                              self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS.gff")
                    self.gene_predict.add_gene_predict_specimen(gene_id, sample,
                                                                self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS_statistics.xls")
                    self.gene_predict.add_gene_predict_seq(gene_id, sample,
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS.fnn",
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/CDS_predict/' + sample + "_whole_genome_CDS.faa")
                    self.gene_predict.add_gene_predict_bar(gene_id, sample,
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/' + "length_distribute.txt")
                self.trna_predict.add_trna_predict_detail(trna_id, sample,
                                                          self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.gff")
                if os.path.exists(
                                                                self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.fnn"):
                    self.trna_predict.add_trna_predict_seq(trna_id, sample,
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.fnn",
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/tRNA/' + sample + "_tRNA.struc")
                if os.path.exists(
                                                                self.output_dir + '/' + sample + '/assembly_predict/predict/rRNA/' + sample + "_rRNA.fnn"):
                    self.rrna_predict.add_rrna_predict_seq(rrna_id, sample,
                                                           self.output_dir + '/' + sample + '/assembly_predict/predict/rRNA/' + sample + "_rRNA.fnn")
                self.rrna_predict.add_rrna_predict_detail(rrna_id, sample,
                                                          self.output_dir + '/' + sample + '/assembly_predict/predict/rRNA/' + sample + "_rRNA.gff")
                if self.option("analysis") in ["uncomplete"]:
                    self.repeat_predict.add_repeat_predict_detail(repeat_id, sample,
                                                                  self.output_dir + '/' + sample + '/assembly_predict/predict/repeats/' + sample + "_TRF.gff")
                elif self.option("analysis") in ["complete"]:
                    self.repeat_predict.add_repeat_predict_detail(repeat_id, sample,
                                                                  self.output_dir + '/' + sample + '/assembly_predict/predict/repeats/' + sample + "_whole_genome_TRF.gff")

            ### 原流程和注释流程共用
            if self.option("analysis") in ["uncomplete"]:
                self.anno_nr.add_anno_nr_detail(nr_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/NR/" + sample + "_anno_nr.xls")
            elif self.option("analysis") in ["complete"]:
                self.anno_nr.add_anno_nr_detail(nr_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/NR/" + sample + "_whole_genome_anno_nr.xls")
            self.anno_cog.add_anno_cog_detail(cog_id, sample,
                                              self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls")
            self.anno_go.add_anno_go_specimen(go_id, sample,
                                              self.output_dir + '/' + sample + "/annotation/GO/" + sample + "_go_statistics.xls") ### 去掉！！！_go_level2_statistics.xls
            self.anno_go.add_anno_go_detail(go_id, sample,
                                              self.output_dir + '/' + sample + "/annotation/GO/" + sample + "_go_anno.xls")
            self.anno_cazy.add_anno_cazy_detail(cazy_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/CAZy/" + sample + "_anno_cazy.xls")
            self.anno_kegg.add_anno_kegg_detail(kegg_id, sample,
                                                self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls")
            self.anno_kegg.add_anno_kegg_pic(kegg_id, sample, self.output_dir + '/' + sample + "/annotation/KEGG/{}_kegg_pathway_img".format(sample))

            self.anno_kegg.add_anno_kegg_level(kegg_id, sample,
                                               self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_level_stat.xls")


            if self.option("analysis") in ["uncomplete"]:
                self.anno_pfam.add_anno_pfam_detail(pfam_id, sample,
                                                    self.output_dir + '/' + sample + "/annotation/Pfam/" + sample + "_anno_pfam.xls")
            elif self.option("analysis") in ["complete"]:
                self.anno_pfam.add_anno_pfam_detail(pfam_id, sample,
                                                    self.output_dir + '/' + sample + "/annotation/Pfam/" + sample + "_whole_genome_anno_pfam.xls")
            if self.option("analysis") in ["uncomplete"]:
                self.anno_swissprot.add_anno_swissprot_detail(swissport_id, sample,
                                                              self.output_dir + '/' + sample + "/annotation/Swissprot/" + sample + "_anno_swissprot.xls")
            elif self.option("analysis") in ["complete"]:
                self.anno_swissprot.add_anno_swissprot_detail(swissport_id, sample,
                                                              self.output_dir + '/' + sample + "/annotation/Swissprot/" + sample + "_whole_genome_anno_swissprot.xls")
            if not self.is_anno_pip:
                dict13= self.get_seq_dict3(self.output_dir + '/' + sample + '/plasmid_predict/all.stat.xls')
            else:
                dict13 = {}
            if os.path.exists(self.output_dir + '/' + sample + "/" + sample + "_summary.xls"):
                self.anno_summary.add_anno_summary_detail(summary_id, sample, self.output_dir + '/' + sample + "/" + sample + "_summary.xls",dict13)

            if self.option("analysis") in ["complete"]:
                if os.path.exists(self.output_dir + '/' + sample + "/metabolic_system/antiSMASH"):
                    files = os.listdir(self.output_dir + '/' + sample + "/metabolic_system/antiSMASH")
                    if len(files) != 0:
                        for file in files:
                            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            if file.endswith('antismash_anno.xls'):
                                self.anno_antismash.add_anno_antismash_stat(antismash_id, sample,self.output_dir + '/' + sample + '/metabolic_system/antiSMASH/' + file)

            elif self.option("analysis") in ["uncomplete"]:
                if os.path.exists(self.output_dir + '/' + sample + "/metabolic_system/antiSMASH/antismash_anno.xls"):
                    self.anno_antismash.add_anno_antismash_stat(antismash_id, sample,self.output_dir + '/' + sample + "/metabolic_system/antiSMASH/antismash_anno.xls")
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_summary.xls"):
                self.island.add_island_detail(island_id,
                                              self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_summary.xls",sample)


            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  add 20190429
            stat_path = self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/sample_stat.xls'
            if os.path.exists(stat_path):
                island_stat_table = pd.read_table(stat_path,sep="\t",header=0)
                if len(island_stat_table) >0:
                    self.island.add_island_stat(island_id, stat_path, sample)

            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_detail.xls"):
                self.island.add_island_gene(island_id,
                                            self.output_dir + '/' + sample + '/mobile_elements/Genomic_Islands/' + sample + "_GI_detail.xls",
                                            sample)
            if os.path.exists(self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + ".stat.xls"):
                self.prephage.add_prephage_stat(prephage_id,
                                                self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + ".stat.xls",
                                                sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_summary.xls"):
                self.prephage.add_prephage_detail(prephage_id,
                                                  self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_summary.xls",
                                                  sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_detail.xls"):
                self.prephage.add_prephage_gene(prephage_id,
                                                self.output_dir + '/' + sample + '/mobile_elements/prephage/' + sample + "_prephage_detail.xls",
                                                sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_summary.xls"):
                self.crispr.add_crispr_stat(crispr_id,
                                            self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_summary.xls",
                                            sample)
                self.crispr.add_crispr_detail(crispr_id,
                                              self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_summary.xls",
                                              sample)
            if os.path.exists(
                                                            self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_detail.xls"):
                self.crispr.add_crispr_psa(crispr_id,
                                           self.output_dir + '/' + sample + '/mobile_elements/CRISPR_Cas/' + sample + "_CRISPR_Cas_detail.xls",
                                           sample)

            ## is预测
            self.logger.info("start is>>>>>>>>>>>>")
            sample_dir = self.output_dir + '/' + sample + '/mobile_elements/Is_Predict'
            sample_path = os.path.join(sample_dir, sample + ".stat.xls")
            if os.path.exists(sample_path):
                self.api_is.add_anno_is_stat(sample_path, main_id=is_id, sample=sample)
            stat_file = os.path.join(sample_dir, sample + ".raw")
            seq_file = os.path.join(sample_dir, sample + "_sequence.fna")
            gff_file = os.path.join(sample_dir, sample + ".is.xls")
            enzyme_file = os.path.join(sample_dir, sample + ".enzyme.xls")
            gene_file = os.path.join(sample_dir, sample + ".gene.ffn")
            blast_file = os.path.join(sample_dir, sample + ".blast.xls")
            if os.path.exists(gff_file) and os.path.getsize(gff_file) > 0:
                if os.path.exists(enzyme_file) and os.path.exists(gene_file):
                    self.api_is.add_anno_is_detail(stat_file,seq_file, main_id=is_id, sample=sample, gff_file=gff_file, transposon_file=enzyme_file,blast_file=blast_file)

                    self.api_is.add_anno_is_summary(gff_file, seq_file, sample, main_id=is_id,transposon_file=enzyme_file, gene_file=gene_file)
                else:
                    self.api_is.add_anno_is_detail(stat_file,seq_file, main_id=is_id, sample=sample,blast_file=blast_file)
                    self.api_is.add_anno_is_summary(gff_file, seq_file, sample, main_id=is_id)
            self.logger.info("end is>>>>>>>>>>>>")

            ## integron预测
            self.logger.info("start integron>>>>>>>>>>>>")
            sample_dir = self.output_dir + '/' + sample + '/mobile_elements/Integron'
            sample_path = os.path.join(sample_dir, sample + ".sample.xls")
            stat_file = os.path.join(sample_dir, sample + ".stat.xls")
            seq_file = os.path.join(sample_dir, sample + ".integron.fna")
            if os.path.exists(sample_path):
                self.api_integron.add_anno_integron_stat(sample_path, main_id=integron_id)
            if os.path.exists(stat_file) and os.path.exists(seq_file):
                self.api_integron.add_anno_integron_detail(stat_file,seq_file, main_id=integron_id)
            summary_path = os.path.join(sample_dir, sample + ".integrons")
            sequence_path = os.path.join(sample_dir, sample + "_sequence.fna")
            if os.path.exists(summary_path) and os.path.exists(sequence_path):
                self.api_integron.add_anno_integron_summary(summary_path, sequence_path, sample, main_id=integron_id)
            self.logger.info("end integron>>>>>>>>>>>>")

            ##散在重复序列预测
            self.logger.info("start interpersed repeat>>>>>>>>>>>>")
            sample_dir = os.path.join(self.output_dir, sample, "assembly_predict/predict/Interpersed_repeat")
            sample_path = os.path.join(sample_dir, sample + ".gff")
            if not self.is_anno_pip:
                scaffoled_path = self.output_dir+"/"+ sample+"/"+ sample+ ".all.fna"
            else:
                scaffoled_path = self.output_dir + "/" + sample + "/assembly_predict/predict/Interpersed_repeat/all.fasta"
            if os.path.exists(sample_path):
                self.up_repeatmasker.add_repeat_predict_detail(interpersed_id, sample, sample_path,scaffoled_path)
            self.logger.info("end interpersed repeat>>>>>>>>>>>>")

            ## 耐药基因Resfinder预测
            self.logger.info(sample)
            sample_dir = os.path.join(self.output_dir, sample, "pathogenic_system/Resfinder")
            self.api_resfinder.add_resfinder_dir(sample_dir, sample, main_id=res_id)

            ##新增双组分调控系统导表  zouguanqing
            mongo_key1 = 'gene_id,type,pfam_id,domain,domain_desc,,,,,,,location,'
            two_com1 = self.output_dir + '/'+sample + '/annotation/Two_component/'+ sample+'.senser_regulator.xls'
            two_com2 = self.output_dir + '/' + sample + '/annotation/Two_component/'+sample+ '.senser_regulator.stat'
            self.two_component_api.add_main_detail(two_com1,'anno_regulator_detail', two_component_id, mongo_key1, has_head =True,
                                                   main_name='regulator_id',other_dic={'specimen_id':sample})
            mongo_key2 = 'senser,regulator,hybrid'
            self.two_component_api.add_main_detail(two_com2,'anno_regulator_stat', two_component_id, mongo_key2, has_head =True,
                        main_name='regulator_id',other_dic={'specimen_id':sample},main_table='anno_regulator', update_dic={"main_id":two_component_id})


            if n == 1:
                prom_id = self.promote.add_promote(
                    self.output_dir + '/' + sample + "/structral_genome/promoter_predict",
                    main=True, specimen_id=sample, update_id=summary_id)
                self.logger.info("promoter end")
                self.api_gene_graph = self.api.api("bacgenome.gene_graph")
                if self.option("analysis") in ["uncomplete"]:
                    if not self.is_anno_pip:
                        cds_fnn = self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_CDS.fnn"
                    else:
                        cds_fnn = self.output_dir + '/' + sample + '/' + sample + '.cds.ffn'
                    gene_graph_id = self.gene_graph.add_gene_graph(
                        cds_fnn,
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls", main=True,
                        specimen_id=sample)
                    self.logger.info("gene_graph end")
                elif self.option("analysis") in ["complete"]:
                    if not self.is_anno_pip:
                        cds_fnn = self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_whole_genome_CDS.fnn"
                    else:
                        cds_fnn = self.output_dir + '/' + sample + '/' + sample + '.cds.ffn'

                    gene_graph_id = self.gene_graph.add_gene_graph(
                        cds_fnn,
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls",
                        main=True, specimen_id=sample)
                    self.logger.info("gene_graph end")
                vfdb_id = self.anno_vfdb.add_anno_vfdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/VFDB", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("vfdb end")
                phi_id = self.anno_phi.add_anno_phi(
                    self.output_dir + '/' + sample + "/pathogenic_system/PHI", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("phi end")
                tcdb_id = self.anno_tcdb.add_anno_tcdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/TCDB", main=True, specimen_id=sample,
                    update_id=summary_id,analysis=self.option('analysis'))
                self.logger.info("tcdb end")
                self.api_card = self.api.api("bacgenome.anno_card")
                card_id = self.anno_card.add_anno_card(
                    self.output_dir + '/' + sample + "/pathogenic_system/CARD", main=True, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("card end")
                tmhmm_id = self.anno_tmhmm.add_anno_tmhmm(
                    self.output_dir + '/' + sample + "/pathogenic_system/TMHMM", main=True, specimen_id=sample,
                    update_id=summary_id,analysis=self.option('analysis'))
                self.logger.info("tmhmm end")
                secretory_id = self.secretory.add_secretory(
                    self.output_dir + '/' + sample + "/pathogenic_system/secretion_system",
                    anno_summary_id=summary_id, kegg_id=kegg_id, main=True, specimen_id=sample)
                n +=1
            else:
                self.logger.info(prom_id)
                self.promote.add_promote(
                    self.output_dir + '/' + sample + "/structral_genome/promoter_predict",
                    main_id=prom_id, specimen_id=sample, update_id=summary_id)
                self.logger.info("promoter end")
                if self.option("analysis") in ["uncomplete"]:
                    if not self.is_anno_pip:
                        cds_fnn = self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_CDS.fnn"
                    else:
                        cds_fnn = self.output_dir + '/' + sample + '/' + sample + '.cds.ffn'

                    self.gene_graph.add_gene_graph(
                        cds_fnn,
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls", main_id=gene_graph_id,
                        specimen_id= sample)
                    self.logger.info("gene_graph end")
                elif self.option("analysis") in ["complete"]:
                    if not self.is_anno_pip:
                        cds_fnn = self.output_dir + '/' + sample + "/assembly_predict/predict/CDS_predict/" + sample + "_whole_genome_CDS.fnn"
                    else:
                        cds_fnn = self.output_dir + '/' + sample + '/' + sample + '.cds.ffn'
                    self.gene_graph.add_gene_graph(
                        cds_fnn,
                        self.output_dir + '/' + sample + "/annotation/KEGG/" + sample + "_kegg_anno.xls",
                        self.output_dir + '/' + sample + "/annotation/COG/" + sample + "_cog_anno.xls",
                        main_id=gene_graph_id, specimen_id=sample)
                    self.logger.info("gene_graph end")
                self.anno_vfdb.add_anno_vfdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/VFDB", main_id=vfdb_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("vfdb end")
                self.anno_phi.add_anno_phi(
                    self.output_dir + '/' + sample + "/pathogenic_system/PHI", main_id=phi_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("phi end")
                self.anno_tcdb.add_anno_tcdb(
                    self.output_dir + '/' + sample + "/pathogenic_system/TCDB", main_id=tcdb_id, specimen_id=sample,
                    update_id=summary_id, analysis=self.option('analysis'))
                self.logger.info("tcdb end")
                self.anno_card.add_anno_card(
                    self.output_dir + '/' + sample + "/pathogenic_system/CARD", main_id=card_id, specimen_id=sample,
                    update_id=summary_id)
                self.logger.info("card end")
                self.anno_tmhmm.add_anno_tmhmm(
                    self.output_dir + '/' + sample + "/pathogenic_system/TMHMM", main_id=tmhmm_id, specimen_id=sample,
                    update_id=summary_id, analysis=self.option('analysis'))
                self.logger.info("tmhmm end")
                self.secretory.add_secretory(self.output_dir + '/' + sample + "/pathogenic_system/secretion_system",
                                             anno_summary_id=summary_id, kegg_id=kegg_id, main_id=secretory_id,
                                             specimen_id=sample)
        self.logger.info("secretory end")
        self.api_circos_table = self.api.api("bacgenome.circos_table")
        self.api_circos_table.add_circos_table(project_sn=project_sn, task_id=task_id)
        self.logger.info("circos_table end")
        # self.summary_map = self.api.api("bacgenome.summary_map")
        # self.summary_map.add_map(task_id=task_id)
        ## circos
        self.circos = self.api.api("bacgenome.circos")
        if self.option("analysis") in ["complete"]:
            if not self.is_anno_pip:
                for sample in self.samples.keys():
                    seq = self.get_seq_type(self.output_dir + "/" + sample +"/"+sample+".all.fna", sample)
                    for key, value in seq.iteritems():
                        locatin_list = value
                        for location in locatin_list.split(","):
                            params = {'specimen_id': key, 'task_type': 2, 'submit_location': 'circos',
                                      'task_id': task_id, 'seq_type': 'Circular',
                                      'location': location, 'para1': 'ncrna'}  # by zzg 圈图反选
                            name = key + '_Circos_' + '_' + location + '_' + datetime.datetime.now().strftime(
                                "%Y%m%d_%H%M%S%f")[:-3]
                            self.circos.add_circos(
                                self.output_dir + '/' + key + '/circos/' + location, project_sn=project_sn, name=name,
                                location=location, task_id=task_id, main=True, specimen_id=key, params=params,
                                link_path=self.remote_dir + key + '/circos/' + location + '/')

            else:
                #20190619
                for sample in self.samples.keys():
                    locs = os.listdir(self.output_dir + '/' + sample + '/circos/')

                    for location in locs:
                        location_ori = location
                        for ori_seq_name in self.samples_seq_map[sample]['map']:
                            if self.samples_seq_map[sample]['map'][ori_seq_name] == location:
                                location_ori = ori_seq_name
                                break
                        params = {'specimen_id': sample, 'task_type': 2, 'submit_location': 'circos',
                                      'task_id': task_id, 'seq_type':'Circular',
                                      'location': location_ori, 'para1': 'ncrna'} # by zzg 圈图反选
                        name = sample + '_Circos_' + '_' + location_ori + '_' + datetime.datetime.now().strftime(
                            "%Y%m%d_%H%M%S%f")[:-3]
                        self.circos.add_circos(
                            self.output_dir + '/' + sample + '/circos/' + location, project_sn=project_sn, name=name,
                            location=location_ori, task_id=task_id, main=True, specimen_id=sample, params=params,
                            link_path=self.remote_dir + sample + '/circos/'+ location + '/')


        if self.option("analysis") in ["uncomplete"]:
            for sample in self.samples.keys():
                params = {'specimen_id': sample, 'seq_type':'Circular', 'task_type': 2, 'submit_location': 'circos', 'task_id': task_id, 'para1': 'ncrna'} # by zzg 圈图反选
                name = sample + '_Circos_' + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
                self.circos.add_circos(
                    self.output_dir + '/' + sample + '/circos/Scaffold', project_sn=project_sn, name=name,
                    task_id=task_id, main=True, specimen_id=sample,
                    params=params, link_path=self.remote_dir + sample + '/circos/Scaffold/')
        ## cgview
        self.cgview = self.api.api("bacgenome.cgview")
        if self.option("analysis") in ["uncomplete"]:
            for sample in self.samples.keys():
                params = {'genome_type': 'Scaffold', 'species_name': '', 'specimen_id': sample,
                          "submit_location": 'cgview', 'seq_type':'Circular',
                          "task_type": 2}
                file = self.output_dir + '/' + sample + '/cgview'
                name = sample + '_Cgview_' + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
                self.cgview.add_cgview(file, params=params, project_sn=project_sn,specimen_id=sample,genome_type= 'scaffold',
                                       task_id=task_id, main='true', name=name, link_path=self.remote_dir + sample + '/cgview/')
        if self.option("analysis") in ["complete"]:
            for sample in self.samples.keys():
                file = self.output_dir + '/' + sample + '/cgview'
                files = os.listdir(file)
                for fi in files:
                    name = fi
                    fil = self.output_dir + '/' + sample + '/cgview/' + fi
                    params = {'genome_type': name, 'species_name': '', 'specimen_id': sample,
                              "submit_location": 'cgview', 'seq_type':'Circular',
                              "task_type": 2}
                    name2 = sample + '_Cgview_' + '_' + fi + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[
                                                                   :-3]
                    self.cgview.add_cgview(fil, genome_type=name, specimen_id=sample, name=name2, params=params,
                                           project_sn=project_sn,
                                           task_id=task_id, main='true',
                                           link_path=self.remote_dir + sample + '/cgview/' + fi + '/')
        ## 核心基因进化树
        if self.option("analysis") in ["uncomplete"]:
            if os.path.exists(self.hgene_tree.work_dir + '/phylo_tree.nwk'):
                    self.datastat.add_datastat_tree_detail(datastat_id,self.hgene_tree.work_dir + '/phylo_tree.nwk',
                                                       '', 'hgene', 'all')
        if self.option("analysis") in ["complete"]:
            if os.path.exists(self.hgene_tree.work_dir + '/phylo_tree.nwk'):
                self.datastat.add_datastat_tree_detail(datastat_id, self.hgene_tree.work_dir + '/phylo_tree.nwk',
                                                       '', 'hgene', 'all')
            if os.path.exists(self.tree.work_dir + '/phylo_tree.nwk'):
                self.datastat.add_datastat_tree_detail(datastat_id,self.tree.work_dir + '/phylo_tree.nwk',
                                                       '', '16s', 'all')
        ## 基因组质量评估
        genome_quality = self.api.api('bacgenome.genome_quality')
        genome_quality.add_genome_quality_detail(task_id, self.option("analysis"),project_sn)


    def re_location(self,location_list):
        rename = []
        for one in location_list.split(','):
            name = one
            rename.append(name)
        return rename


    def send_files(self):
        """
        结果放置到/upload_results
        """
        dir_o = self.output_dir
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        repaths = []
        regexps = []
        for sample in self.samples.keys():
            if self.is_anno_pip:
                self.move_dir(os.path.join(dir_o, sample +'/gene'),os.path.join(dir_up, sample + "/gene"))

            if not self.is_anno_pip:
                if self.option("analysis") in ["uncomplete"]:
                    files = os.listdir(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly/"))
                    for file in files:
                        if re.search(r'.fna.index.', file):
                            os.remove(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly/" + file))
                if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                    if os.path.exists(os.path.join(dir_o, sample + "/data_QC")):
                        self.move_dir(os.path.join(dir_o, sample + "/data_QC"), os.path.join(dir_up, sample + "/data_QC"))
                elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                    if os.path.exists(os.path.join(dir_o, sample + "/data_QC")):
                        self.move_dir(os.path.join(dir_o, sample + "/data_QC"), os.path.join(dir_up, sample + "/data_QC"))
                if self.option("analysis") in ["uncomplete"] and self.option('raw_dir').is_set:
                    if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency")):
                        self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency"),os.path.join(dir_up, sample + "/genomic_assessment/kmer_frequency"))
                        for i in [1000,3000,5000,8000,10000]:
                            if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/depth_gc_" + str(i))):
                                self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/depth_gc_" + str(i)),os.path.join(dir_up, sample + "/genomic_assessment/depth_gc_" + str(i)))
                elif self.option("analysis") in ["complete"] and self.option('raw_dir').is_set:
                    if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency")):
                        self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/kmer_frequency"),os.path.join(dir_up, sample + "/genomic_assessment/kmer_frequency"))
                        for i in [1000,3000,5000,8000,10000]:
                            if os.path.exists(os.path.join(dir_o, sample + "/genomic_assessment/depth_gc_" + str(i))):
                                self.move_dir(os.path.join(dir_o, sample + "/genomic_assessment/depth_gc_" + str(i)),os.path.join(dir_up, sample + "/genomic_assessment/depth_gc_" + str(i)))
                if self.option("analysis") in ["uncomplete"]:
                    if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly")):
                        self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/assembly"),os.path.join(dir_up, sample + "/assembly_predict/assembly"))
                elif self.option("analysis") in ["complete"]:
                    if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/")):
                        self.move_file(os.path.join(dir_o, sample + "/assembly_predict/assembly/" + sample + '_assembly_summary.xls'),os.path.join(dir_up, sample + "/assembly_predict/assembly/"  + sample + '_assembly_summary.xls'))
                        self.move_file(os.path.join(dir_o, sample + "/assembly_predict/assembly/" + sample + '_assembly_details.xls'),os.path.join(dir_up, sample + "/assembly_predict/assembly/" + sample + '_assembly_details.xls'))
                    if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/assembly/seq_dir")):
                        self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/assembly/seq_dir"),os.path.join(dir_up, sample + "/assembly_predict/assembly/seq_dir"))
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/CDS_predict")):
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/CDS_predict"),os.path.join(dir_up, sample + "/assembly_predict/predict/CDS_predict"))
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/tRNA")):
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/tRNA"),os.path.join(dir_up, sample + "/assembly_predict/predict/tRNA"))
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/rRNA")):
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/rRNA"),os.path.join(dir_up, sample + "/assembly_predict/predict/rRNA"))
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/repeats")):
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/repeats"),os.path.join(dir_up, sample + "/assembly_predict/predict/repeats"))
                if os.path.exists(os.path.join(dir_o, sample + "/assembly_predict/predict/Interpersed_repeat")):
                    self.move_dir(os.path.join(dir_o, sample + "/assembly_predict/predict/Interpersed_repeat"),os.path.join(dir_up, sample + "/assembly_predict/predict/Interpersed_repeat"))

                if self.option("analysis") in ["complete"]:
                    if os.path.exists(os.path.join(dir_o, sample + "/methylation")):
                        self.move_dir(os.path.join(dir_o, sample + "/methylation"), os.path.join(dir_up, sample + "/methylation"))
                    if os.path.exists(os.path.join(dir_up, sample + "/project_temp")):
                        shutil.rmtree(os.path.join(dir_up, sample + "/project_temp"))
                    os.mkdir(os.path.join(dir_up, sample + "/project_temp"))
                    if os.path.exists(os.path.join(dir_o, sample + "/Unicycler", sample +".raw.assembly.log")):
                        self.move_file(os.path.join(dir_o, sample + "/Unicycler", sample +".raw.assembly.log"), os.path.join(dir_up, sample + "/project_temp", sample +".raw.assembly.log"))
                    if os.path.exists(os.path.join(dir_o, sample + "/Unicycler", sample + ".scaffold.fna")):
                        self.move_file(os.path.join(dir_o, sample + "/Unicycler", sample + ".scaffold.fna"),
                                           os.path.join(dir_up, sample + "/project_temp", sample + ".raw.assembly.fa"))
                if os.path.exists(os.path.join(dir_up, sample + "/assembly_predict/Plasmid_identification")):
                    shutil.rmtree(os.path.join(dir_up, sample + "/assembly_predict/Plasmid_identification"))
                os.makedirs(os.path.join(dir_up, sample + "/assembly_predict/Plasmid_identification"))
                if os.path.exists(os.path.join(dir_o, sample + "/plasmid_predict", sample +".detail.xls")):
                    self.move_file(os.path.join(dir_o, sample + "/plasmid_predict", sample +".detail.xls"), os.path.join(dir_up, sample + "/assembly_predict/Plasmid_identification", sample +".detail.xls"))
                if os.path.exists(os.path.join(dir_o, sample + "/plasmid_predict", sample +".plasmid_anno.xls")):
                    self.move_file(os.path.join(dir_o, sample + "/plasmid_predict", sample +".plasmid_anno.xls"), os.path.join(dir_up, sample + "/assembly_predict/Plasmid_identification", sample +".plasmid_anno.xls"))
                if os.path.exists(os.path.join(dir_o, sample + "/plasmid_predict", "all.stat.xls")):
                    self.move_file(os.path.join(dir_o, sample + "/plasmid_predict", "all.stat.xls"),
                               os.path.join(dir_up, sample + "/assembly_predict/Plasmid_identification",
                                            "all.stat.xls"))
                if os.path.exists(os.path.join(dir_up, sample + "/assembly_predict/predict/sRNA_predict")):
                    shutil.rmtree(os.path.join(dir_up, sample + "/assembly_predict/predict/sRNA_predict"))
                os.makedirs(os.path.join(dir_up, sample + "/assembly_predict/predict/sRNA_predict"))
                if os.path.exists(os.path.join(dir_o, sample + "/srna_predict", sample +".sRNA.xls")):
                    self.move_file(os.path.join(dir_o, sample + "/srna_predict", sample +".sRNA.xls"), os.path.join(dir_up, sample + "/assembly_predict/predict/sRNA_predict", sample +".sRNA.xls"))

            ##ori和anno流程共用
            if self.option("analysis") in ["uncomplete"]:
                if os.path.exists(os.path.join(dir_o, sample + "/gbk")):
                    self.move_dir(os.path.join(dir_o, sample + "/gbk/gbk/Scaffold"), os.path.join(dir_up, sample + "/project_overview/Analysis_files"))
                    self.move_file(os.path.join(dir_o, sample + "/gbk/" +sample + ".ptt"), os.path.join(dir_up, sample + "/project_overview/Analysis_files/"+sample+".ptt"))
                    self.move_file(os.path.join(dir_o, sample + "/gbk/" +sample + ".gff"), os.path.join(dir_up, sample + "/project_overview/Analysis_files/"+sample+".gff"))
            elif self.option("analysis") in ["complete"]:
                if os.path.exists(os.path.join(dir_up, sample + "/project_overview/Analysis_files")):
                    shutil.rmtree(os.path.join(dir_up, sample + "/project_overview/Analysis_files"))
                os.makedirs(os.path.join(dir_up, sample + "/project_overview/Analysis_files"))
                if os.path.exists(os.path.join(dir_o, sample + "/gbk")):
                    ptts = glob.glob(os.path.join(dir_o, sample +"/gbk/*.ptt"))
                    gffs = glob.glob(os.path.join(dir_o, sample +"/gbk/*.gff"))
                    for p in ptts:
                        f_name = os.path.basename(p)
                        self.move_file(p, os.path.join(dir_up, sample + "/project_overview/Analysis_files/"+f_name))
                    for g in gffs:
                        f_name = os.path.basename(g)
                        self.move_file(g, os.path.join(dir_up, sample + "/project_overview/Analysis_files/"+f_name))
                    files = glob.glob(os.path.join(dir_o, sample + "/gbk/seq_gbk/*/*.gbk"))
                    for file in files:
                        f_name = os.path.basename(file)
                        self.move_file(file, os.path.join(dir_up, sample + "/project_overview/Analysis_files/"  + f_name))

            ## send  tree  core.fa
            core_gene_fa = os.path.join(dir_o, sample+'/tree/hgene/'+sample+'.cor_gene.fa')
            tar_core_gene_fa = os.path.join(dir_up, sample + "/project_overview/core_faa/" +sample+'.cor_gene.fa')
            if os.path.exists(core_gene_fa):
                self.move_file(core_gene_fa, tar_core_gene_fa)
            core_gene_fa1 = os.path.join(dir_o, sample+'/tree/hgene/'+sample+'.cor_gene.faa')
            tar_core_gene_fa1 = os.path.join(dir_up, sample + "/project_overview/core_faa/" +sample+'.cor_gene.faa')
            if os.path.exists(core_gene_fa1):
                self.move_file(core_gene_fa1, tar_core_gene_fa1)
            core_gene_fa2 = os.path.join(dir_o, sample + '/tree/hgene/' + sample + '.cor_gene.fnn')
            tar_core_gene_fa2 = os.path.join(dir_up, sample + "/project_overview/core_faa/" + sample + '.cor_gene.fnn')
            if os.path.exists(core_gene_fa2):
                self.move_file(core_gene_fa2, tar_core_gene_fa2)
            core_gene_blast = os.path.join(dir_o, sample + '/tree/hgene/' + sample + '.cor_blast.xls')
            tar_core_gene_blast = os.path.join(dir_up, sample + "/project_overview/core_faa/" + sample + '.cor_gene.xls')
            if os.path.exists(core_gene_blast):
                self.move_file(core_gene_blast, tar_core_gene_blast)
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Swissprot")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Swissprot"),os.path.join(dir_up, sample + "/annotation/Swissprot"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Summary")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Summary"),os.path.join(dir_up, sample + "/annotation/Summary"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/Pfam")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/Pfam"),os.path.join(dir_up, sample + "/annotation/Pfam"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/KEGG")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/KEGG"),os.path.join(dir_up, sample + "/annotation/KEGG"))
            files = os.listdir(os.path.join(dir_up, sample + "/annotation/KEGG/" + sample + "_kegg_pathway_img"))
            for file in files:
                if re.search(r'.png$',file) or re.search(r'.html$',file):
                    os.remove(os.path.join(dir_up, sample + "/annotation/KEGG/" + sample + "_kegg_pathway_img/" + file))

            kegg_pathway_dir = os.path.join(dir_up, sample + "/annotation/KEGG/" + sample + "_kegg_pathway_img")
            if os.path.exists(kegg_pathway_dir):
                if os.path.exists(kegg_pathway_dir+'.tar.gz'):
                    os.remove(kegg_pathway_dir+'.tar.gz')
                os.system('tar -zcf %s.tar.gz %s'%(kegg_pathway_dir, kegg_pathway_dir))
                shutil.rmtree(kegg_pathway_dir)
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/GO")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/GO"),os.path.join(dir_up, sample + "/annotation/GO"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/NR")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/NR"),os.path.join(dir_up, sample + "/annotation/NR"))  #guanqing.zou 20180903
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/COG")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/COG"),os.path.join(dir_up, sample + "/annotation/COG"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/CAZy")):
                self.move_dir(os.path.join(dir_o, sample + "/annotation/CAZy"),os.path.join(dir_up, sample + "/metabolic_system/CAZy"))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/gene_cazy_family_stat.xls")):
                os.link(os.path.join(dir_o, sample + "/annotation/gene_cazy_family_stat.xls"),
                              os.path.join(dir_up, sample + "/metabolic_system/CAZy/" + sample + '_cazy_family_stat.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/annotation/gene_cazy_class_stat.xls")):
                os.link(os.path.join(dir_o, sample + "/annotation/gene_cazy_class_stat.xls"),
                              os.path.join(dir_up, sample + "/metabolic_system/CAZy/" + sample + '_cazy_class_stat.xls'))
            if os.path.exists(os.path.join(dir_o, sample + "/structral_genome/promoter_predict")):
                self.move_dir(os.path.join(dir_o, sample + "/structral_genome/promoter_predict"), os.path.join(dir_up, sample + "/structral_genome/promoter_predict"))
            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/Genomic_Islands")):
                self.move_dir(os.path.join(dir_o, sample + "/mobile_elements/Genomic_Islands"), os.path.join(dir_up, sample + "/mobile_elements/Genomic_Islands"))
            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/CRISPR_Cas")):
                self.move_dir(os.path.join(dir_o, sample + "/mobile_elements/CRISPR_Cas"), os.path.join(dir_up, sample + "/mobile_elements/CRISPR_Cas"))
            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/Is_Predict")):
                self.move_dir(os.path.join(dir_o, sample + "/mobile_elements/Is_Predict"), os.path.join(dir_up, sample + "/mobile_elements/Is_Predict"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".gene.ffn")):
                    os.remove(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".gene.ffn"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".raw")):
                    os.remove(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".raw"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".sum")):
                    os.remove(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".sum"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".gff")):
                    os.remove(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".gff"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".orf.faa")):
                    os.remove(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".orf.faa"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".orf.fna")):
                    os.remove(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + ".orf.fna"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + "_sequence.fna")):
                    os.remove(os.path.join(dir_up, sample + "/mobile_elements/Is_Predict", sample + "_sequence.fna"))

            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/Integron")):
                self.move_dir(os.path.join(dir_o, sample + "/mobile_elements/Integron"), os.path.join(dir_up, sample + "/mobile_elements/Integron"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Integron", sample + ".stat.xls")):
                    os.rename(os.path.join(dir_up, sample + "/mobile_elements/Integron", sample + ".stat.xls"), os.path.join(dir_up, sample + "/mobile_elements/Integron", sample + "_Integron_detail.xls"))
                if os.path.exists(os.path.join(dir_up, sample + "/mobile_elements/Integron", sample + ".sample.xls")):
                    os.rename(os.path.join(dir_up, sample + "/mobile_elements/Integron", sample + ".sample.xls"), os.path.join(dir_up, sample + "/mobile_elements/Integron", sample + "_stat.xls"))

            if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prephage")):
                if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prephage/" + sample + '_prephage_summary.xls')):
                    self.move_file(os.path.join(dir_o, sample + "/mobile_elements/prephage/" + sample + '_prephage_summary.xls'),
                      os.path.join(dir_up, sample + "/mobile_elements/prephage/" + sample + '_prephage_summary.xls'))
                if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prephage/" + sample + '_prephage_detail.xls')):
                    self.move_file(os.path.join(dir_o, sample + "/mobile_elements/prephage/" + sample + '_prephage_detail.xls'),
                    os.path.join(dir_up, sample + "/mobile_elements/prephage/" + sample + '_prophage_detail.xls'))
                if os.path.exists(os.path.join(dir_o, sample + "/mobile_elements/prephage/" + sample + '_prephage.fna')):
                    self.move_file(os.path.join(dir_o, sample + "/mobile_elements/prephage/" + sample + '_prephage.fna'),
                    os.path.join(dir_up, sample + "/mobile_elements/prephage/" + sample + '_prephage.fna'))
            if os.path.exists(os.path.join(dir_o, sample + "/metabolic_system/antiSMASH")):
                self.move_dir(os.path.join(dir_o, sample + "/metabolic_system/antiSMASH"), os.path.join(dir_up, sample + "/metabolic_system/antiSMASH"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/VFDB")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/VFDB"), os.path.join(dir_up, sample + "/pathogenic_system/VFDB"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/TMHMM")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/TMHMM"), os.path.join(dir_up, sample + "/pathogenic_system/TMHMM"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/TCDB")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/TCDB"), os.path.join(dir_up, sample + "/pathogenic_system/TCDB"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/PHI")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/PHI"), os.path.join(dir_up, sample + "/pathogenic_system/PHI"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/CARD")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/CARD"), os.path.join(dir_up, sample + "/pathogenic_system/CARD"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/secretion_system")):
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/secretion_system"), os.path.join(dir_up, sample + "/pathogenic_system/secretion_system"))
            if os.path.exists(os.path.join(dir_o, sample + "/pathogenic_system/Resfinder")): ## 新增resfinder软件结果
                self.move_dir(os.path.join(dir_o, sample + "/pathogenic_system/Resfinder"), os.path.join(dir_up, sample + "/pathogenic_system/Resfinder"))
            if os.path.exists(os.path.join(dir_o, sample + "/cgview")):
                self.move_dir(os.path.join(dir_o, sample + "/cgview"), os.path.join(dir_up, sample + "/cgview"))
            if os.path.exists(os.path.join(dir_o, sample + "/circos")):
                self.move_dir(os.path.join(dir_o, sample + "/circos"),os.path.join(dir_up, sample + "/circos"))

            ##生成新的gbk文件，用于结果文件的生成
            if self.option("analysis") in ["uncomplete"]:
                if os.path.exists(os.path.join(dir_o, sample + "/real_gbk")):
                    self.move_dir(os.path.join(dir_o, sample + "/real_gbk/gbk/Scaffold"), os.path.join(dir_up, sample + "/project_overview/All_predict_info"))
                    self.move_file(os.path.join(dir_o, sample + "/real_gbk/" +sample + ".ptt"), os.path.join(dir_up, sample + "/project_overview/All_predict_info/"+sample+".ptt"))
                    self.move_file(os.path.join(dir_o, sample + "/real_gbk/" +sample + ".gff"), os.path.join(dir_up, sample + "/project_overview/All_predict_info/"+sample+".gff"))
            elif self.option("analysis") in ["complete"]:
                if os.path.exists(os.path.join(dir_up, sample + "/project_overview/All_predict_info")):
                    shutil.rmtree(os.path.join(dir_up, sample + "/project_overview/All_predict_info"))
                os.makedirs(os.path.join(dir_up, sample + "/project_overview/All_predict_info"))
                if os.path.exists(os.path.join(dir_o, sample + "/real_gbk")):
                    ptts = glob.glob(os.path.join(dir_o, sample +"/real_gbk/*.ptt"))
                    gffs = glob.glob(os.path.join(dir_o, sample +"/real_gbk/*.gff"))
                    for p in ptts:
                        f_name = os.path.basename(p)
                        self.move_file(p, os.path.join(dir_up, sample + "/project_overview/All_predict_info/"+f_name))
                    for g in gffs:
                        f_name = os.path.basename(g)
                        self.move_file(g, os.path.join(dir_up, sample + "/project_overview/All_predict_info/"+f_name))
                    files = glob.glob(os.path.join(dir_o, sample + "/real_gbk/seq_gbk/*/*.gbk"))
                    for file in files:
                        f_name = os.path.basename(file)
                        self.move_file(file, os.path.join(dir_up, sample + "/project_overview/All_predict_info/"  + f_name))
            genome_path = self.output_dir + '/' + sample + '/' + sample + '.all.fna'
            gff_path = self.output_dir + '/' + sample + '/' + sample + '.all.gff'
            if os.path.exists(genome_path):## 用于打通小工具
                shutil.copyfile(genome_path, os.path.join(dir_up, sample, sample + ".all.fna"))
            if os.path.exists(gff_path):
                shutil.copyfile(gff_path, os.path.join(dir_up, sample, sample+ ".all.gff"))

            if not self.is_anno_pip:
                repaths += [
                    ["%s/data_QC" % sample, "", "样品%s的测序数据质控统计目录" % sample,0,'130004'],
                    ["%s/project_temp" % sample, "", "样品%s的项目中间结果目录" % sample, 0],
                    ["%s/genomic_assessment" % sample, "", "样品%s的基因组评估结果目录" % sample,0,'130005'],
                    ["%s/genomic_assessment/depth_gc_.*" % sample, "", "样品%s的基因组评估gc_depth结果目录" % sample,0,"130134"],
                    ["%s/genomic_assessment/kmer_frequency" % sample, "", "样品%s的基因组kmer评估目录" % sample,0,'130006'],
                    ["%s/assembly_predict" % sample, "", "样品%s的基因组组装与预测结果目录" % sample,0,'130007'],
                    ["%s/assembly_predict/assembly" % sample, "", "样品%s的基因组组装结果目录" % sample,0,'130008'],
                    ["%s/assembly_predict/assembly/seq_dir" % sample, "", "样品%s的基因组完成图seq_dir结果目录" % sample,0,"130135"],
                    ["%s/assembly_predict/predict" % sample, "", "样品%s的基因组基因预测目录" % sample,0,'130009'],
                    ["%s/assembly_predict/predict/CDS_predict" % sample, "", "样品%s的编码基因预测目录" % sample,0,'130010'],
                    ["%s/assembly_predict/predict/repeats" % sample, "", "样品%s的基因组串联重复序列预测结果目录" % sample,0,'130011'],
                    ["%s/assembly_predict/predict/Interpersed_repeat" % sample, "", "样品%s的基因组散在重复序列预测结果目录" % sample,0,"130175"],
                    ["%s/assembly_predict/predict/rRNA" % sample, "", "样品%s的rRNA预测结果目录" % sample,0,'130012'],
                    ["%s/assembly_predict/predict/tRNA" % sample, "", "样品%s的tRNA预测结果目录" % sample,0,'130013'],
                    ["%s/methylation" % sample, "", "样品%s的甲基化结果目录" % sample,0,"130136"],
                    ["%s/assembly_predict/Plasmid_identification" % sample, "", "样品%s的基因组组装结果目录基因组质粒鉴定与注释结果目录" % sample, 0],
                    ["%s/assembly_predict/predict/sRNA_predict" % sample, "", "样品%s的基因组sRNA预测结果目录" % sample, 0],
                ]
            else:
                repaths += [
                    ["%s/gene" % sample, "", "样品%s的基因目录" % sample,0,"130165"] #anno pip add gene

                ]
            ##ori和anno流程共用
            repaths += [
                [".", "", "基础分析结果文件夹", 0, "130000"],
                ["%s" % sample, "", "样品%s流程分析结果目录" % sample,0,'130001'],
                ["%s/project_overview" % sample, "", "样品%s的基因组总览目录" % sample,0,'130002'],
                ["%s/project_overview/Analysis_files" % sample, "", "样品%s的基因组的gbk文件目录" % sample,1,'130003'],
                ["%s/project_overview/All_predict_info" % sample, "", "样品%s的基因组的gbk文件目录" % sample,1,''],
                ["%s/project_overview/core_faa" % sample, "", "样品%s的核心基因序列" % sample,0,"130166"], #core gene fa

                ["%s/annotation" % sample, "", "样品%s的基因注释结果目录" % sample,0,'130014'],
                ["%s/annotation/NR" % sample, "", "样品%s的基因NR注释结果目录" % sample,0,'130015'],
                ["%s/annotation/Swissprot" % sample, "", "样品%s的基因Swiss-Prot注释结果目录" % sample,0,'130016'],
                ["%s/annotation/Pfam" % sample, "", "样品%s的基因Pfam注释结果目录" % sample,0,'130017'],
                ["%s/annotation/COG" % sample, "", "样品%s的基因COG注释结果目录" % sample,0,'130018'],
                ["%s/annotation/GO" % sample, "", "样品%s的基因GO注释结果目录" % sample,0,'130019'],
                ["%s/annotation/KEGG" % sample, "", "样品%s的基因KEGG注释结果目录" % sample,0,'130020'],
                ["%s/annotation/Summary" % sample, "", "样品%s的基因组注释结果汇总目录" % sample,0,'130021'],
                ["%s/annotation/Two_component" % sample, "", "样品%s的基因组注释结果汇总目录" % sample,0,"130167"], #新增Two_component

                ["%s/structral_genome" % sample, "", "样品%s的结构基因组分析目录" % sample,0,'130022'],
                ["%s/structral_genome/promoter_predict" % sample, "", "样品%s的启动子预测结果目录" % sample,0,'130023'],
                ["%s/mobile_elements" % sample, "", "样品%s的可移动元件分析结果目录" % sample,0,'130024'],
                ["%s/mobile_elements/Genomic_Islands" % sample, "", "样品%s的基因组岛预测结果目录" % sample,0,'130025'],
                ["%s/mobile_elements/prephage" % sample, "", "样品%s的前噬菌体预测结果目录" % sample,0,'130026'],
                ["%s/mobile_elements/CRISPR_Cas" % sample, "", "样品%s的CRISPR_Cas系统预测结果目录" % sample,0,'130027'],
                ["%s/mobile_elements/Integron" % sample, "", "样品%s的整合子预测结果目录" % sample,0,"130176"],
                ["%s/mobile_elements/Is_Predict" % sample, "", "样品%s的插入序列预测结果目录" % sample,0,"130177"],
                ["%s/metabolic_system" % sample, "", "样品%s的代谢系统分析目录" % sample,0,'130028'],
                ["%s/metabolic_system/CAZy" % sample, "", "样品%s的碳水化合物活性酶注释结果目录" % sample,0,'130029'],
                ["%s/metabolic_system/antiSMASH" % sample, "", "样品%s的次级代谢产物合成基因簇分析结果目录" % sample,0,'130030'],
                ["%s/pathogenic_system" % sample, "", "样品%s的致病系统分析目录" % sample,0,'130031'],
                ["%s/pathogenic_system/VFDB" % sample, "", "样品%s的毒力基因预测结果目录" % sample,0,'130032'],
                ["%s/pathogenic_system/CARD" % sample, "", "样品%s的耐药基因预测结果目录" % sample,0,'130033'],
                ["%s/pathogenic_system/PHI" % sample, "", "样品%s的病原菌与宿主互作分析结果目录" % sample,0,'130034'],
                ["%s/pathogenic_system/TCDB" % sample, "", "样品%s的转运蛋白分析结果目录" % sample,0,'130035'],
                ["%s/pathogenic_system/TMHMM" % sample, "", "样品%s的跨膜蛋白分析结果目录" % sample,0,'130036'],
                ["%s/pathogenic_system/secretion_system" % sample, "", "样品%s的分泌系统分析结果目录" % sample,0,'130037'],
                ["%s/pathogenic_system/Resfinder" % sample, "", "样品%s的耐药基因Resfinder预测结果目录" % sample,0,''],
                ["%s/structral_genome" % sample, "", "样品%s的结构基因组分析目录" % sample, 0, '130022'],
                ["%s/circos" % sample, "", "样品%s的基因圈图分析circos目录" % sample,0,"130140"],
                ["%s/cgview" % sample, "", "样品%s的基因圈图分析cgview目录" % sample,0,"130141"],

            ]
            if self.option("analysis") in ["uncomplete"]:
                if not self.is_anno_pip:
                    regexps += [
                        [r"%s/data_QC/.+_Illumina_statistics.xls" % sample, "xls", "二代测序数据质控统计表",0,'130039'],
                        [r"%s/data_QC/.+_PacBio_statistics.xls" % sample, "xls", "三代测序PacBio数据质控统计表",0,'130040'],
                        [r"%s/data_QC/.+_Nanopore_statistics.xls" % sample, "xls", "三代测序Nanopore数据质控统计表", 0, '130040'],
                        [r"%s/genomic_assessment/kmer_frequency/.+_Kmer_frequency.xls" % sample, "xls", "基因组评估kmer频率表",0,'130041'],
                        [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'130042'],
                        [r"%s/assembly_predict/assembly/.+_assembly_details.xls" % sample, "xls", "组装结果详情表",0,'130043'],
                        [r"%s/assembly_predict/assembly/.+_assembly_scaffold_details.xls" % sample, "xls","组装结果scaffolds详情表",0,'130044'],
                        [r"%s/assembly_predict/assembly/.+_assembly_contig_details.xls" % sample, "xls", "组装结果contigs详情表",0,'130045'],
                        [r"%s/assembly_predict/assembly/.+.agp" % sample, "", "基因组组装的scaffold与contig对应关系文件",0,'130046'],
                        [r"%s/assembly_predict/assembly/.+_scaf.fna" % sample, "", "基因组组装的scaffold文件",0,'130047'],
                        [r"%s/assembly_predict/assembly/.+_ctg.fna" % sample, "", "基因组组装的contigs文件",0,'130048'],
                        [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'130042'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.faa" % sample, "", "编码基因预测氨基酸序列",0,'130049'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.fnn" % sample, "", "编码基因预测核苷酸序列",0,'130050'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_CDS.gff" % sample, "", "编码基因预测gff格式统计表",0,'130051'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_CDS_statistics.xls" % sample, "xls", "编码基因预测统计表",0,'130052'],
                        [r"%s/assembly_predict/predict/repeats/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0,'130053'],
                        [r"%s/assembly_predict/predict/repeats/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0,'130054'],
                        [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0,'130055'],
                        [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0,'130056'],
                        [r"%s/assembly_predict/predict/rRNA/.+_rRNA.fnn" % sample, "", "rRNA预测序列文件",0,'130057'],
                        [r"%s/assembly_predict/predict/rRNA/.+_rRNA.gff" % sample, "", "rRNA预测gff统计文件",0,'130058'],
                        [r"%s/assembly_predict/predict/tRNA/.+_tRNA.fnn" % sample, "", "tRNA预测序列文件",0,'130059'],
                        [r"%s/assembly_predict/predict/tRNA/.+_tRNA.gff" % sample, "", "tRNA预测结果文件",0,'130060'],
                        [r"%s/assembly_predict/predict/tRNA/.+_tRNA.struc" % sample, "", "tRNA预测二级结构文件",0,'130061'],
                        [r"%s/assembly_predict/predict/Interpersed_repeat/.+.gff" % sample, "", "散在重复序列预测gff结果文件",0,"130178"],
                        [r"%s/assembly_predict/predict/Interpersed_repeat/.+.out" % sample, "", "散在重复序列预测out结果文件",0,"130179"],
                        [r"%s/assembly_predict/predict/Interpersed_repeat/.+/.+.tbl" % sample, "", "散在重复序列预测tbl结果文件",0,"130180"],
                                                [r"%s/assembly_predict/Plasmid_identification/.+\.detail.xls" % sample, "xls", "质粒鉴定详情表", 0],
                        [r"%s/assembly_predict/Plasmid_identification/.+\.plasmid_anno.xls" % sample, "xls", "质粒注释详情表",
                         0],
                        [r"%s/assembly_predict/Plasmid_identification/all.stat.xls" % sample, "xls", "质粒鉴定中间文件", 1],
                        [r"%s/assembly_predict/predict/sRNA_predict/.+\.sRNA.xls" % sample, "", "sRNA预测详情表", 0],
                    ]

                else:
                     regexps += [
                        [r"%s/gene/.+cds.faa" % sample, "", "编码基因预测氨基酸序列",0,'130049'],
                        [r"%s/gene/.+cds.fnn" % sample, "", "编码基因预测核苷酸序列",0,'130050'],
                    ]

                ##
                ####ori和anno流程共用
                regexps += [
                    [r"%s/project_overview/Analysis_files/.+\.gbk" % sample, "", "基因组的gbk文件",1,'130038'],
                    [r"%s/project_overview/Analysis_files/.+\.ptt" % sample, "", "基因组的ptt文件",1,"130142"],
                    [r"%s/project_overview/Analysis_files/.+\.gff" % sample, "", "基因组的gff文件",1,"130143"],
                    [r"%s/project_overview/All_predict_info/.+\.gbk" % sample, "", "基因组的gbk文件",1,''],
                    [r"%s/project_overview/All_predict_info/.+\.ptt" % sample, "", "基因组的ptt文件",1,""],
                    [r"%s/project_overview/All_predict_info/.+\.gff" % sample, "", "基因组的gff文件",1,""],
                    [r"%s/project_overview/core_faa/.+\.cor_gene.fa" % sample, "", "基因组的核心基因串联文件",0,"130168"], #core gene fa
                    [r"%s/project_overview/core_faa/.+\.cor_gene.xls" % sample, "", "持家基因预测详情表", 0, "130168"],
                    [r"%s/project_overview/core_faa/.+\.cor_gene.fnn" % sample, "", "各个持家基因核苷酸序列", 0, "130168"],
                    [r"%s/project_overview/core_faa/.+\.cor_gene.faa" % sample, "", "各个持家基因氨基酸序列",0,"130168"], #core gene fa
                    [r"%s/annotation/NR/.+_anno_nr.xls" % sample, "xls", "基因NR注释结果详情表",0,'130062'],
                    [r"%s/annotation/Swissprot/.+_anno_swissprot.xls" % sample, "xls", "基因Swiss-Prot注释结果详情表",0,'130063'],
                    [r"%s/annotation/Pfam/.+_anno_pfam.xls" % sample, "xls", "基因Pfam注释结果详情表",0,'130064'],
                    [r"%s/annotation/COG/.+_cog_anno.xls" % sample, "xls", "基因COG注释结果详情表",0,'130065'],
                    [r"%s/annotation/COG/.+_cog_summary.xls" % sample, "xls", "COG注释结果汇总统计表",0,'130066'],
                    [r"%s/annotation/GO/.+_go_anno.xls" % sample, "xls", "基因GO注释结果详情表",0,'130067'],
                    [r"%s/annotation/GO/.+_go_statistics.xls" % sample, "xls", "GO注释统计文件",0,"130145"],
                    [r"%s/annotation/GO/.+_go_list.xls" % sample, "xls", "基因与GO注释ID对应表",0,'130069'],
                    [r"%s/annotation/KEGG/.+_kegg_anno.xls" % sample, "xls", "KEGG注释结果文件",0,'130070'],
                    [r"%s/annotation/KEGG/.+_kegg_level_stat.xls" % sample, "xls", "KEGG注释level统计文件",0,'130071'],
                    [r"%s/annotation/KEGG/kegg_pathway_img.tar.gz" % sample, "", "KEGG注释pathway通路图片的压缩文件",0,'130072'],
                    [r"%s/annotation/Summary/.+_anno_summary.xls" % sample, "xls", "基因组注释结果汇总表",0,'130073'],
                    [r"%s/annotation/Two_component/.+senser_regulator.xls" % sample, "", "样品%s双组分调控系统分析结果文件" % sample,0,"130169"], #新增Two_component

                    [r"%s/annotation/Two_component/.+senser_regulator.stat" % sample, "", "样品%s双组分调控系统分析统计结果文件" % sample,0,"130170"], #新增Two_component

                    [r"%s/structral_genome/promoter_predict/.+_promoter_result.xls" % sample, "xls", "全基因组水平上启动子预测结果表",0,'130074'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_summary.xls" % sample, "xls", "全基因组水平上基因组岛预测结果文件",0,'130075'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_detail.xls" % sample, "xls", "全基因组水平上基因组岛预测结果详情表",0,'130076'],
                    [r"%s/mobile_elements/prephage/.+_prephage_summary.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果统计表",0,'130077'],
                    [r"%s/mobile_elements/prephage/.+_prephage_detail.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果详情表",0,'130078'],
                    [r"%s/mobile_elements/prephage/.+_prephage.fna" % sample, "", "全基因组水平上前噬菌体预测结果序列文件",0,'130079'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_summary.xls" % sample, "xls", "CRISPR_Cas系统预测结果统计表",0,'130080'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_detail.xls" % sample, "xls", "CRISPR_Cas系统预测结果详情表",0,'130081'],
                    [r"%s/mobile_elements/Integron/.+.integron.fna" % sample, "xls", "整合子序列文件",0,"130181"],
                    [r"%s/mobile_elements/Integron/.+.integrons" % sample, "xls", "整合子元素表",0,"130182"],
                    [r"%s/mobile_elements/Integron/.+_stat.xls" % sample, "xls", "整合子分析统计表",0,"130183"],
                    [r"%s/mobile_elements/Integron/._sequence.fna" % sample, "xls", "整合子元素序列文件",0,"130184"],
                    [r"%s/mobile_elements/Integron/.+_Integron_detail.xls" % sample, "xls", "整合子分析详情表",0,"130185"],
                    [r"%s/mobile_elements/Integron/.+.summary" % sample, "xls", "整合子预测结果type类型统计表",0,"130186"],
                    [r"%s/mobile_elements/Is_Predict/.+.enzyme.xls" % sample, "xls", "插入序列预测转座酶详情表",0,"130187"],
                    [r"%s/mobile_elements/Is_Predict/.+.is.xls" % sample, "xls", "插入序列预测is统计表",0,"130188"],
                    [r"%s/mobile_elements/Is_Predict/.+.out" % sample, "xls", "插入序列预测详情表",0,"130189"],
                    [r"%s/mobile_elements/Is_Predict/.+.blast.xls" % sample, "xls", "插入序列预测blast结果表",0,"130190"],
                    [r"%s/mobile_elements/Is_Predict/.+.stat.xls" % sample, "xls", "插入序列预测样品统计表",0,"130191"],

                    [r"%s/metabolic_system/CAZy/.+_cazy_anno.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释结果表",0,'130082'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_family_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释family统计表",0,'130083'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_class_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释class统计表",0,'130084'],
                    [r"%s/metabolic_system/antiSMASH/.+_antismash_anno.xls" % sample, "xls", "次级代谢产物合成基因簇预测结果",0,'130085'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_align.xls" % sample, "xls", "全基因组水平上毒力基因预测表",0,'130086'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_anno.xls" % sample, "xls", "全基因组水平上毒力基因注释表",0,'130087'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_level.xls" % sample, "xls", "全基因组水平上毒力基因分级统计表",0,'130088'],
                    [r"%s/pathogenic_system/CARD/.+_card_align.xls" % sample, "xls", "全基因组水平上耐药基因预测表",0,'130089'],
                    [r"%s/pathogenic_system/CARD/.+_card_anno.xls" % sample, "xls", "全基因组水平上耐药基因注释表",0,'130090'],
                    [r"%s/pathogenic_system/CARD/.+_card_category.xls" % sample, "xls", "全基因组水平上耐药基因分类统计表",0,'130091'],
                    [r"%s/pathogenic_system/PHI/.+_phi_align.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因预测表",0,'130092'],
                    [r"%s/pathogenic_system/PHI/.+_phi_anno.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因注释表",0,'130093'],
                    [r"%s/pathogenic_system/PHI/.+_phi_stat.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析PHI ID统计表",0,'130094'],
                    [r"%s/pathogenic_system/PHI/.+_phi_phenotype.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析表型统计表",0,'130095'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align.xls" % sample, "xls", "转运蛋白预测表",0,'130096'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align_top1.xls" % sample, "xls", "转运蛋白预测表(TOP1)",0,'130097'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_anno.xls" % sample, "xls", "转运蛋白注释表",0,'130098'],
                    [r"%s/pathogenic_system/TCDB/.+_whole_genome_tcdb_anno.xls" % sample, "xls", "全基因组水平上转运蛋白注释表",0,'130099'],
                    [r"%s/pathogenic_system/TMHMM/.+_tmhmm_anno.xls" % sample, "xls", "跨膜蛋白分析结果详情表",0,'130100'],
                    [r"%s/pathogenic_system/TMHMM/.+_whole_genome_tmhmm_anno.xls" % sample, "xls", "全基因组水平上跨膜蛋白注释表",0,'130101'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_genes.xls" % sample, "xls","分泌系统相关基因统计表",0,'130102'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_type.xls" % sample, "xls","分泌系统的类型统计表",0,'130103'],
                    [r"%s/pathogenic_system/Resfinder/.+_disinfinder.class.xls" % sample, "xls","Desinfinder分类统计表",0,''],
                    [r"%s/pathogenic_system/Resfinder/.+_disinfinder.detail.xls" % sample, "xls","Desinfinder分类详情表",0,''],
                    [r"%s/pathogenic_system/Resfinder/.+_resfinder.class.xls" % sample, "xls","Resfinder分类统计表",0,''],
                    [r"%s/pathogenic_system/Resfinder/.+_resfinder.detail.xls" % sample, "xls","Resfinder分类详情表",0,''],
                    [r"%s/pathogenic_system/Resfinder/.+.stat.xls" % sample, "xls","耐药基因Resfinder预测结果统计表",0,''],
                ]
            elif self.option("analysis") in ["complete"]:


                if not self.is_anno_pip:
                    regexps += [
                        [r"%s/data_QC/.+_Illumina_statistics.xls" % sample, "xls", "二代测序数据质控统计表",0,'130039'],
                        [r"%s/data_QC/.+_PacBio_statistics.xls" % sample, "xls", "三代测序数据质控统计表",0,'130040'],
                        [r"%s/data_QC/.+_Nanopore_statistics.xls" % sample, "xls", "三代测序Nanopore数据质控统计表", 0, '130040'],
                        [r"%s/genomic_assessment/kmer_frequency/.+_Kmer_frequency.xls" % sample, "xls", "基因组评估kmer频率表",0,'130041'],
                        [r"%s/assembly_predict/assembly/.+_assembly_summary.xls" % sample, "xls", "基因组组装统计表",0,'130042'],
                        [r"%s/assembly_predict/assembly/.+_assembly_details.xls" % sample, "xls", "组装结果详情表",0,'130043'],
                        [r"%s/assembly_predict/assembly/.+_chromosome[0-9]\.fna" % sample, "", "基因组完成图组装的染色体文件",0,'130106'],
                        [r"%s/assembly_predict/assembly/.+_plasmid[A-Z]\.fna" % sample, "", "基因组完成图组装的质粒文件",0,'130107'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS_statistics.xls" % sample, "xls","全基因组水平编码基因预测统计表",0,'130108'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS_statistics.xls" % sample, "xls","核染色体水平编码基因预测统计表",0,'130109'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS_statistics.xls" % sample, "xls","质粒水平编码基因预测统计表",0,'130110'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS.gff" % sample, "","全基因组水平编码基因预测gff格式统计表",0,'130111'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS.gff" % sample, "","核染色体水平编码基因预测gff格式统计表",0,'130112'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS.gff" % sample, "","质粒水平编码基因预测gff格式统计表",0,'130113'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS.faa" % sample, "","全基因组水平编码基因预测氨基酸序列",0,'130114'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS.faa" % sample, "","核染色体水平编码基因预测氨基酸序列",0,'130115'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS.faa" % sample, "","质粒水平编码基因预测氨基酸序列",0,'130116'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_whole_genome_CDS.fnn" % sample, "","全基因组水平编码基因预测核苷酸序列",0,'130117'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_chromosome[0-9]_CDS.fnn" % sample, "","染色体水平编码基因预测核苷酸序列",0,'130118'],
                        [r"%s/assembly_predict/predict/CDS_predict/.+_plasmid[A-Z]_CDS.fnn" % sample, "","质粒水平编码基因预测核苷酸序列",0,'130119'],
                        [r"%s/assembly_predict/predict/repeats/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0, '130053'],
                        [r"%s/assembly_predict/predict/repeats/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0, '130054'],
                        [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.dat" % sample, "", "串联重复序列预测dat结果文件",0,'130055'],
                        [r"%s/assembly_predict/predict/repeats/.+/.+_TRF.gff" % sample, "", "串联重复序列预测gff结果文件",0,'130056'],
                        [r"%s/assembly_predict/predict/rRNA/.+_rRNA.fnn" % sample, "", "rRNA预测序列文件",0,'130057'],
                        [r"%s/assembly_predict/predict/rRNA/.+_rRNA.gff" % sample, "", "rRNA预测gff统计文件",0,'130058'],
                        [r"%s/assembly_predict/predict/tRNA/.+_tRNA.fnn" % sample, "", "tRNA预测序列文件",0,'130059'],
                        [r"%s/assembly_predict/predict/tRNA/.+_tRNA.gff" % sample, "", "tRNA预测结果文件",0,'130060'],
                        [r"%s/assembly_predict/predict/tRNA/.+_tRNA.struc" % sample, "", "tRNA预测二级结构文件",0,'130061'],
                        [r"%s/methylation/.+motif_detail.xls" % sample, "", "甲基化预测结果文件",0,"130148"],
                        [r"%s/methylation/.+motifs.csv" % sample, "", "甲基化预测结果文件",0,"130149"],
                        [r"%s/project_temp/.+\.raw.assembly.fa" % sample, "", "样本初始组装序列", 0, "130149"],
                        [r"%s/project_temp/.+\.raw.assembly.log" % sample, "", "样本初始组装日志", 0, "130149"],
                        [r"%s/assembly_predict/Plasmid_identification/.+\.detail.xls" % sample, "xls", "质粒鉴定详情表", 0],
                        [r"%s/assembly_predict/Plasmid_identification/.+\.plasmid_anno.xls" % sample, "xls", "质粒注释详情表", 0],
                        [r"%s/assembly_predict/Plasmid_identification/all.stat.xls" % sample, "xls", "质粒鉴定中间文件", 1],
                        [r"%s/assembly_predict/predict/sRNA_predict/.+\.sRNA.xls" % sample, "", "sRNA预测详情表", 0],
                    ]
                ##
                ##ori和anno流程共用
                regexps += [
                    [r"%s/project_overview/Analysis_files/.+_chromosome[0-9]\.gbk" % sample, "", "基因组染色体的gbk文件",1,'130104'],
                    [r"%s/project_overview/Analysis_files/.+_plasmid[A-Z]\.gbk" % sample, "", "基因组质粒的gbk文件",1,'130105'],
                    [r"%s/project_overview/Analysis_files/.+\.ptt" % sample, "", "基因组的ptt文件",1,"130171"], #增加ptt文件
                    [r"%s/project_overview/Analysis_files/.+\.gff" % sample, "", "基因组的gff文件",1,"130172"], #增加gff文件
                    [r"%s/project_overview/core_faa/.+\.cor_gene.fa" % sample, "", "基因组的核心基因串联文件",0,"130168"], #core gene fa
                    [r"%s/project_overview/core_faa/.+\.cor_gene.xls" % sample, "", "持家基因预测详情表", 0, "130168"],
                    [r"%s/project_overview/core_faa/.+\.cor_gene.fnn" % sample, "", "各个持家基因核苷酸序列", 0, "130168"],
                    [r"%s/project_overview/core_faa/.+\.cor_gene.faa" % sample, "", "各个持家基因氨基酸序列",0,"130168"], #core gene
                    [r"%s/project_overview/All_predict_info/.+_chromosome[0-9]\.gbk" % sample, "", "基因组染色体的gbk文件",1,''],
                    [r"%s/project_overview/All_predict_info/.+_plasmid[A-Z]\.gbk" % sample, "", "基因组质粒的gbk文件",1,''],
                    [r"%s/project_overview/All_predict_info/.+\.ptt" % sample, "", "基因组的ptt文件",1,""], #增加ptt文件
                    [r"%s/project_overview/All_predict_info/.+\.gff" % sample, "", "基因组的gff文件",1,""], #增加gff文件

                    [r"%s/annotation/NR/.+_whole_genome_anno_nr.xls" % sample, "xls", "全基因组水平上基因NR注释结果详情表",0,'130120'],
                    [r"%s/annotation/NR/.+_chromosome[0-9]_anno_nr.xls" % sample, "xls", "核染色体水平上基因NR注释结果详情表",0,'130121'],
                    [r"%s/annotation/NR/.+_plasmid[A-Z]_anno_nr.xls" % sample, "xls", "质粒水平上基因NR注释结果详情表",0,'130122'],
                    [r"%s/annotation/Swissprot/.+_chromosome[0-9]_anno_swissprot.xls" % sample, "xls","核染色体水平上基因Swiss-Prot注释结果详情表",0,'130123'],
                    [r"%s/annotation/Swissprot/.+_plasmid[A-Z]_anno_swissprot.xls" % sample, "xls","质粒水平上基因Swiss-Prot注释结果详情表",0,'130124'],
                    [r"%s/annotation/Pfam/.+_whole_genome_anno_pfam.xls" % sample, "xls", "全基因组水平上基因Pfam注释结果详情表",0,'130125'],
                    [r"%s/annotation/Pfam/.+_chromosome[0-9]_anno_pfam.xls" % sample, "xls", "核染色体水平上基因Pfam注释结果详情表",0,'130126'],
                    [r"%s/annotation/Pfam/.+_plasmid[A-Z]_anno_pfam.xls" % sample, "xls", "质粒水平上基因Pfam注释结果详情表",0,'130127'],
                    [r"%s/annotation/COG/.+_cog_anno.xls" % sample, "xls", "基因COG注释结果详情表",0,'130065'],
                    [r"%s/annotation/COG/.+_cog_summary.xls" % sample, "xls", "COG注释结果汇总统计表",0,'130066'],
                    [r"%s/annotation/GO/.+_go_anno.xls" % sample, "xls", "基因GO注释结果详情表",0,'130067'],
                    [r"%s/annotation/GO/.+_go_statistics.xls" % sample, "xls", "GO注释统计文件",0,"130152"],
                    [r"%s/annotation/GO/.+_go_list.xls" % sample, "xls", "基因与GO注释ID对应表",0,'130069'],
                    [r"%s/annotation/KEGG/.+_kegg_anno.xls" % sample, "xls", "KEGG注释结果文件",0,'130070'],
                    [r"%s/annotation/KEGG/.+_kegg_level_stat.xls" % sample, "xls", "KEGG注释level统计文件",0,'130071'],
                    [r"%s/annotation/KEGG/kegg_pathway_img.tar.gz" % sample, "", "KEGG注释pathway通路图片的压缩文件",0,'130072'],
                    [r"%s/annotation/Summary/.+_anno_summary.xls" % sample, "xls", "基因组注释结果汇总表",0,'130073'],
                    [r"%s/annotation/Two_component/.+senser_regulator.xls" % sample, "", "样品%s双组分调控系统分析结果文件" % sample,0,"130173"], #新增Two_component

                    [r"%s/annotation/Two_component/.+senser_regulator.stat" % sample, "", "样品%s双组分调控系统分析统计结果文件" % sample,0,"130174"], #新增Two_component

                    [r"%s/structral_genome/promoter_predict/.+_promoter_result.xls" % sample, "xls", "全基因组水平上启动子预测结果表",0,'130074'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_summary.xls" % sample, "xls", "全基因组水平上基因组岛预测结果文件",0,'130075'],
                    [r"%s/mobile_elements/Genomic_Islands/.+_GI_detail.xls" % sample, "xls", "全基因组水平上基因组岛预测结果详情表",0,'130076'],
                    [r"%s/mobile_elements/prephage/.+_prephage_summary.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果统计表",0,'130077'],
                    [r"%s/mobile_elements/prephage/.+_prephage_detail.xls" % sample, "xls", "全基因组水平上前噬菌体预测结果详情表",0,'130078'],
                    [r"%s/mobile_elements/prephage/.+_prephage.fna" % sample, "", "全基因组水平上前噬菌体预测结果序列文件",0,'130079'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_summary.xls" % sample, "xls", "CRISPR_Cas系统预测结果统计表",0,'130080'],
                    [r"%s/mobile_elements/CRISPR_Cas/.+_CRISPR_Cas_detail.xls" % sample, "xls", "CRISPR_Cas系统预测结果详情表",0,'130081'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_anno.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释结果表",0,'130082'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_family_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释family统计表",0,'130083'],
                    [r"%s/metabolic_system/CAZy/.+_cazy_class_stat.xls" % sample, "xls", "全基因组水平上碳水化合物活性酶注释class统计表",0,'130084'],
                    [r"%s/metabolic_system/antiSMASH/.+_antismash_anno.xls" % sample, "xls", "次级代谢产物合成基因簇预测结果",0,'130085'],
                    [r"%s/metabolic_system/antiSMASH/.+_chromosome[0-9]_antismash_anno.xls" % sample, "xls","核染色体水平上次级代谢产物合成基因簇预测结果",0,'130128'],
                    [r"%s/metabolic_system/antiSMASH/.+_plasmid[A-Z]_antismash_anno.xls" % sample, "xls","质粒水平上次级代谢产物合成基因簇预测结果",0,'130129'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_align.xls" % sample, "xls", "全基因组水平上毒力基因预测表",0,'130086'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_anno.xls" % sample, "xls", "全基因组水平上毒力基因注释表",0,'130087'],
                    [r"%s/pathogenic_system/VFDB/.+_vfdb_level.xls" % sample, "xls", "全基因组水平上毒力基因分级统计表",0,'130088'],
                    [r"%s/pathogenic_system/CARD/.+_card_align.xls" % sample, "xls", "全基因组水平上耐药基因预测表",0,'130089'],
                    [r"%s/pathogenic_system/CARD/.+_card_anno.xls" % sample, "xls", "全基因组水平上耐药基因注释表",0,'130090'],
                    [r"%s/pathogenic_system/CARD/.+_card_category.xls" % sample, "xls", "全基因组水平上耐药基因分类统计表",0,'130091'],
                    [r"%s/pathogenic_system/PHI/.+_phi_align.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因预测表",0,'130092'],
                    [r"%s/pathogenic_system/PHI/.+_phi_anno.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作相关基因注释表",0,'130093'],
                    [r"%s/pathogenic_system/PHI/.+_phi_stat.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析PHI ID统计表",0,'130094'],
                    [r"%s/pathogenic_system/PHI/.+_phi_phenotype.xls" % sample, "xls", "全基因组水平上病原菌与宿主互作分析表型统计表",0,'130095'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align.xls" % sample, "xls", "转运蛋白预测表",0,'130096'],
                    [r"%s/pathogenic_system/TCDB/.+_tcdb_align_top1.xls" % sample, "xls", "转运蛋白预测表(TOP1)",0,'130097'],
                    [r"%s/pathogenic_system/TCDB/.+_whole_genome_tcdb_anno.xls" % sample, "xls", "全基因组水平上转运蛋白注释表",0,'130099'],
                    [r"%s/pathogenic_system/TCDB/.+_chromosome[0-9]_tcdb_anno.xls" % sample, "xls", "核染色体水平上转运蛋白注释表",0,'130130'],
                    [r"%s/pathogenic_system/TCDB/.+_plasmid[A-Z]_tcdb_anno.xls" % sample, "xls", "质粒水平上转运蛋白注释表",0,'130131'],
                    [r"%s/pathogenic_system/TMHMM/.+_whole_genome_tmhmm_anno.xls" % sample, "xls", "全基因组水平上跨膜蛋白注释表",0,'130101'],
                    [r"%s/pathogenic_system/TMHMM/.+_chromosome[0-9]_tmhmm_anno.xls" % sample, "xls", "核染色体水平上跨膜蛋白注释表",0,'130132'],
                    [r"%s/pathogenic_system/TMHMM/.+_plasmid[A-Z]_tmhmm_anno.xls" % sample, "xls", "质粒水平上跨膜蛋白注释表",0,'130133'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_genes.xls" % sample, "xls","分泌系统相关基因统计表",0,'130102'],
                    [r"%s/pathogenic_system/secretion_system/.+_secretion_system_type.xls" % sample, "xls","分泌系统的类型统计表",0,'130103'],
                ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

    #  change_name_map {'pre':"seq",'max_id':1,"map":{"EC_111":"seq1"},'chr_max_id': 1, 'pls_max_id':1}
    def change_atcgn(self,ori_fna,new_fna,change_name_map=None,type=None,from_ncbi=None):  #type : plasmid  chromosome

        with open(ori_fna) as f, open(new_fna,'w') as fw:
            for line in f:
                if line.startswith(">"):
                    line = line.strip()
                    if change_name_map:

                        sline = line.split()
                        seq_name = sline[0][1:]
                        if from_ncbi:
                            type = from_ncbi[seq_name]
                        if seq_name not in change_name_map["map"].keys():
                            if type == 'chromosome':
                                if 'chr_max_id' not in change_name_map.keys():
                                    change_name_map['chr_max_id'] = 1
                                new_name = 'chromosome' + str(change_name_map['chr_max_id'])
                                change_name_map['chr_max_id'] +=1
                            elif type == 'plasmid':
                                if 'pls_max_id' not in change_name_map.keys():
                                    change_name_map['pls_max_id'] = 1
                                new_name = 'plasmid' + str(change_name_map['pls_max_id'])
                                change_name_map['pls_max_id'] +=1
                            else:
                                new_name = change_name_map['pre']+str(change_name_map['max_id'])
                            change_name_map['max_id'] +=1
                            change_name_map["map"][seq_name] = new_name
                        else:
                            raise 'name repeat %s' % seq_name
                        sline[0] = '>'+new_name
                        fw.write(' '.join(sline)+'\n')
                    else:
                        fw.write(line)
                else:
                    line=line.strip()
                    nline = re.sub('[^atcgnATCGN\s]','n',line)
                    fw.write(nline+'\n')

    def change_gff_seq_name(self,ori_gff, change_name_map):
        with open(ori_gff) as fr, open(ori_gff+'_new','w') as fw:
            for line in fr:
                if line.startswith('#'):
                    fw.write(line)
                else:
                    spline = line.split('\t',1)
                    seq_name = spline[0]
                    if seq_name not in change_name_map['map'].keys():
                        raise "gff file %s not in fasta " % seq_name
                    else:
                        spline[0] = change_name_map['map'][seq_name]
                    fw.write('\t'.join(spline))

    def change_name(self,sample):
        sample = sample.replace('_','-')
        return  sample

    ##zouguanqing 201904
    def pre_anno_pipeline(self):
        ########t 测试
        #genome_dir = '/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/bac/bac_update/test_anno/asse_dir'
        #gff_dir = '/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/bac/bac_update/test_anno/gff_dir'

        genome_dir = self.option('fna_dir').prop['path']
        gff_dir = self.option('gff_dir').prop['path']
        self.com_chr_pls_num = {}
        sample_list = []
        fna_info = {}
        fna_type_info = {}
        pat = re.compile('\s+')
        with open(genome_dir + '/list.txt') as fr1:
            fr1.readline()
            for line in fr1:
                spline = pat.split(line)
                sample = spline[0]
                sample = self.change_name(sample)
                if sample not in self.com_chr_pls_num.keys():
                    self.com_chr_pls_num[sample] = {'chr_num':0, 'plasmid_num':0}
                if len(spline) > 2 and self.option('analysis') in ['complete']:
                    if spline[2] == 'chromosome':
                        self.com_chr_pls_num[sample]['chr_num'] +=1
                    elif spline[2] == 'plasmid':
                        self.com_chr_pls_num[sample]['plasmid_num'] +=1
                    fa_file = spline[1]
                    if sample not in fna_info.keys():
                        fna_info[sample] = [fa_file]
                        fna_type_info[sample] = {fa_file:spline[2]}
                    else:
                        fna_info[sample].append(fa_file)
                        fna_type_info[sample][fa_file] = spline[2]
                elif self.option('analysis') in ['uncomplete']:
                    fa_file = spline[1]
                    if sample not in fna_info.keys():
                        fna_info[sample] = [fa_file]
                    else:
                        fna_info[sample].append(fa_file)
                else:
                    raise OptionError("完成图的list文件必须三列！")

        with open(gff_dir + '/list.txt') as fr:
            fr.readline()
            for line in fr:
                spline = pat.split(line)
                sample = spline[0]
                sample = self.change_name(sample)
                if '.' in sample:
                    raise OptionError('样本名不能有"."')
                # if '_' in  sample:
                #     self.logger.error('样本名不能有下划线')
                if sample in sample_list:
                    raise OptionError('一个样本一个gff文件')
                sample_list.append(sample)
        if sorted(fna_info.keys()) != sorted(sample_list):
            raise OptionError('基因组序列的样本名和gff的样本名不一致')
        ##生成self.sample，为了兼容
        self.samples = {}
        for s in sample_list:
            self.samples[s] = s

        with open(gff_dir + '/list.txt') as fr:
            fr.readline()
            for line in fr:
                spline = pat.split(line)
                sample = spline[0]
                sample = self.change_name(sample)
                gff_file = spline[1]
                if sample not in self.samples_seq_map.keys():
                    self.samples_seq_map[sample] = {"pre":"seq","max_id":1,"map":{}}
                #一个样本建一个文件夹
                sample_dir = self.work_dir + '/' + sample
                if  os.path.exists(sample_dir):
                    shutil.rmtree(sample_dir)
                os.mkdir(sample_dir)
                #样本gff文件link到文件夹
                gff_path  = gff_dir + '/'+ gff_file
                new_gff = sample+'.gff'
                os.link(gff_path, sample_dir+'/'+new_gff)  ##修改gff文件名
                if sample not in fna_info.keys():
                    raise OptionError('上传的gff和fna的list.txt的样本名没有对应上')
                seq_dir = sample_dir +'/seq_dir'
                os.mkdir(seq_dir)
                if self.option("analysis") in ["complete"]:
                    for each_fa in fna_info[sample]:
                        each_fa_type = fna_type_info[sample][each_fa]
                        #os.link(genome_dir+'/'+each_fa, seq_dir+'/'+each_fa)
                        self.change_atcgn(genome_dir+'/'+each_fa, seq_dir+'/'+each_fa, change_name_map=self.samples_seq_map[sample],type=each_fa_type)
                        with open(seq_dir+'/'+each_fa) as f_tmp:
                            line1 = f_tmp.readline()
                            seq_name = line1.strip().split()[0][1:]
                        os.rename(seq_dir+'/'+each_fa,seq_dir+'/'+seq_name+'.fasta')

                    fa_path_li = [seq_dir+'/'+i  for i in os.listdir(seq_dir)]
                    os.system('cat {} >  {}/all.fasta'.format(' '.join(fa_path_li),sample_dir))
                else:
                    if len(fna_info[sample])>1:
                        raise OptionError('扫描图一个样本只能有一个基因组序列文件')
                    sample_fa = fna_info[sample][0]
                    #os.link(genome_dir+'/'+sample_fa,sample_dir+'/all.fasta')
                    self.change_atcgn(genome_dir+'/'+sample_fa,sample_dir+'/all.fasta',change_name_map=self.samples_seq_map[sample])
                    fa_dic = {}
                    with open(sample_dir+'/all.fasta') as f:
                        for line in f:
                            if line[0] == '>':
                                line = line.strip()
                                tmp_s = pat.split(line)
                                each = tmp_s[0][1:]
                                if each not in fa_dic.keys():
                                    fa_dic[each] = ''
                                else:
                                    raise OptionError('{} 的{} 序列名有多条'.format(sample_fa, each))
                            else:
                                fa_dic[each] += line
                    for k in fa_dic.keys():
                        with open(sample_dir+'/seq_dir/'+ k+'.fasta','w') as fw:
                            fw.write('>'+k+'\n'+fa_dic[k])

                ##修改gff 序列名称
                self.change_gff_seq_name(sample_dir+'/'+new_gff,change_name_map=self.samples_seq_map[sample])
                os.remove(sample_dir+'/'+new_gff)
                os.rename(sample_dir+'/'+new_gff+"_new",sample_dir+'/'+new_gff)
                if self.option("analysis") in ["complete"]:
                    with open(sample_dir + '/' + new_gff, 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            if re.search("#", line):
                                pass
                            else:
                                lin = line.strip('\n').split('\t')
                                try:
                                    m = re.search("Parent=\s*(.*\D)[0-9]+;", lin[8]) ## 在开头或者是中间的情况 ## fix by qingchen.zhang @20200918
                                except:
                                    m = re.search("Parent=\s*(.*\D)[0-9]+$", lin[8]) ###在结尾## fix by qingchen.zhang @20200918
                                # else:
                                #     m = re.search("ID=([A-Za-z_]*)[0-9]*", lin[8])
                                #     des = m.group(1)
                                if m:
                                    des = m.group(1)
                                    if sample not in self.samples_gene_map:
                                        self.samples_gene_map[sample] = {lin[0]: des}
                                    else:
                                        if lin[0] not in self.samples_gene_map[sample]:
                                            self.samples_gene_map[sample][lin[0]] = des
                else:
                    with open(sample_dir + '/' + new_gff, 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            if re.search("#", line):
                                pass
                            else:
                                lin = line.strip('\n').split('\t')
                                try:
                                    m = re.search("Parent=\s*(.*\D)[0-9]+;", lin[8])## 在开头或者是中间的情况 ## fix by qingchen.zhang @20200918
                                except:
                                    m = re.search("Parent=\s*(.*\D)[0-9]+$", lin[8])###在结尾## fix by qingchen.zhang @20200918
                                # else:
                                #     m = re.search("ID=([A-Za-z_]*)[0-9]*", lin[8])
                                #     des = m.group(1)
                                if m:
                                    des = m.group(1)
                                    if sample not in self.samples_gene_map:
                                        self.samples_gene_map[sample] = des

    def pre_anno_pipline_from_ncbi(self):
        self.com_chr_pls_num = {}
        ##生成self.sample，为了兼容
        self.samples = {}
        self.assemble_id = self.option('assemble_id')
        samples = self.assemble_id.split(',')
        for sample in samples:
            sample_type ={}
            #在工作目录以样本名建文件夹，并放入fna和gff文件，并改名
            new_name = sample.split('.')[0]
            new_name = new_name.replace('_','-')
            if new_name not in self.com_chr_pls_num.keys():
                self.com_chr_pls_num[new_name] = {'chr_num':0, 'plasmid_num':0}
            s_file = self.config.SOFTWARE_DIR + '/database/bacgenome/complete/bacgenome_seq_type.xls'
            with open (s_file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    lin =line.strip().split("\t")
                    if lin[0] == sample:
                        if 'plasmid' in lin[2]:
                            self.com_chr_pls_num[new_name]['plasmid_num'] += 1
                        else:
                            self.com_chr_pls_num[new_name]['chr_num'] += 1
                        sample_type[lin[1]] = lin[2]
            if new_name not in self.samples_seq_map.keys():
                self.samples_seq_map[new_name] = {"pre":"seq","max_id":1,"map":{}}
            else:
                raise  "repeat sample name: %s" %new_name
            self.samples[new_name] = new_name
            s_dir = os.path.join(self.work_dir,new_name)
            s_seqs = s_dir + '/seq_dir'
            if os.path.exists(s_dir):
                shutil.rmtree(s_dir)
            os.mkdir(s_dir)
            os.mkdir(s_seqs)
            s_fna = self.config.SOFTWARE_DIR + '/database/bacgenome/complete/fna/'+sample+'_genomic.fna.gz'
            s_gff = self.config.SOFTWARE_DIR + '/database/bacgenome/complete/gff/'+sample+'_genomic.gff.gz'
            if not os.path.exists(s_fna):
                raise Exception(s_fna + ' not exist')
            if not os.path.exists(s_gff):
                raise Exception(s_gff + ' not exist')

            work_fna =  '{}/{}.fasta'.format(s_dir,'all')
            work_gff =  '{}/{}.gff'.format(s_dir,new_name)
            os.system('zcat {} > {}'.format(s_fna,work_fna+"_tmp"))
            os.system('zcat {} > {}'.format(s_gff, work_gff))
            self.change_atcgn(work_fna+"_tmp", work_fna,change_name_map=self.samples_seq_map[new_name],from_ncbi=sample_type)
            os.remove(work_fna+"_tmp")
            ##

            #拆分fna到seq_dir文件夹下
            tmp_fna_dic = {}
            sp_pat = re.compile('\s*')
            with open(work_fna) as fr:
                for line in fr:
                    if line[0] ==  '>':
                        seq_name = sp_pat.split(line)[0][1:]
                        tmp_fna_dic[seq_name] = line
                    else:
                        tmp_fna_dic[seq_name] += line
            for k in tmp_fna_dic.keys():
                with open(s_seqs+'/'+ k+'.fasta','w') as fw:
                    fw.write(tmp_fna_dic[k])
            #修改gff序列名称
            self.change_gff_seq_name(work_gff,change_name_map=self.samples_seq_map[new_name])
            os.remove(work_gff)
            os.rename(work_gff+"_new",work_gff)
            with open(work_gff, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if re.search("#", line):
                        pass
                    else:
                        lin = line.strip('\n').split('\t')
                        try:
                            m = re.search("Parent=\s*(.*\D)[0-9]+;", lin[8])## 在开头或者是中间的情况 ## fix by qingchen.zhang @20200918
                            des = m.group(1)
                        except:
                            m = re.search("Parent=\s*(.*\D)[0-9]+$", lin[8])###在结尾## fix by qingchen.zhang @20200918
                            des = m.group(1)
                        # else:
                        #     m = re.search("ID=([A-Za-z_]*)[0-9]*", lin[8])
                        #     des = m.group(1)
                        if sample not in self.samples_gene_map:
                            self.samples_gene_map[sample] = {lin[0]: des}
                        else:
                            if lin[0] not in self.samples_gene_map[sample]:
                                self.samples_gene_map[sample][lin[0]] = des

    def this_project_need_change_files(self,sample,pre_path):
        files = [
            ("%s.cds.ffn","spe_fa"), #spe_fa
            ("%s_summary.xls", "csv","Location"),
            #gene
            ("gene/%s.gene.gff","sp_column","Sequence id"), #sp_column
            ("gene/%s.trna.gff","sp_column","Sequence id"), #sp_column
            ("gene/%s.rrna.gff","sp_column","Sequence id"), #sp_column
            ("gene/%s.cds.ffn","spe_fa"),  #spe_fa
            ("gene/%s.cds.faa","spe_fa"),
            ("gene/%s.sequence_len.xls",'csv','Sequence id'),
            #annotation
            ("annotation/CAZy/%s_anno_cazy.xls","csv","Location"),
            ("annotation/COG/%s_cog_anno.xls","csv","Location"),
            ("annotation/GO/%s_go_anno.xls","csv","Location"),
            ("annotation/KEGG/%s_kegg_anno.xls","csv","Location"),
            ("annotation/NR/%s_anno_nr.xls","csv","Location"),
            ("annotation/NR/%s_whole_genome_anno_nr.xls","csv","Location"),
            ("annotation/Pfam/%s_anno_pfam.xls","csv","Location"),
            ("annotation/Pfam/%s_whole_genome_anno_pfam.xls","csv","Location"),
            ("annotation/Summary/%s_anno_summary.xls","csv","Location"),
            ("annotation/Swissprot/%s_anno_swissprot.xls","csv","Location"),
            ("annotation/Swissprot/%s_whole_genome_anno_swissprot.xls","csv","Location"),
            ("annotation/Two_component/%s.senser_regulator.xls","csv","Location"),
            #cgview
            #circos
            #gbk/analysis_gbk/
            #gbk/gbk/
            #gbk/
            ##("gbk/%s_seq1.gff","gff"),
            ##("gbk/%s_seq1.ptt","ptt"),
            # metabolic_system
            ("metabolic_system/antiSMASH/%s_gene_antismash.xls","csv","Location"),
            #("metabolic_system/antiSMASH/%s_antismash_anno.xls","sp_column","Cluster ID"),  #sp_column
            #mobile_elements/
            ("mobile_elements/CRISPR_Cas/%s_CRISPR_Cas_summary.xls","csv","Sample_chr"),
            ("mobile_elements/CRISPR_Cas/%s_CRISPR_Cas_detail.xls","csv","Sample_chr"),
            ("mobile_elements/Genomic_Islands/%s_GI_summary.xls","csv","Location"),
            ("mobile_elements/Genomic_Islands/%s_GI_detail.xls","csv","Location"),
            ("mobile_elements/prephage/%s_prephage_detail.xls","csv","Location"),
            ("mobile_elements/prephage/%s_prephage_summary.xls","csv","Location"),
            #pathogenic_system
            ("pathogenic_system/CARD/%s_card_anno.xls","csv","Location"),
            ("pathogenic_system/PHI/%s_phi_anno.xls","csv","Location"),
            ("pathogenic_system/SIGNALP/%s_Gram-_SignalP.txt","csv","Location"),
            ("pathogenic_system/SIGNALP/%s_Gram+_SignalP.txt","csv","Location"),
            ("pathogenic_system/TCDB/%s_whole_genome_tcdb_anno.xls","csv","Location"),
            ("pathogenic_system/TCDB/%s_tcdb_anno.xls","csv","Location"),
            ("pathogenic_system/TMHMM/%s_whole_genome_tmhmm_anno.xls","csv","Location"),
            ("pathogenic_system/TMHMM/%s_tmhmm_anno.xls","csv","Location"),
            ("pathogenic_system/VFDB/%s_vfdb_anno.xls","csv","Location"),
            ("structral_genome/promoter_predict/%s_promoter_result.xls","csv","Location")
            #tree
            ]

        new_files=[]
        for f in files:
            lf = list(f)
            lf[0]=pre_path + '/'+lf[0]%sample
            new_files.append(lf)

        antismash_anno_files = glob.glob(pre_path+'/metabolic_system/antiSMASH/*_antismash_anno.xls')
        for f in antismash_anno_files:
            new_files.append([f,"sp_column","Cluster ID"])

        return new_files

    #zouguanqing
    def change_back_result_files_name(self):
        if self.samples_seq_map == {}:
            return
        for sample in self.samples_seq_map.keys():
            pre_path = self.output_dir + '/' +sample
            deal_files = self.this_project_need_change_files(sample, pre_path)
            map = self.samples_seq_map[sample]['map']
            self.logger.info("self.samples_seq_map: {}".format(self.samples_seq_map))
            self.logger.info("deal_files: {}".format(deal_files))

            new_map =  {v:k for k, v in map.items()}

            long_map_short = {     #解决有些tool改成缩写导致后面没法将名字改回
                "chromosome" : 'Chr',
                "chromosome1" : "Chr1",
                "chromosome2" : "Chr2",
                "chromosome3" : "Chr3",
                "plasmid" :'p',
                "plasmid1" : 'p1',
                "plasmid2" :'p2',
                "plasmid3" : 'p3',
                "plasmid4" :'p4',
                "plasmid5" : 'p5',
                "plasmid6" :'p6',
                "plasmid7" : 'p7'
            }
            klist = new_map.keys()
            for k in klist:
                if k in long_map_short:
                    new_map[long_map_short[k]] = new_map[k]


            CR = ChangeResultSeqName(new_map)
            CR.change_pip(deal_files)

    def write_map_file(self):
        if self.samples_seq_map == {}:
            return
        with open(self.work_dir+ '/samples_seqs_map.txt','w') as fw:
            for s in self.samples_seq_map.keys():
                for seq in self.samples_seq_map[s]['map'].keys():
                    fw.write('\t'.join([s,seq,self.samples_seq_map[s]['map'][seq]])+'\n')


    def get_lib_type(self,file):
        lib_type = {}
        with open(file,'r') as f:
            lines =f.readlines()
            for line in lines[1:]:
                lin =line.rstrip('\r\n').split('\t')
                if re.search(r'PE',lin[5]):
                    des = lin[5] + str(lin[2])
                    if lin[0] in lib_type.keys():
                        if des not in lib_type.values():
                            lib_type[lin[0]][des] = lin[5]
                    else:
                        lib_type[lin[0]] = {des: lin[5]}
        return lib_type

    def get_assemble_type(self,file):
        dict ={}
        type =[]
        with open(file, "rb") as l:
            raw_lines = l.readlines()
            for line in raw_lines[1:]:
                line2 = line.strip('\r\n').split("\t")
                dict[line2[5]] = line2[5]
        for k in dict.iterkeys():
            type.append(k)
        if len(type) ==1:
            sequence_type=type[0]
        else:
            sequence_type = ','.join(type)
        return sequence_type


    def get_list(self):
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                if os.path.exists(list_path):
                    self.logger.info(list_path)
                with open(list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 6 or len(line) == 7:
                            self.samples[line[0]] = line[0]
                        else:
                            self.set_error('raw_dir的list.txt文件格式有误', code="11400101")
            elif self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                if os.path.exists(list_path):
                    self.logger.info(list_path)
                with open(list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 2 or len(line) == 3:
                            self.samples[line[0]] = line[0]
                        else:
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400102")
            elif self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                raw_sample ={}
                assemble_sample ={}
                with open(raw_list_path, "rb") as f:
                    lines=f.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 6 or len(line) == 7:
                            raw_sample[line[0]] = line[0]
                        else:
                            self.set_error('raw_dir的list.txt文件格式有误', code="11400103")
                with open(assemble_list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 2 or len(line) == 3:
                            assemble_sample[line[0]] = line[0]
                        else:
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400104")
                for key in raw_sample.keys():
                    if key in assemble_sample.keys():
                        self.samples[key]=key
                    else:
                        self.set_error('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致', code="11400105")

        elif self.option("analysis") in ["complete"]:
            if self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                if os.path.exists(list_path):
                    self.logger.info(list_path)
                with open(list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 3 or len(line) == 2:
                            self.samples[line[0]] = line[0]
                        else:
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400106")
            elif self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                raw_sample ={}
                assemble_sample ={}
                with open(raw_list_path, "rb") as f:
                    lines=f.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 6 or len(line) == 7:
                            raw_sample[line[0]] = line[0]
                        else:
                            self.set_error('raw_dir的list.txt文件格式有误', code="11400108")
                        if line[5] not in ['PE', "pe", "MP","mp", "pacbio","Pacbio","Nanopore", "nanopore"]:
                            self.set_error('数据类型不是PE、MP、Pacbio、Nanopore其中一种！')
                        else:
                            if line[5] in ["pacbio","Pacbio"]:
                                if not line[1].endswith(".bam"):
                                    self.set_error('Pacbio原始数据不是bam文件格式！')
                            elif line[5] in ["Nanopore", "nanopore"]:
                                if not line[1].endswith((".fastq", ".gz", ".fq", ".tar.gz")):
                                    self.set_error('Nanopore原始数据不是fastq文件格式！ or 压缩文件不是gz或tar.gz')
                with open(assemble_list_path, "rb") as l:
                    lines=l.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 3 or len(line) == 2:
                            assemble_sample[line[0]] = line[0]
                        else:
                            self.set_error('asse_dir的list.txt文件格式有误', code="11400109")
                for key in raw_sample.keys():
                    if key in assemble_sample.keys():
                        self.samples[key]=key
                    else:
                        self.set_error('raw_dir的list.txt文件和asse_dir的list.txt文件的样品名称不一致', code="11400111")
            elif not self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                raw_sample ={}
                with open(raw_list_path, "rb") as f:
                    lines=f.readlines()
                    for line in lines[1:]:
                        line = line.strip().split('\t')
                        if len(line) == 6 or len(line) == 7:
                            self.samples[line[0]] = line[0]
                        else:
                            self.set_error('raw_dir的list.txt文件格式有误', code="11400108")
                        if line[5] not in ['PE', "pe", "MP", "mp", "pacbio", "Pacbio", "Nanopore", "nanopore"]:
                            self.set_error('数据类型不是PE、MP、Pacbio、Nanopore其中一种！')
                        else:
                            if line[5] in ["pacbio","Pacbio"]:
                                if not re.search("bam", line[1]) :
                                    self.set_error('Pacbio原始数据不是bam文件格式！')
                            elif line[5] in ["Nanopore", "nanopore"]:
                                if not line[1].endswith((".fastq",".gz", ".fq", ".tar.gz")):
                                    self.set_error('Nanopore原始数据不是fastq文件格式！ or 压缩文件不是gz或tar.gz')

    def run_split_dir(self):
        self.logger.info("rr1")
        self.get_list()
        self.logger.info("rr1")
        samples = self.samples
        if self.option("analysis") in ["uncomplete"]:
            if self.option("raw_dir").is_set and not self.option("asse_dir").is_set:
                list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path = os.path.join(self.work_dir,sample)
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    reslut_path = os.path.join(reslut_path, 'data')
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    output = reslut_path + "/list.txt"
                    file = open(output, 'w')
                    with open(list_path, "rb") as l:
                        lines = l.readlines()
                        file.write(lines[0])
                        for line in lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file.write(line )
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files_path = raw.split(',')
                                        for file_path in files_path:
                                            if os.path.exists(reslut_path + "/" + file_path):
                                                os.remove(reslut_path + "/" + file_path)

                                            if not os.path.exists(self.option("raw_dir").prop['path'] + "/" + file_path):
                                                self.set_error('%s 文件没有找到，请检查输入文件'% file_path)
                                            os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                         reslut_path + "/" + file_path)
                                else:
                                    files_path = line2[1].split(',')
                                    for file_path in files_path:
                                        if os.path.exists(reslut_path + "/" + file_path):
                                            os.remove(reslut_path + "/" + file_path)
                                        if not os.path.exists(self.option("raw_dir").prop['path'] + "/" + file_path):
                                                self.set_error('%s 文件没有找到，请检查输入文件'% file_path)
                                        os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                     reslut_path + "/" + file_path)

            if self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    reslut_path = os.path.join(reslut_path, 'assemble')
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    output = reslut_path + "/list.txt"
                    file = open(output, 'w')
                    with open(list_path, "rb") as l:
                        lines = l.readlines()
                        file.write("Sample\tfile\n")
                        list_ass = []
                        for line in lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                if not os.path.exists(self.option("asse_dir").prop['path'] + "/" + line2[1]):
                                    self.set_error('%s 文件没有找到，请检查输入文件'% line2[1])
                                else:
                                    list_ass.append(self.option("asse_dir").prop['path'] + "/" + line2[1])
                        os.system("cat {} >{}".format(" ".join(list_ass), reslut_path +"/"+sample+".fasta"))
                        file.write("{}\t{}\n".format(sample, sample+".fasta"))

            if self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path1 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path1 = os.path.join(reslut_path1, 'data')
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path2 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    reslut_path2 = os.path.join(reslut_path2, 'assemble')
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    output1 = reslut_path1 + "/list.txt"
                    file1 = open(output1, 'w')
                    output2 = reslut_path2 + "/list.txt"
                    file2 = open(output2, 'w')
                    with open(raw_list_path, "rb") as l, open(assemble_list_path, "rb") as f:
                        raw_lines = l.readlines()
                        ass_lines = f.readlines()
                        file1.write(raw_lines[0])
                        file2.write("Sample\tfile\n")
                        for line in raw_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file1.write(line)
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files_path = raw.split(',')
                                        for file_path in files_path:
                                            if os.path.exists(reslut_path1 + "/" + file_path):
                                                os.remove(reslut_path1 + "/" + file_path)
                                            if not os.path.exists(self.option("raw_dir").prop['path'] + "/" + file_path):
                                                self.set_error('%s 文件没有找到，请检查输入文件'% file_path)
                                            os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                         reslut_path1 + "/" + file_path)
                                else:
                                    files_path = line2[1].split(',')
                                    for file_path in files_path:
                                        if os.path.exists(reslut_path1 + "/" + file_path):
                                            os.remove(reslut_path1 + "/" + file_path)
                                        if not os.path.exists(self.option("raw_dir").prop['path'] + "/" + file_path):
                                            self.set_error('%s 文件没有找到，请检查输入文件'% file_path)
                                        os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                     reslut_path1 + "/" + file_path)
                        list_ass =[]
                        for line in ass_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                if not os.path.exists(self.option("asse_dir").prop['path'] + "/" + line2[1]):
                                    self.set_error('%s 文件没有找到，请检查输入文件' % line2[1])
                                else:
                                    list_ass.append(self.option("asse_dir").prop['path'] + "/" + line2[1])
                        os.system("cat {} >{}".format(" ".join(list_ass), reslut_path2 + "/" + sample + ".fasta"))
                        file2.write("{}\t{}\n".format(sample, sample + ".fasta"))

        elif self.option("analysis") in ["complete"]:
            if self.option("asse_dir").is_set and not self.option("raw_dir").is_set:
                list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    reslut_path = os.path.join(reslut_path, 'assemble')
                    if not os.path.exists(reslut_path):
                        os.mkdir(reslut_path)
                    output = reslut_path + "/list.txt"
                    file = open(output, 'w')
                    with open(list_path, "rb") as l:
                        lines = l.readlines()
                        file.write("Sample\tfile\n")
                        list_ass = []
                        for line in lines[1:]:
                            line2 = line.strip('\n\r').split("\t")
                            if line2[0] == sample:
                                if not os.path.exists(self.option("asse_dir").prop['path'] + "/" + line2[1]):
                                    self.set_error('%s 文件没有找到，请检查输入文件' % line2[1])
                                else:
                                    list_ass.append(self.option("asse_dir").prop['path'] + "/" + line2[1])
                        os.system("cat {} >{}".format(" ".join(list_ass), reslut_path + "/" + sample + ".fasta"))
                        file.write("{}\t{}\n".format(sample, sample + ".fasta"))

            if self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                assemble_list_path = os.path.join(self.option("asse_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path1 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path1 = os.path.join(reslut_path1, 'data')
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path2 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    reslut_path2 = os.path.join(reslut_path2, 'assemble')
                    if not os.path.exists(reslut_path2):
                        os.mkdir(reslut_path2)
                    output1 = reslut_path1 + "/list.txt"
                    file1 = open(output1, 'w')
                    output2 = reslut_path2 + "/list.txt"
                    file2 = open(output2, 'w')
                    with open(raw_list_path, "rb") as l, open(assemble_list_path, "rb") as f:
                        raw_lines = l.readlines()
                        ass_lines = f.readlines()
                        file1.write(raw_lines[0])
                        file2.write("Sample\tfile\n")
                        for line in raw_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file1.write(line)
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files_path = raw.split(',')
                                        for file_path in files_path:
                                            if os.path.exists(reslut_path1 + "/" + file_path):
                                                os.remove(reslut_path1 + "/" + file_path)
                                            if not os.path.exists(self.option("raw_dir").prop['path'] + "/" + file_path):
                                                self.set_error('%s 文件没有找到，请检查输入文件'% file_path)
                                            os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                        reslut_path1 + "/" + file_path)
                                else:
                                    files_path = line2[1].split(',')
                                    for file_path in files_path:
                                        if os.path.exists(reslut_path1 + "/" + file_path):
                                            os.remove(reslut_path1 + "/" + file_path)
                                        if not os.path.exists(self.option("raw_dir").prop['path'] + "/" + file_path):
                                            self.set_error('%s 文件没有找到，请检查输入文件'% file_path)

                                        os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                    reslut_path1 + "/" + file_path)
                        list_ass = []
                        for line in ass_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                if not os.path.exists(self.option("asse_dir").prop['path'] + "/" + line2[1]):
                                    self.set_error('%s 文件没有找到，请检查输入文件' % line2[1])
                                else:
                                    list_ass.append(self.option("asse_dir").prop['path'] + "/" + line2[1])
                        os.system("cat {} >{}".format(" ".join(list_ass), reslut_path2 + "/" + sample + ".fasta"))
                        file2.write("{}\t{}\n".format(sample, sample + ".fasta"))
            elif not self.option("asse_dir").is_set and self.option("raw_dir").is_set:
                raw_list_path = os.path.join(self.option("raw_dir").prop['path'], "list.txt")
                for sample in samples.keys():
                    reslut_path1 = os.path.join(self.work_dir, sample)
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    reslut_path1 = os.path.join(reslut_path1, 'data')
                    if not os.path.exists(reslut_path1):
                        os.mkdir(reslut_path1)
                    output1 = reslut_path1 + "/list.txt"
                    file1 = open(output1, 'w')
                    with open(raw_list_path, "rb") as l:
                        raw_lines = l.readlines()
                        file1.write(raw_lines[0])
                        for line in raw_lines[1:]:
                            line2 = line.strip().split("\t")
                            if line2[0] == sample:
                                file1.write(line)
                                if re.search(';', line2[1]):
                                    raw_path = line2[1].split(';')
                                    for raw in raw_path:
                                        files_path = raw.split(',')
                                        for file_path in files_path:
                                            if os.path.exists(reslut_path1 + "/" + file_path):
                                                os.remove(reslut_path1 + "/" + file_path)
                                            if not os.path.exists(
                                                    self.option("raw_dir").prop['path'] + "/" + file_path):
                                                self.set_error('%s 文件没有找到，请检查输入文件' % file_path)
                                            os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                    reslut_path1 + "/" + file_path)
                                else:
                                    files_path = line2[1].split(',')
                                    for file_path in files_path:
                                        if os.path.exists(reslut_path1 + "/" + file_path):
                                            os.remove(reslut_path1 + "/" + file_path)
                                        if not os.path.exists(self.option("raw_dir").prop['path'] + "/" + file_path):
                                            self.set_error('%s 文件没有找到，请检查输入文件' % file_path)
                                        os.link(self.option("raw_dir").prop['path'] + "/" + file_path,
                                                reslut_path1 + "/" + file_path)

    def get_seq_type(self,file, sample):
        dict = {}
        if self.option("analysis") in ["complete"]:
            list = []
            for iterator in SeqIO.parse(file, "fasta"):
                id = iterator.id
                list.append(id)
        dict[sample] = ",".join(list)
        return dict

    def get_seq_dict(self,file):
        dict ={}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                dict[line[0]] = line[-1]
        return dict

    def get_seq_dict3(self,file):
        dict ={}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                dict[line[-1]] = line[-2]
        return dict

    def get_seq_dict2(self,file):
        dict ={}
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                dict[line[-1]] = [line[0],line[-2]]
        return dict

    def move_dir(self, olddir, newdir):  # 原函数名move2outputdir
        """
        移动一个目录下所有文件/文件夹到workflow输出路径下，供set_output调用
        """
        start = time.time()
        if not os.path.isdir(olddir):
            #raise Exception('需要移动到output目录的文件夹不存在。')
            self.set_error('需要移动到output目录的文件夹不存在。', code="11400112")
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}移动到{},耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        """
        递归移动文件夹的内容，供move_dir调用
        """
        if os.path.isfile(old_file):
            if not os.path.isdir(os.path.dirname(new_file)):
                os.makedirs(os.path.dirname(new_file))
            os.link(old_file, new_file)
        elif os.path.isdir(old_file):
            os.makedirs(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)
        else:
            self.logger.info("导出失败：请检查{}".format(old_file))