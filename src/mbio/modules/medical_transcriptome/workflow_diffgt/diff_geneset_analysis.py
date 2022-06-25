#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@time    : 2020/10/27
@file    : geneset_analysis.py
"""
import datetime
import glob
import gevent
from biocluster.module import Module
import os
import types
import pandas as pd
import json
from biocluster.core.exceptions import OptionError
import unittest
import shutil
from collections import OrderedDict
from mbio.packages.medical_transcriptome.copy_file import CopyFile

class DiffGenesetAnalysisModule(Module):
    def __init__(self, work_id):
        """
        差异基因集详细分析的模块,对单个基因集进行基因集各模块
        以下几个分析项:
        go注释    go富集    kegg注释  kegg富集
        do注释    do富集    reactome注释  reactome富集
        其中分析模块按一下逻辑进行分布
        物种  层级  go kegg do reactome
        人   T   yes   yes   no  no
        人   G   yes   yes   yes  yes
        非人 T   yes  yes   no   no
        非人 G   yes   yes  no   yes
        """
        super(DiffGenesetAnalysisModule, self).__init__(work_id)
        options = [
            #公共文件
            {'name': 'level', 'type': 'string', 'default': "G"},
            {"name": "annot_result", "type": "string", "default": None},
            {"name": "species", "type": "string", "default": "Homo_sapiens"},
            {'name': 'go_list', 'type': 'infile', 'format': 'medical_transcriptome.go_list'},
            {'name': 'go_version', 'type': 'string', 'default': '20210918'},
            {'name': 'kegg_table', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'kegg_table2', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'all_list', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'add_info', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'reactome_annot', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'do_list', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'kegg_version', 'type': 'string', 'default': None},
            {'name': 'reactome_version', 'type': 'string', 'default': None},
            # 基因集文件
            {'name': 'geneset_name', 'type': 'string', 'default': None},
            {'name': 'genes_num', 'type': 'int', 'default': None},
            {'name': 'geneset_path', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': None},
            {'name': 'regulate', 'type': 'string', 'default': None},
            {'name': 'gene_list', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'gene_multi', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'go_class', 'type': 'infile', 'format': 'medical_transcriptome.common'},
            {'name': 'do_class', 'type': 'infile', 'format': 'medical_transcriptome.common'},
        ]
        self.add_option(options)
        self.go_class_tool = None
        self.kegg_class_tool = self.add_tool("medical_transcriptome.diffgt4work.diff_kegg_class")
        self.go_enrich_tool = self.add_tool("medical_transcriptome.diff_geneset.diff_go_enrich")
        self.kegg_enrich_tool = self.add_module("medical_transcriptome.workflow_diffgt.diff_kegg_enrich")
        self.do_class_tool = None
        self.do_enrich_tool = None
        self.reactome_class_tool = None
        self.reactome_enrich_tool = None
        self.tools_names = ["diff_go_class","diff_kegg_class","diff_go_enrich","diff_kegg_enrich"]
        self.analysis_tools = [self.kegg_class_tool,self.go_enrich_tool,self.kegg_enrich_tool]
        self.analysis_funcs = [self.run_go_class,self.run_kegg_class,self.run_go_enrich,self.run_kegg_enrich]
        self.result_dict = OrderedDict()


    def check_options(self):
        self.logger.info("本次分析参数如下")
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        # if self.option("genes_num") >= 50:
        if self.option("genes_num") > 0:
            if self.option("level") == "G":
                self.logger.info("水平为基因水平,进行reactome分析")
                self.reactome_class_tool = self.add_tool("medical_transcriptome.diff_geneset.diff_reactome_class")
                self.tools_names.append('diff_reactome_class')
                self.analysis_tools.append(self.reactome_class_tool)
                self.analysis_funcs.append(self.run_reactome_class)
                self.reactome_enrich_tool = self.add_module("medical_transcriptome.diff_geneset.diff_reactome_enrich")
                self.tools_names.append('diff_reactome_enrich')
                self.analysis_tools.append(self.reactome_enrich_tool)
                self.analysis_funcs.append(self.run_reactome_enrich)
                if self.option("species") == "Homo_sapiens":
                    self.logger.info("物种为人,且分析水平为基因水平,进行DO分析")
                    self.do_class_tool = self.add_tool("medical_transcriptome.diff_geneset.diff_do_class")
                    self.tools_names.append('diff_do_class')
                    self.analysis_tools.append(self.do_class_tool)
                    self.analysis_funcs.append(self.run_do_class)
                    self.do_enrich_tool = self.add_tool("medical_transcriptome.diff_geneset.diff_do_enrich")
                    self.tools_names.append('diff_do_enrich')
                    self.analysis_tools.append(self.do_enrich_tool)
                    self.analysis_funcs.append(self.run_do_enrich)
                else:
                    self.logger.info("物种不是人,不进行DO分析")
            else:
                self.logger.info("分析水平为转录本水平,不进行reactome分析")
        return True

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def prepare_file_json(self):
        option_names = ["level","species","go_list","kegg_table",
                        "kegg_table2","all_list","add_info","reactome_annot","do_list","kegg_version",
                        "reactome_version","geneset_name","genes_num","geneset_path","level","regulate",
                        "gene_list","gene_multi","go_class","do_class"]
        with open(os.path.join(self.work_dir,"file_list"),"w") as f:
            for option in option_names:
                if option in self._options:
                    try:
                        path = self.option(option).prop["path"]
                        f.write(option+"\t"+path+"\n")
                    except:
                        f.write(option + "\t" + str(self.option(option))+"\n")


    def geneset_analysis(self):
        """

        :return:
        """
        for analysis_tool in self.analysis_funcs:
            gevent.sleep(1)
            analysis_tool()

    def run_go_class(self):
        pass


    def run_go_enrich(self):
        """
           diff_list
               """
        self.logger.info("开始进行go富集分析")

        self.go_enrich_tool.set_options({
            'diff_list': self.option('gene_list').prop["path"],
            'go_list': self.option('go_list').prop["path"],
            'go_version': self.option('go_version'),
        })
        self.go_enrich_tool.run()

    def run_kegg_class(self):
        """
        """
        self.logger.info("开始进行kegg注释分析")
        opts = {
            "geneset_kegg": self.option("gene_multi").prop["path"],
            "kegg_table": self.option("kegg_table").prop["path"],
            "kegg_table2": self.option("kegg_table2").prop["path"],
            "level":self.option("level"),
            "background_links": self.option("add_info").prop["path"],
            "annot_result": self.option("annot_result"),
            "source": "diff_exp",
            "regulate":self.option("regulate"),
            "kegg_version": self.option("kegg_version"),
        }
        self.kegg_class_tool.set_options(opts)
        self.kegg_class_tool.run()

    def run_kegg_enrich(self):
        """
        """
        self.logger.info("开始进行kegg富集分析")
        opts = {
            "geneset_kegg": self.option("gene_multi").prop["path"],
            "kegg_table": self.option("kegg_table").prop["path"],
            "kegg_table2": self.option("kegg_table2").prop["path"],
            "regulate":self.option("regulate"),
            "geneset_list": self.option("gene_list").prop["path"],
            "all_list": self.option("all_list").prop["path"],
            "add_info": self.option("add_info").prop["path"],
            "level":self.option("level"),
            "annot_result": self.option("annot_result"),
            # "task_id": self.option("task_id"),
            "kegg_version": self.option("kegg_version"),
        }
        self.kegg_enrich_tool.set_options(opts)
        self.kegg_enrich_tool.run()

    def run_do_class(self):
        """

        """
        self.logger.info("开始进行do注释分析")
        opts = {
            'do_ids': self.option("do_class").prop["path"],
            'geneset_names': self.option("geneset_name"),
        }
        self.do_class_tool.set_options(opts)
        self.do_class_tool.run()

    def run_do_enrich(self):
        """
        """
        self.logger.info("开始进行do富集分析")
        opts = {
            'diff_list': self.option("gene_list").prop["path"],
            'do_list': self.option("do_list").prop["path"],
        }
        self.do_enrich_tool.set_options(opts)
        self.do_enrich_tool.run()

    def run_reactome_class(self):
        """

        """
        self.logger.info("开始进行reactome注释分析")
        opts = {
            'geneset_ids': self.option('gene_multi').prop["path"],
            'reactome_annot': self.option('reactome_annot').prop["path"],
            'reactome_version': self.option('reactome_version')
        }
        self.reactome_class_tool.set_options(opts)
        self.reactome_class_tool.run()

    def run_reactome_enrich(self):
        """

        """
        self.logger.info("开始进行reactome富集分析")
        opts = {
            'geneset_reactome': self.option('gene_multi').prop["path"],
            'reactome_annot': self.option('reactome_annot').prop["path"],
            'geneset_list': self.option('gene_list').prop["path"],
            'all_list': self.option('all_list').prop["path"],
            # "task_id": self.option("task_id"),
            'level':self.option("level"),
            'reactome_version': self.option('reactome_version')
        }
        self.reactome_enrich_tool.set_options(opts)
        self.reactome_enrich_tool.run()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    try:
                        os.system('rm -r %s' % newfile)
                    except:
                        raise Exception("{}旧文件删除失败".format(newfile))
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                self.logger.info('ln -s %s %s' % (oldfiles[i], newdir))
                try:
                    os.system('ln -s %s %s' % (oldfiles[i], newdir))
                except:
                    raise Exception("{}文件夹移动失败".format(newdir))
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                # try:
                #     os.system('cp -r %s %s' % (oldfiles[i], newdir))
                # except:
                #     raise Exception("{}文件夹移动失败".format(newdir))

    def set_output(self):
        self.logger.info("已完成基因集分析,现在进行文件整理")
        self.result_dict["geneset_name"] = self.option("geneset_name")
        self.result_dict["level"] = self.option("level")
        # if self.option("genes_num") < 50:
        if self.option("genes_num") == 0 :
            self.result_dict["gene_num"] = self.option("genes_num")
            self.result_dict["regulate"] = self.option("regulate")
            self.result_dict["geneset_path"] = self.option("geneset_path")
        else:
            self.result_dict["gene_num"] = self.option("genes_num")
            self.result_dict["regulate"] = self.option("regulate")
            self.result_dict["geneset_path"] = self.option("geneset_path")
            self.result_dict["geneset_list"] = self.option("gene_list").prop["path"]
            self.result_dict["all_list"] = self.option('all_list').prop["path"]
            self.result_dict["tool_names"] =self.tools_names
            self.result_dict["geneset_results"] = OrderedDict()
            if "diff_go_class" in self.tools_names:
                self.logger.info("开始整理go注释相关文件")
                if os.path.exists(os.path.join(self.output_dir,"diff_go_class")):
                    shutil.rmtree(os.path.join(self.output_dir,"diff_go_class"))
                os.makedirs(os.path.join(self.output_dir,"diff_go_class"))
                go_class_path = self.option("go_class").prop["path"]
                try:
                    os.link(go_class_path,os.path.join(self.output_dir,"diff_go_class",os.path.basename(go_class_path)))
                except:
                    raise Exception("diff_go_class文件整理失败")
                self.result_dict["geneset_results"]["diff_go_class"] = os.path.join(self.output_dir,"diff_go_class")
                self.logger.info("go注释相关文件整理完成")
            if "diff_go_enrich" in self.tools_names:
                self.logger.info("开始整理go富集相关文件")
                if os.path.exists(os.path.join(self.output_dir,"diff_go_enrich")):
                    shutil.rmtree(os.path.join(self.output_dir,"diff_go_enrich"))
                os.makedirs(os.path.join(self.output_dir,"diff_go_enrich"))
                go_enrich_result = self.go_enrich_tool.output_dir
                try:
                    # self.linkdir(go_enrich_result,"diff_go_enrich")
                    #add by fwy 20210125
                    CopyFile().linkdir(go_enrich_result, os.path.join(self.output_dir, "diff_go_enrich"))

                except:
                    raise Exception("diff_go_enrich文件整理失败")
                self.result_dict["geneset_results"]["diff_go_enrich"] = os.path.join(self.output_dir, "diff_go_enrich")
                self.logger.info("go富集相关文件整理完成")
            if "diff_kegg_class" in self.tools_names:
                self.logger.info("开始整理kegg注释相关文件")
                if os.path.exists(os.path.join(self.output_dir, "diff_kegg_class")):
                    shutil.rmtree(os.path.join(self.output_dir, "diff_kegg_class"))
                os.makedirs(os.path.join(self.output_dir, "diff_kegg_class"))
                kegg_class_result = self.kegg_class_tool.output_dir
                try:
                    # self.linkdir(kegg_class_result, "diff_kegg_class")
                    # add by fwy 20210125
                    CopyFile().linkdir(kegg_class_result, os.path.join(self.output_dir, "diff_kegg_class"))

                except:
                    raise Exception("diff_kegg_class文件整理失败")
                self.result_dict["geneset_results"]["diff_kegg_class"] = os.path.join(self.output_dir, "diff_kegg_class")
                self.logger.info("kegg注释相关文件整理完成")
            if "diff_kegg_enrich" in self.tools_names:
                self.logger.info("开始整理kegg富集相关文件")
                if os.path.exists(os.path.join(self.output_dir, "diff_kegg_enrich")):
                    shutil.rmtree(os.path.join(self.output_dir, "diff_kegg_enrich"))
                os.makedirs(os.path.join(self.output_dir, "diff_kegg_enrich"))
                kegg_enrich_result = self.kegg_enrich_tool.output_dir
                try:
                    # self.linkdir(kegg_enrich_result, "diff_kegg_enrich")
                    # add by fwy 20210125
                    CopyFile().linkdir(kegg_enrich_result, os.path.join(self.output_dir, "diff_kegg_enrich"))
                except:
                    raise Exception("diff_kegg_enrich文件整理失败")
                self.result_dict["geneset_results"]["diff_kegg_enrich"] = os.path.join(self.output_dir, "diff_kegg_enrich")
                self.logger.info("kegg富集相关文件整理完成")
            if "diff_do_class" in self.tools_names:
                self.logger.info("开始整理do注释相关文件")
                if os.path.exists(os.path.join(self.output_dir, "diff_do_class")):
                    shutil.rmtree(os.path.join(self.output_dir, "diff_do_class"))
                os.makedirs(os.path.join(self.output_dir, "diff_do_class"))
                do_class_result = self.do_class_tool.output_dir
                try:
                    # self.linkdir(do_class_result, "diff_do_class")
                    # add by fwy 20210125
                    CopyFile().linkdir(do_class_result, os.path.join(self.output_dir, "diff_do_class"))
                except:
                    raise Exception("diff_do_class文件整理失败")
                self.result_dict["geneset_results"]["diff_do_class"] = os.path.join(self.output_dir, "diff_do_class")
                self.logger.info("do注释相关文件整理完成")
            if "diff_do_enrich" in self.tools_names:
                self.logger.info("开始整理do富集相关文件")
                if os.path.exists(os.path.join(self.output_dir, "diff_do_enrich")):
                    shutil.rmtree(os.path.join(self.output_dir, "diff_do_enrich"))
                os.makedirs(os.path.join(self.output_dir, "diff_do_enrich"))
                do_enrich_result = self.do_enrich_tool.output_dir
                try:
                    # self.linkdir(do_enrich_result, "diff_do_enrich")
                    # add by fwy 20210125
                    CopyFile().linkdir(do_enrich_result, os.path.join(self.output_dir, "diff_do_enrich"))

                except:
                    raise Exception("diff_do_enrich文件整理失败")
                self.result_dict["geneset_results"]["diff_do_enrich"] = os.path.join(self.output_dir, "diff_do_enrich")
                self.logger.info("do富集相关文件整理完成")
            if "diff_reactome_class" in self.tools_names:
                self.logger.info("开始整理reactome注释相关文件")
                if os.path.exists(os.path.join(self.output_dir, "diff_reactome_class")):
                    shutil.rmtree(os.path.join(self.output_dir, "diff_reactome_class"))
                os.makedirs(os.path.join(self.output_dir, "diff_reactome_class"))
                reactome_class_result = self.reactome_class_tool.output_dir
                try:
                    # self.linkdir(reactome_class_result, "diff_reactome_class")
                    # add by fwy 20210125
                    CopyFile().linkdir(reactome_class_result, os.path.join(self.output_dir, "diff_reactome_class"))
                except:
                    raise Exception("diff_reactome_class文件整理失败")
                self.result_dict["geneset_results"]["diff_reactome_class"] = os.path.join(self.output_dir, "diff_reactome_class")
                self.logger.info("reactome注释相关文件整理完成")
            if "diff_reactome_enrich" in self.tools_names:
                self.logger.info("开始整理reactome富集相关文件")
                if os.path.exists(os.path.join(self.output_dir, "diff_reactome_enrich")):
                    shutil.rmtree(os.path.join(self.output_dir, "diff_reactome_enrich"))
                os.makedirs(os.path.join(self.output_dir, "diff_reactome_enrich"))
                reactome_enrich_result = self.reactome_enrich_tool.output_dir
                try:
                    # self.linkdir(reactome_enrich_result, "diff_reactome_enrich")
                    # add by fwy 20210125
                    CopyFile().linkdir(reactome_enrich_result, os.path.join(self.output_dir, "diff_reactome_enrich"))
                except:
                    raise Exception("diff_reactome_class文件整理失败")
                self.result_dict["geneset_results"]["diff_reactome_enrich"] = os.path.join(self.output_dir, "diff_reactome_enrich")
                self.logger.info("reactome富集相关文件整理完成")
        with open(os.path.join(self.output_dir,"analysis_json"),"w") as f:
            json.dump(self.result_dict, f, indent=2)
        self.end()

    def run(self):
        super(DiffGenesetAnalysisModule, self).run()
        self.prepare_file_json()
        # if self.option("genes_num") >= 50:
        if self.option("genes_num") > 0 :
            rely_tools = self.analysis_tools
            self.on_rely(rely_tools, self.set_output)
            self.geneset_analysis()
        else:
            self.set_output()
        # _ = [self.on_rely(self.basic_tool, tool) for tool in self.predict_tools]



if __name__ == '__main__':
    class TestFunction(unittest.TestCase):
        """
        This is test for the tool. Just run this script to do test.
        """

        def test(self):
            import random
            from mbio.workflows.single import SingleWorkflow
            from biocluster.wsheet import Sheet
            data = {
                "id": "diff_geneset_analysis" + str(random.randint(1, 10000)),
                "type": "module",
                "name": "medical_transcriptome.diff_geneset.diff_geneset_analysis",
                "instant": False,
                "options": dict(
                    task_id='medical_transcriptome',
                    level='G',
                    species="Homo_sapiens",
                    go_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/GO.list",
                    kegg_table="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/gene_kegg_table.xls",
                    kegg_table2="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/gene_kegg_level_table.xls",
                    all_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/all_gene.list",
                    add_info="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/add_info.txt",
                    reactome_annot="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/all_reactome.list",
                    do_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/AnnotPrepare/output/all_do.list",
                    kegg_version="202007",
                    reactome_version="202007",
                    geneset_name="H1_vs_H3_all_12",
                    genes_num=33,
                    geneset_path="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/DiffGenesetPrepare4/output/H1_vs_H3_all_12",
                    regulate="all",
                    gene_list="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/DiffGenesetPrepare4/output/H1_vs_H3_all_12_gene.list",
                    gene_multi = "/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/DiffGenesetPrepare4/output/kegg_class/multi_geneset_list",
                    go_class="/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/DiffGenesetPrepare4/output/go_class/go_class_table.xls",
                    do_class = "/mnt/ilustre/users/sanger-dev/workspace/20200914/Single_geneset_prepare9564/DiffGenesetPrepare/DiffGenesetPrepare4/output/do_class/H1_vs_H3_all_12.do.list"

                )
            }
            wsheet = Sheet(data=data)
            wf = SingleWorkflow(wsheet)
            wf.run()


    unittest.main()
