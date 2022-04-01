# -*- coding: utf-8 -*-
# __author__ = 'zhujuan'
# __last_modify_: 20180419

from biocluster.workflow import Workflow
import os
import re
import glob
from bson import ObjectId
from mbio.packages.beta_diversity.filter_newick import get_level_newicktree
from bson.objectid import ObjectId
from mbio.packages.meta.common_function import get_save_pdf_status,get_name
from mbio.packages.meta.save_params import save_params


class MicropitaWorkflow(Workflow):
    """
    MicroPITA（microbiomes: Picking Interesting Taxonomic Abundance） 是一种在分级研究中挑选样本的计算工具， 能够更有效的分配资源， 降低研究成本， 最大化的利用样本。
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(MicropitaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_file", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的OTU id
            {"name": "level", "type": "int"},
            {"name": "otu_id", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "group_file", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "submit_location", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "target_feature", "type": "string", "default": ""},  # 筛选某些特定物种，物种间以“；”连接
            {"name": "target_metrics", "type": "string", "default": "rank"},  # target_feature的度量方法，[rank/abundance]
            {"name": "distance_method", "type": "string"},  # 距离算法
            {"name": "diversity_index", "type": "string"},  # 多样性指数
            {"name": "filter_nu", "type": "int"},  # 挑选的样品量
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.micropita = self.add_tool('meta.micropita')

    def change_otuname(self, tablepath, file_name):
        newtable = self.work_dir + "/" + file_name + "_input_abund.xls"
        with open(tablepath, "r") as f, open(newtable, "w") as g:
            head = f.readline()
            g.write(head)
            for line in f:
                lines = line.split("\t", 1)
                specimen = re.subn("^.*; ", "", lines[0])[0]
                g.write(specimen + "\t" + lines[1])
        return  newtable

    def run_micropita(self):
        newtable = self.change_otuname(self.option('otu_file').prop['path'],"otutable")
        options = {
            'otutable': newtable,
            'distance_method': self.option('distance_method'),
            'diversity_index': self.option('diversity_index'),
            'filter_nu': self.option('filter_nu')
        }
        if self.option('group_file').is_set:
            options['group_file'] = self.option('group_file')
        if self.option('target_feature') != "":
            options['target_feature'] = self.option('target_feature')
        self.micropita.set_options(options)
        self.micropita.on('end', self.set_db)
        self.micropita.run()

    def run(self):
        self.run_micropita()
        super(MicropitaWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_micropita = self.api.micropita
        micropita_id = ObjectId(self.option("main_id"))
        pcoa_file = self.micropita.output_dir + "/pcoa.xls"
        select_file = self.micropita.output_dir + "/microPITA_select.xls"
        api_micropita.add_micropita_detail(micropita_id, pcoa_file, select_file)
        #self.end()
        self.save_pdf()

    def save_pdf(self):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "sg_micropita")
            self.figsave = self.add_tool("meta.fig_save")
            self.figsave.on('end', self.end)
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "meta",
                "submit_loc": "micropita",
                "interaction": 1,
                "main_table": "sg_micropita",
            })
            self.figsave.run()
        else:
            self.end()

    def end(self):
        self.output_dir = self.micropita.output_dir
        if self.pdf_status:
            pdf_outs = self.figsave.output_dir
            os.system("cp -r {}/* {}/".format(pdf_outs, self.output_dir))
        save_params(self.output_dir, self.id)
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Micropita分析结果目录", 0, "110234"],
            ["./MicroPITA_select.xls", "xls", "MicroPIT分析结果表", 0, "110235"],
            ["./MicroPITA分析图.pdf", "pdf", "MicroPITA分析图", 0, ""],
        ])
        super(MicropitaWorkflow, self).end()