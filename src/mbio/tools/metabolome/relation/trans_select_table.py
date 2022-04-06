# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last modify date: 2021.11.24
# last modified: zhaoyuzhuo

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import numpy as np
import commands
from mbio.packages.metabolome.get_data import dump_trans_data
from biocluster.config import Config


class TransSelectTableAgent(Agent):
    """
    代谢转录关联分析--筛选转录丰度表
    """
    def __init__(self, parent):
        super(TransSelectTableAgent, self).__init__(parent)
        options = [
            {"name": "trans_exp_main_id", "type": "string"},
            {"name": "trans_geneset_main_id", "type": "string", 'defult': ""},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},  # 样本分组名
            {"name": "group_method", "type": "int", "default": 0},
            {"name": "top", "type": "int", "default": 0},
            {"name": "select_table", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("group").is_set:
            raise OptionError("请设置group分组", code="34002502")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(TransSelectTableAgent, self).end()


class TransSelectTableTool(Tool):
    def __init__(self, config):
        super(TransSelectTableTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/metabolome/scripts/profile_select.py'

    def run(self):
        super(TransSelectTableTool, self).run()
        self.run_geneset()
        self.set_output()
        self.end()

    def run_geneset(self):
        # 从关联表中取数据
        metab_client = Config().get_mongo_client(mtype="metabolome")
        relation_db = metab_client[Config().get_mongo_dbname("metabolome")]
        relation_info = relation_db['sg_relation_analysis'].find_one({"task_id": self.option('task_id'), "delete_by": ""})
        relate_task_id = relation_info["relate_task_id"]
        relate_project_type = relation_info["relate_project_type"]
        try:
            db_version = relation_info["relate_db_version"]
        except:
            db_version = 1
        matab_sample_list = relation_info['sp_name']
        trans_sample_list = relation_info['relate_sp_name']
        sample_dict = dict(zip(trans_sample_list, matab_sample_list))
        self.logger.info("matab_sample_list为{}".format(matab_sample_list))
        self.logger.info("trans_sample_list为{}".format(trans_sample_list))
        # 获取基因集表
        self.logger.info("geneset_id为{}".format(self.option("trans_geneset_main_id")))
        if self.option("trans_geneset_main_id") == None:
            pass
        else:
            select_genes_file = self.work_dir + "/select_genes.xls"
            dump_trans_data(outfile=select_genes_file, proj_type=relate_project_type, task_id=relate_task_id,
                            col_type="geneset", db_version=db_version, main_id=self.option("trans_geneset_main_id"))
        # 获取基因name表
        genename_table = self.work_dir + "/gene_id2name.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="gene_name",
                        db_version=db_version, outfile=genename_table)
        #获取转录表达量表
        origin_table = self.work_dir + "/trans_exp_table.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="exp",
                        db_version=db_version, main_id=self.option("trans_exp_main_id"), outfile=origin_table)
        # 根据转录代谢样本名对应关系，修改转录表的样本名称
        df_origin = pd.read_table(origin_table, '\t')
        df_trans_exp = df_origin.ix[:,trans_sample_list]
        df_final = df_trans_exp.rename(sample_dict, axis="columns")
        df_final.insert(0, "gene_id", df_origin["gene_id"])
        trans_exp_table = self.output_dir + "/trans_table.xls"
        df_final.to_csv(trans_exp_table, '\t', header=True, index=False)
        # 筛选转录表达量数据表
        outfile = self.output_dir + "/trans_select_table.xls"
        cmd = self.python_path + ' {} -i {} -o {} -sc gene_id -top {} -merge nomerge'.format(self.script, trans_exp_table, outfile, str(self.option("top")))
        if not self.option("trans_geneset_main_id") == None:
            cmd += ' -s ' + select_genes_file
        if self.option("group").is_set:
            cmd += " -g " + self.option("group").prop["path"] + " -gm " + str(self.option("group_method"))
        self.logger.info(cmd)
        command = self.add_command('trans_table_select', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("table_select succeed")
        elif command.return_code in [-9, 1]:  # modified return_code by guhaidong @ 20180628
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("table_select failed", code="34002501")

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        outfile = self.output_dir + "/trans_select_table.xls"
        cmd = '/usr/bin/head -n 10 {} | /usr/bin/wc -l'.format(outfile)
        _, lines = commands.getstatusoutput(cmd)
        if int(lines) > 1:
            self.option("select_table").set_path(outfile)
            self.logger.info("设置输出结果文件成功")
            self.logger.info(self.option("select_table").path)
        else:
            self.set_error("原始表格在该基因集筛选下为空！", code="34002502")
