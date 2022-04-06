# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
from mbio.packages.metabolome.common import Relation
import subprocess
import re


class TopologyAgent(Agent):
    """
    代谢集富集分析
    last_modify: 2018.6.11
    """

    def __init__(self, parent):
        super(TopologyAgent, self).__init__(parent)
        options = [
            #{"name": "pathway_table", "type": "infile", "format": ""},  # 代谢集文件
            {"name": "compound_table", "type": "string"}, #
            #{"name": "result", "type": "outfile", "format": ""},  # 代谢集总览表
            {"name": "species", "type": "string", "default": "ko"},
            {"name": "method", "type": "string", "default": "rbc"}, # rbc or  rod
            {"name": "backgroup", "type": "string", "default":""},
            {"name": "version", "type": "string", "default": ""}, # 用于区分新老版本
        ]
        self.add_option(options)
        self.step.add_steps("topology")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.topology.start()
        self.step.update()

    def stepfinish(self):
        self.step.topology.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        # if not self.option('pathway_table').is_set:
        #     raise OptionError("必须设置", code="34701001")
        # if not self.option('anno_overview').is_set:
        #     raise OptionError("必须设置代谢总览表", code="34701002")
        # if self.option("correct") not in ["BH", "BY", "bonferroni", "holm"]:
        #     raise OptionError("矫正参数不在范围内，错误参数值：%s", variables=(self.option('correct')), code="34701003")
        return True


    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '8G'

    def end(self):
        super(TopologyAgent, self).end()


class TopologyTool(Tool):
    def __init__(self, config):
        super(TopologyTool, self).__init__(config)
        self.python = '/miniconda2/bin/python'
        #self._path = self.config.PACKAGE_DIR + "/metabolome/impact_value.py"
        if self.option("version") in ['kegg']:
            self.db_path = self.config.SOFTWARE_DIR + "/database/metabolome/topo_json/"
        else:
            self.db_path = self.config.SOFTWARE_DIR + "/database/metabolome/topo_json_{}/".format(self.option("version"))
        self.src = self.config.PACKAGE_DIR + "/metabolome/cal_path_impact.py"
        self.spe = self.option('species')
        self.compound_table = self.option("compound_table")


    def run(self):
        """
        运行
        :return:
        """
        super(TopologyTool, self).run()
        with open(self.compound_table) as fr:
            lines = fr.readlines()
        if len(lines) > 1:
            self.run_impact()
            self.run_plot()
            self.set_output()
            self.end()
        else:
            self.end()


    def run_impact(self):

        #pathway_table = self.option("pathway_table").path
        name2cid = self.compound_table
        out_file = self.work_dir + '/kegg_topology.xls'
        cmd = self.python
        #cmd += ' ' + self.src + ' -p_class {} -cid {} -org {} -method {} -out {} '.format(pathway_table,name2cid,self.option('species'),self.option("method"), out_file)
        cmd += ' ' + self.src + ' -cid {} -org {} -method {} -out {} -db_json {}'.format(name2cid,self.spe,self.option("method"), out_file, self.db_path)
        if self.option("backgroup"):
            cmd = "{} -backgroup {}".format(cmd, self.option("backgroup"))
        # command = self.add_command('cal_impact',cmd)
        # self.logger.info(command)
        # command.run()
        # self.wait()
        # if command.return_code == 0:
        #     self.logger.info('运行成功： %s'%cmd)
        # else:
        #     self.logger.set_error('运行失败： %s'%cmd)
        try:
            cmd =self.config.SOFTWARE_DIR + "/"+ cmd
            self.logger.info('开始运行： %s' %cmd)
            subprocess.check_output(cmd,shell=True)
        except Exception as e:
            self.logger.info('运行失败： %s'%cmd)
            self.logger.info(e)


    def run_plot(self):
        self.plot_src =  self.config.PACKAGE_DIR + "/metabolome/draw_cpd_dag_as_pathways.py"
        db_dir = self.db_path
        spe = self.spe
        pic_out = self.work_dir + '/pic'
        try:
            cmd = self.config.SOFTWARE_DIR + "/%s %s %s %s %s %s"%(self.python, self.plot_src, self.compound_table, pic_out, db_dir,spe)
            if self.option("backgroup"):
                cmd += ' ' +  self.option("backgroup")
            self.logger.info('开始运行： %s' %cmd)
            subprocess.check_output(cmd,shell=True)
        except:
            self.logger.info('运行失败： %s'%cmd)


        # cmd = "%s %s %s %s %s %s"%(self.python, self.plot_src, self.compound_table, pic_out, db_dir,spe)
        # command = self.add_command('plot_fun',cmd)
        # self.logger.info(command)
        # command.run()
        # self.wait()
        # if command.return_code == 0:
        #     self.logger.info('运行成功： %s'%cmd)
        # else:
        #     self.logger.set_error('运行失败： %s'%cmd)
    def change_ko_to_map(self):
        kegg_xls = pd.read_table(self.work_dir +'/kegg_topology.xls',sep='\t',header=0)
        new_ko = kegg_xls['ko'].apply(lambda x: re.sub('^ko','map',x))
        kegg_xls['ko'] = new_ko
        kegg_xls.to_csv(self.work_dir +'/kegg_topology.xls_new',sep='\t',index=False)

    def set_output(self):
        all_files = ['kegg_topology.xls']
        for each in all_files:
            if not os.path.exists(self.work_dir + '/'+each):
                continue
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if each == 'kegg_topology.xls':
                if self.spe == 'ko':
                    self.change_ko_to_map()
                    each = 'kegg_topology.xls_new'
            if os.path.exists(link):
                os.remove(link)
            os.link(self.work_dir + '/'+each, link)

        if os.path.exists(self.work_dir+'/pic'):
            pics = os.listdir(self.work_dir+'/pic')
            for f in pics:
                if f.endswith('svg'):
                    continue
                ori_f = self.work_dir+'/pic/' + f
                tar_f = self.output_dir + '/' + f
                if self.spe == 'ko':
                    nf = re.sub('^ko','map',f)
                    tar_f = self.output_dir + '/' + nf
                if os.path.exists(tar_f):
                    os.remove(tar_f)
                os.link(ori_f,tar_f)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "topo",
            "type": "tool",
            "name": "metabolome.metabset.topology",
            "instant": True,
            "options": dict(
                anno_overview="/mnt/ilustre/users/sanger-dev/workspace/20190416/MetabsetEnrich_tsg_33827_84988_847666/anno_overview_input.xls",
                metabset="/mnt/ilustre/users/sanger-dev/workspace/20190416/MetabsetEnrich_tsg_33827_84988_847666/metabset_input.set.xls",
                correct="BH",
                #bg="project",
                bg="species",
                species="All"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
