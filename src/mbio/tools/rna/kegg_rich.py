# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class KeggRichAgent(Agent):
    """
    Kegg富集分析
    version v1.0.1
    author: qiuping
    last_modify: 2016.11.23
    """
    def __init__(self, parent):
        super(KeggRichAgent, self).__init__(parent)
        options = [
            {"name": "kegg_table", "type": "infile", "format": "annotation.kegg.kegg_table"},  # 只含有基因的kegg table结果文件
            {"name": "all_list", "type": "infile", "format": "rna.gene_list"},  # gene名字文件
            {"name": "diff_list", "type": "infile", "format": "rna.gene_list"},
            {"name": "correct", "type": "string", "default": "BH"}  # 多重检验校正方法
        ]
        self.add_option(options)
        self.step.add_steps("kegg_rich")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.kegg_rich.start()
        self.step.update()

    def stepfinish(self):
        self.step.kegg_rich.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('kegg_table').is_set:
            raise OptionError('必须设置kegg的pathway输入文件')
        if self.option('correct') not in ['BY', 'BH', 'None', 'QVALUE']:
            raise OptionError('多重检验校正的方法不在提供的范围内')
        if not self.option("diff_list").is_set:
            raise OptionError("必须设置输入文件diff_list")
        if not self.option("all_list").is_set:
            raise OptionError("必须设置输入文件all_list")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '2G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r"kegg_enrichment.xls$", "xls", "kegg富集分析结果"]
        ])
        super(KeggRichAgent, self).end()


class KeggRichTool(Tool):
    def __init__(self, config):
        super(KeggRichTool, self).__init__(config)
        self._version = "v1.0.1"
        self.kobas = '/bioinfo/annotation/kobas-2.1.1/src/kobas/scripts/'
        self.kobas_path = self.config.SOFTWARE_DIR + '/bioinfo/annotation/kobas-2.1.1/src/'
        self.set_environ(PYTHONPATH=self.kobas_path)
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.python = '/program/Python/bin/'
        self.all_list = self.option('all_list').prop['gene_list']
        self.diff_list = self.option('diff_list').prop['gene_list']

    def run(self):
        """
        运行
        :return:
        """
        super(KeggRichTool, self).run()
        self.run_kegg_rich()

    def run_kegg_rich(self):
        """
        运行kobas软件，进行kegg富集分析
        """
        try:
            self.option('kegg_table').get_kegg_list(self.work_dir, self.all_list, self.diff_list)
            self.logger.info("kegg富集第一步运行完成")
            self.run_identify()
        except Exception as e:
            self.set_error("kegg富集第一步运行出错:{}".format(e))
            self.logger.info("kegg富集第一步运行出错:{}".format(e))

    def run_identify(self):
        kofile = os.path.splitext(os.path.basename(self.option('diff_list').prop['path']))[0]
        cmd_2 = self.python + 'python {}identify.py -f {} -n {} -b {} -o {}.kegg_enrichment.xls'.format(self.config.SOFTWARE_DIR + self.kobas, self.work_dir + '/kofile', self.option('correct'), self.work_dir + '/all_kofile', kofile)
        self.logger.info('开始运行kegg富集第二步：进行kegg富集分析')
        command_2 = self.add_command("cmd_2", cmd_2).run()
        self.wait(command_2)
        if command_2.return_code == 0:
            self.logger.info("kegg富集分析运行完成")
            self.set_output(kofile + '.kegg_enrichment.xls')
            self.end()
        else:
            self.set_error("kegg富集分析运行出错!")

    def set_output(self, linkfile):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        try:
            os.link(linkfile, self.output_dir + '/{}'.format(linkfile))
            self.logger.info("设置kegg富集分析结果目录成功")
        except Exception as e:
            self.logger.info("设置kegg富集分析结果目录失败{}".format(e))
            self.set_error("设置kegg富集分析结果目录失败{}".format(e))
