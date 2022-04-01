# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class MetagenKeggStatAgent(Agent):
    """
    宏基因kegg注释结果统计tool 调用脚本 meta_kegg_stat.py
    author: zhouxuan
    last_modify: 2017.0608
    """

    def __init__(self, parent):
        super(MetagenKeggStatAgent, self).__init__(parent)
        options = [
            {"name": "kegg_result_dir", "type": "infile", "format": "meta_genomic.kegg_dir"},
            {"name": "reads_profile", "type": "infile", "format": "sequence.profile_table"},
            {"name": "kegg_profile_dir", "type": "outfile", "format": "meta_genomic.kegg_dir"},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("kegg_result_dir").is_set:
            raise OptionError("必须设置输入文件夹")
        if not self.option("reads_profile").is_set:
            raise OptionError("必须设置基因丰度表")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
            ])
        super(MetagenKeggStatAgent, self).end()


class MetagenKeggStatTool(Tool):
    def __init__(self, config):
        super(MetagenKeggStatTool, self).__init__(config)
        self._version = "1.0"
        self.python_path = "program/Python/bin/python"
        # self.python_path = self.config.SOFTWARE_DIR + "/program/Python/bin/python"
        self.python_script = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/meta_kegg_stat.py'

    def run(self):
        """
        运行
        :return:
        """
        super(MetagenKeggStatTool, self).run()
        self.run_kegg_stat()
        self.set_output()
        self.end()

    def run_kegg_stat(self):
        kegg_anno = os.path.join(self.option('kegg_result_dir').prop['path'], "kegg_anno.xls")
        enzyme_list = os.path.join(self.option('kegg_result_dir').prop['path'], "kegg_enzyme_list.xls")
        module_list = os.path.join(self.option('kegg_result_dir').prop['path'], "kegg_module_list.xls")
        pathway_list = os.path.join(self.option('kegg_result_dir').prop['path'], "kegg_pathway_list.xls")
        cmd = self.python_path + ' {} -k {} -e {} -p {}  -m {} -r {} -o {} '.\
            format(self.python_script, kegg_anno, enzyme_list, pathway_list, module_list,
                   self.option('reads_profile').prop['path'], self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('kegg_stat', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kegg_stat succeed")
        else:
            self.set_error("kegg_stat failed")
            raise Exception("kegg_stat failed")

    def set_output(self):
        self.logger.info("set_output")
        if len(os.listdir(self.output_dir)) == 5:
            try:
                self.option("kegg_profile_dir", self.output_dir)
            except Exception as e:
                raise Exception("SET_OUTFILE FAILED {}".format(e))
            self.logger.info("OUTPUT RIGHT")
        else:
            raise Exception("OUTPUT WRONG")