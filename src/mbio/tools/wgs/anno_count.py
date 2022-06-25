# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# last modify 20180408

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class AnnoCountAgent(Agent):
    """
    snp与indel合并后的vcf进行 注释统计
    """
    def __init__(self, parent):
        super(AnnoCountAgent, self).__init__(parent)
        options = [
            {"name": "snp_anno_genes", "type": "string"},
            {"name": "indel_anno_genes", "type": "string"},
            {"name": "anno_summary", "type": "string"}   # 该文件每个基因组版本对应一个该文件，用于生成pop summary文件
        ]
        self.add_option(options)
        self.step.add_steps('anncount')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.anncount.start()
        self.step.update()

    def step_end(self):
        self.step.anncount.finish()
        self.step.update()
        
    def check_options(self):
        if not self.option("snp_anno_genes"):
            raise OptionError("缺少snp_anno_genes参数", code="34500201")
        if not self.option("indel_anno_genes"):
            raise OptionError("缺少indel_anno_genes参数", code="34500202")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '10G'
        
    def end(self):
        super(AnnoCountAgent, self).end()


class AnnoCountTool(Tool):
    def __init__(self, config):
        super(AnnoCountTool, self).__init__(config)
        self.scrpit_path = self.config.PACKAGE_DIR + "/wgs/anno-count.pl"
        self.perl_path = 'miniconda2/bin/perl '
        
    def cnv_diff(self):
        """
        perl anno-count.pl -snp snp.anno.genes.txt -indel indel.anno.genes.txt
        -anno /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/ann
        -out /mnt/ilustre/users/qingmei.cui/newmdt/Project/BSA_Test_TAIR10/2018.3.2.var/10.annovar/pop
        :return:
        """
        cmd = "{}{} -snp {} -indel {} -anno {} -out {}"\
            .format(self.perl_path, self.scrpit_path, self.option("snp_anno_genes"), self.option("indel_anno_genes"),
                    self.option("anno_summary"), os.path.join(self.output_dir, "pop"))
        self.logger.info(cmd)
        self.logger.info("开始进行anno_count")
        command = self.add_command("anno_count", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("anno_count完成！")
        else:
            self.set_error("anno_count出错！", code="34500201")
            self.set_error("anno_count出错！", code="34500204")

    def run(self):
        super(AnnoCountTool, self).run()
        self.cnv_diff()
        self.end()
