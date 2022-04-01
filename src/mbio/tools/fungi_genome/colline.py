
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import shutil
import xml.etree.ElementTree as ET
from biocluster.config import Config
from biocluster.core.exceptions import OptionError


class CollineAgent(Agent):
    """
    author: guanqing.zou
    last_modify: 20180531
    """
    def __init__(self, parent):
        super(CollineAgent, self).__init__(parent)
        options = [
            {"name": "qfna","type": "string","default": ""},
            {"name": "rfna","type": "string", "default": ""},
            {"name":"sample","type": "string", "default": "sample"}

            # # {"name": "bsnout", "type": "infile", },
            # {"name": "bsnou", "type": "string", "default": ""},
            # # {"name": "database", "type": "infile", }
            # #{"name": "database", "type": "string", "default": ""}
            ]
        self.add_option(options)
        self.step.add_steps('colline')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        #self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def step_start(self):
        self.step.colline.start()
        self.step.update()

    def step_end(self):
        self.step.colline.finish()
        self.step.update()

    def check_options(self):
        #if not self.option("bsnout").is_set:
        if not self.option("qfna"):
            raise OptionError("必须设置参数qfna", code="32100601")
        if not self.option("rfna"):
            raise OptionError("必须设置参数rfna", code="32100602")
        return True

    def set_resource(self):
        self._cpu = 6
        self._memory = '30G'

    def end(self):
        super(CollineAgent, self).end()


class CollineTool(Tool):
    def __init__(self, config):
        super(CollineTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/Sibelia-3.0.6-Linux/bin")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/Sibelia-3.0.6-Linux/lib")
        self.sibelia_path = "bioinfo/Genomic/Sofware/Sibelia-3.0.6-Linux/bin/Sibelia"
        self.python_path = "program/Python/bin/python"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.circos_path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/circos-0.69-6/bin/circos"
        self.conf_path = os.path.join(self.config.SOFTWARE_DIR, "database/fungi_circos_conf/circos.conf")


    def run_sibelia(self):
        allfna = os.path.join(self.work_dir, 'All.fna')
        self.qfna_new = os.path.join(self.work_dir, 'q.new.fna')
        self.rfna_new = os.path.join(self.work_dir, 'r.new.fna')
        os.system("sed 's/>/>q_/' {} > {}".format(self.option('qfna'), self.qfna_new ))
        os.system("sed 's/>/>r_/' {} > {}".format(self.option('rfna'), self.rfna_new ))
        os.system("cat {} {} > {}".format(self.qfna_new, self.rfna_new, allfna))
        cmd_str = "{} -s loose {}".format(self.sibelia_path, allfna)

        # outputfile = os.path.join(self.output_dir, "result.out.xls")

        # cmd_path = os.path.join(self.config.PACKAGE_DIR, "annotation/dfvf_merge_info.py")

        self.logger.info("开始运行{}".format(cmd_str))
        command = self.add_command("sibelia",cmd_str)
        command.run()
        self.wait(command)
        if command.return_code == 0:
                self.logger.info("{} 运行完成".format(cmd_str))
                #self.logger.info(outputfile)
        else:
                self.set_error("%s 运行失败", variables=(cmd_str), code="32100601")


    def run_newpos_circos(self):
        fr = open(self.rfna_new)
        rfnaFirst = fr.readline().strip().split(' ')[0][1:]
        fr.close()
        block_coords = os.path.join(self.config.PACKAGE_DIR, "fungi_genome/colline_blocks_coords.py")
        block_file = os.path.join(self.work_dir,"blocks_coords.txt")
        cmd1 = "{} {} {} {}".format(self.python_path, block_coords, block_file, rfnaFirst)
        command1 =self.add_command("block_deal",cmd1)
        self.logger.info("开始运行{}".format(cmd1))
        command1.run()
        self.wait(command1)
        pos_dic =  os.path.join(self.work_dir,"pos.dic")
        highlight_file = os.path.join(self.work_dir, "circos/circos.highlight.txt")
        segdup_file = os.path.join(self.work_dir, "circos/circos.segdup.txt")
        change_pos = os.path.join(self.config.PACKAGE_DIR, "fungi_genome/colline_changePos.py")
        conf_new_path = os.path.join(self.work_dir, "circos/circos.conf")
        if  command1.return_code == 0:
            # os.system("cd circos")
            cmd2 = "{} {} {} {} {}".format(self.python_path, change_pos, pos_dic, segdup_file, highlight_file)
            command2 = self.add_command("change_pos", cmd2)
            command2.run()
            self.wait(command2)

        os.system("cp {} {}".format(self.conf_path, conf_new_path))
        cmd3 = "{} {} -conf {}".format(self.perl_path, self.circos_path, conf_new_path)
        command3 = self.add_command("circos",cmd3)
        self.logger.info("开始运行circos")
        command3.run()
        self.wait(command3)
        if command3.return_code == 0:
            self.logger.info("{} 运行完成".format(cmd3))
            #self.logger.info(outputfile)
        else:
            self.set_error("%s 运行失败", variables=(cmd3), code="32100602")

    def set_output(self):
        olds = "{0}/blocks_coords.txt,{0}/coverage_report.txt,{0}/circos.png,{0}/circos.svg,{0}/block_new.xls".format(self.work_dir).split(',')
        # news = "{0}/{1}_blocks_coords.txt,{0}/{1}_coverage_report.txt,{0}/circos.png,{0}/{1}_circos.svg,{0}/{1}_block_new.xls".format(self.output_dir,self.option('sample')).split(',')
        # circos 的png/svg图片名称前缀不同，导致svg图片下载错误 ghd @20190423
        news = "{0}/{1}_blocks_coords.txt,{0}/{1}_coverage_report.txt,{0}/circos.png,{0}/circos.svg,{0}/{1}_block_new.xls".format(self.output_dir,self.option('sample')).split(',')
        for f in news:
            if os.path.exists(f):
                os.remove(f)
            os.link(olds[news.index(f)],f)

    def run(self):
        """
        运行
        :return:
        """
        super(CollineTool, self).run()
        self.run_sibelia()
        self.run_newpos_circos()
        self.set_output()
        self.end()


