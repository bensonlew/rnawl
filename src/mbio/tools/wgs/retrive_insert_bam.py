# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# lasted modified by hongdong@20180603

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class RetriveInsertBamAgent(Agent):
    """
    软件:perl
    """
    def __init__(self, parent):
        super(RetriveInsertBamAgent, self).__init__(parent)
        options = [
            {"name": "bam", "type": "string"},
            {"name": "id", "type": "string"},
            {"name": "fai", "type": "string"}  # pop.fa.fai
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam"):
            raise OptionError("必须输入bam文件", code="34504601")
        if not self.option("id"):
            raise OptionError("必须输入id", code="34504602")
        if not self.option("fai"):
            raise OptionError("必须输入fai文件", code="34504603")

    def set_resource(self):
        self._cpu = 4
        self._memory = "50G"

    def end(self):
        super(RetriveInsertBamAgent, self).end()


class RetriveInsertBamTool(Tool):
    def __init__(self, config):
        super(RetriveInsertBamTool, self).__init__(config)
        self.perl = 'program/perl-5.24.0/bin/perl5.24.0'
        self.script = self.config.PACKAGE_DIR + "/wgs/retrive.insert.bam.pl"

    def run_retrive_insert_bam(self):
        """
        perl /mnt/ilustre/users/nanshan.yang/newmdt/Pipline/09.insert.v3.0/bin/bin/retrive.insert.bam.pl -id 30211907
        -fai /mnt/ilustre/users/nanshan.yang/newmdt/Pipline/09.insert.v2.0/demo/step01.index/pop.fa.fai
        -bam /mnt/ilustre/users/nanshan.yang/newmdt/Pipline/09.insert/demo/step02.bwa-mapping/30211907.b1.bam
        -out /mnt/ilustre/users/nanshan.yang/newmdt/Pipline/09.insert.v3.0/demo/step03.bam-insert/
        """
        bam_path = os.path.join(self.work_dir, "bams")
        if not os.path.exists(bam_path):
            os.mkdir(bam_path)
        else:
            os.system('rm -r {}'.format(bam_path))
            os.mkdir(bam_path)
        cmd = "{} {} -id {} -fai {} -bam {} -out {}"\
            .format(self.perl, self.script, self.option("id"), self.option("fai"), self.option("bam"), bam_path)
        self.logger.info(cmd)
        command = self.add_command("retrive_insert_bam", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("Retrive.insert.bam运行完成")
        else:
            self.set_error("Retrive.insert.bam运行失败", code="34504601")

    def set_output(self):
        """
        要将产品线的fq名字30211907.sca5-2283205-2284356.fq， 改成JY102.chr1:10000-20000.fastq.gz
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        path = os.path.join(self.work_dir, "bams/")
        results = os.listdir(path)
        if len(results) == 2:  # combine文件夹与qc.xls
            self.set_error("提取出来的fasq文件列表为空！", code="34504602")
        self.gzip_file(path)
        for f in os.listdir(path):
            if f == "combine":
                pass
            elif f == "stat.xls":
                os.link(self.work_dir + "/bams/stat.xls", self.output_dir + "/stat.xls")
            else:
                new_name = '.'.join([f.split(".")[0], ':'.join([f.split(".")[1].split('-')[0],
                                                                '-'.join(f.split(".")[1].split('-')[1:3])]), 'fastq.gz'])
                os.link(path + f, self.output_dir + "/" + new_name)

    def gzip_file(self, dir_path):
        for m in os.listdir(dir_path):
            if m not in ['combine', "stat.xls"]:
                if not re.match(r".*\.gz$", m):
                    code = os.system("gzip {}".format(os.path.join(dir_path, m)))
                    if code != 0:
                        self.set_error("压缩文件失败！", code="34504603")

    def run(self):
        super(RetriveInsertBamTool, self).run()
        self.run_retrive_insert_bam()
        self.set_output()
        self.end()


