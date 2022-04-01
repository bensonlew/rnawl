# -*_ coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modified: 20171218
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os


class Bcl2fastqAgent(Agent):
    """
    bcl2fastq2-v2.17.1.14
    文库拆分
    """
    def __init__(self, parent=None):
        super(Bcl2fastqAgent, self).__init__(parent)
        options = [
            {"name": "data_path", "type": "string"},  # 下机数据路径
            {"name": "sample_sheet", "type": "infile", "format": "datasplit.sample_sheet"},  # csv文件
            {"name": "bases_mask", "type": "string"},  # 测序模式，hiseq4000 为y151,i6nn,y151; miseq为y301,i6,y301
            {"name": "barcode_mismatch", "type": "int", "default": 0},  # --barcode-mismatches，barcode错配数
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        """参数检查"""
        if not self.option("data_path"):
            raise OptionError("缺少下机数据路径")
        if not self.option("sample_sheet").is_set:
            raise OptionError("缺少csv表")
        if not self.option("bases_mask"):
            raise OptionError("缺少测序模式")
        if self.option("barcode_mismatch") < 0:
            raise OptionError("参数barcode_mismatch:{}不能小于0，请检查".format(self.option("barcode_mismatch")))
        if re.match(r"\d+", self.option("bases_mask")):
            raise OptionError("测序模式有问题，请检查")

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 15
        self._memory = "60G"

    def end(self):
        super(Bcl2fastqAgent, self).end()


class Bcl2fastqTool(Tool):
    def __init__(self, config):
        super(Bcl2fastqTool, self).__init__(config)
        self._version = 2.0
        self.bcl2fastq_path = "bioinfo/seq/bcl2fastq2-v2.20/bin/bcl2fastq"
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/lib64')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.4.0/bin')

    def new_split_path(self, split_path):
        if os.path.exists(split_path):
            return split_path
        if "ilustre" in split_path:
            split_path1 = split_path.replace("ilustre", "clustre")
            if os.path.exists(split_path1):
                return split_path1
        if "sglustre" in split_path:
            split_path1 = split_path.replace("sglustre", "ilustre")
            if os.path.exists(split_path1):
                return split_path1
        return split_path

    def bcl2fastq(self):
        """
        运行bcl2fastq
        """
        self.data_path = self.option("data_path")
        # if not os.path.exists(self.data_path):
        #     self.data_path = self.option("data_path").replace("ilustre", "clustre")
        #     if not os.path.exists(self.data_path):
        #         self.set_error("下机数据路径:%s 没有找到，请检查" % self.option("data_path"))
        self.data_path = self.new_split_path(self.data_path)
        if not os.path.exists(self.data_path):
            self.set_error("下机数据路径:%s 没有找到，请检查" % self.option("data_path"))
        self.bcl_path = self.data_path + "/Data/Intensities/BaseCalls/"
        cmd = self.bcl2fastq_path + ' -i ' + self.bcl_path + ' -o ' + self.output_dir + ' --sample-sheet '\
              + self.option("sample_sheet").prop["path"] + ' --use-bases-mask ' + self.option("bases_mask")\
              + ' --ignore-missing-bcl ' + '-R ' + self.data_path + ' -r 8 -w 8 -p 15 '\
              + '--barcode-mismatches ' + str(self.option("barcode_mismatch"))
        command = self.add_command('bcl2fastq', cmd, script_dir=False, default_return_code=0, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行bcl2fastq完成")
            with open(self.work_dir+"/bcl2fastq.o", "rb") as f:
                lines = f.readlines()
                m = re.match(".* and (.*) warnings.*", lines[-2])
                if m and int(m.group(1)) > 1:
                    m1 = re.match(".*(Processing.*warnings).*", lines[-2])
                    if m1:
                        self.logger.info(lines[-2])
                        self.set_error(m1.group(1))
                        raise Exception(m1.group(1))
        else:
            self.set_error("运行bcl2fastq出错,可能是测序模式有问题，请检查")
            raise Exception("运行bcl2fastq出错,可能是测序模式有问题，请检查")

    def run_md5sum(self):
        """
        生成md5校验码
        """
        fastq_dir = os.path.join(self.output_dir, "Fastq")
        for f1 in os.listdir(fastq_dir):
            fq_dir = os.path.join(fastq_dir, f1)
            md5_info = {}
            for f2 in os.listdir(fq_dir):
                fq_path = os.path.join(fq_dir, f2)
                md5 = os.popen("md5sum {}".format(fq_path)).readlines()[0].split(" ")[0]
                md5_info[f2] = md5
            with open(os.path.join(fq_dir, "md5sum.txt"), "wb") as w:
                for f in md5_info.keys():
                    w.write(md5_info[f] + "  " + f + "\n")

    def run(self):
        super(Bcl2fastqTool, self).run()
        self.bcl2fastq()
        self.run_md5sum()
        self.end()
