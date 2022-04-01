# -*- coding: utf-8 -*-
# __author__ = 'Shijin'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from mbio.files.sequence.fastq_dir import FastqDirFile
from biocluster.core.exceptions import FileError
from biocluster.config import Config

class GenListAgent(Agent):
    """
    diamond version: 0.8.35
    version 1.0
    author: shijin
    last_modify: 20170316
    """
    def __init__(self, parent):
        super(GenListAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "fq_type", "type": "string", "default": "PE"},
            {"name": "samplebase_dir", "type": "outfile", "format": "sequence.fastq_dir"},
            {"name": "table_id", "type":"string"}
            ]
        self.add_option(options)
        self.step.add_steps('change_fq_dir')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.change_fq_dir.start()
        self.step.update()

    def step_end(self):
        self.step.change_fq_dir.finish()
        self.step.update()

    def check_options(self):
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = '4G'

    def end(self):
        super(GenListAgent, self).end()


class GenListTool(Tool):
    def __init__(self, config):
        super(GenListTool, self).__init__(config)
        # self.sample_base = Config().SAMPLE_BASE
        self.sample_base = Config().WORK_DIR + "/sample_base"
        if not os.path.exists(self.sample_base):
            os.mkdir(self.sample_base)
        # self.id = self._id
        self.logger.info(self._id)
        self.fastq = list()
        self.samples = list()
        self.map = dict()

    def get_fastq_info(self):
        self.logger.info("进行get_fastq_info")
        dir_path = os.path.join(self.work_dir, "output")  # 直接以output文件夹作为新fastq的存放路径
        if os.path.exists(dir_path):
            os.system("rm -r {}".format(dir_path))
        os.mkdir(dir_path)
        fq_dir = FastqDirFile()
        fq_dir.set_path(self.option("fastq_dir").prop["path"])
        fq_dir.get_full_info(dir_path)
        self.fastqs = fq_dir.prop["unzip_fastqs"]

    def get_pairs(self):
        for fq_full in self.fastqs:
            fq = os.path.basename(os.path.splitext(fq_full)[0])
            if not (fq.endswith("_1") or fq.endswith("_2")):
                raise FileError("PE端测序文件，应以_1.fq或_2.fq结尾")
            else:
                if fq.endswith("_1"):
                    sample = fq[:-2]
                    if sample not in self.samples:
                        self.samples.append(sample)
                        self.map[sample] = {}
                    self.map[sample]["l"] = os.path.basename(fq_full)
                else:
                    sample = fq[:-2]
                    if sample not in self.samples:
                        self.samples.append(sample)
                        self.map[sample] = {}
                    self.map[sample]["r"] = os.path.basename(fq_full)
        return True

    def get_single_sample(self):
        for fq in self.fastqs:
            sample = os.path.splitext(os.path.basename(fq))[0]
            if sample not in self.samples:
                self.samples.append(sample)
            else:
                raise FileError("单端测序文件以不同后缀出现多次，请重新检查输入文件")
            self.map[sample] = fq

    def gen_list(self):
        self.logger.info("进入gen_list")
        list_path = os.path.join(self.output_dir, "list.txt")
        if self.option("fq_type") == "PE":
            with open(list_path, "w") as file:
                for sample in self.samples:
                    try:
                        file.write(self.map[sample]["l"] + "\t" + sample + "\tl\n")
                        file.write(self.map[sample]["r"] + "\t" + sample + "\tr\n")
                    except:
                        self.logger.info(sample)
        else:
            with open(list_path, "w") as file:
                for sample in self.samples:
                    file.write(self.map[sample] + "\t" + sample + "\n")

    def export2database(self):  # 将样本存放于固定的位置
        self.logger.info("开始导入数据库")
        id = self.option("table_id")
        task_dir = os.path.join(self.sample_base, id)  # 使用table_id生成本地文件夹
        if os.path.exists(task_dir):
            os.system("rm -r {}/*".format(task_dir))
        else:
            os.mkdir(task_dir)
        for fq in self.fastqs:
            fq_name = os.path.basename(fq)
            new_fq = os.path.join(task_dir, fq_name)
            os.link(fq, new_fq)
        os.link(self.output_dir + "/list.txt", task_dir + "/list.txt")
        self.option("samplebase_dir").set_path(task_dir)
        return True

    def run(self):
        """
        运行
        :return:
        """
        super(GenListTool, self).run()
        self.get_fastq_info()
        if self.option("fq_type") == "PE":
            self.get_pairs()
        else:
            self.get_single_sample()
        self.gen_list()
        self.export2database()
        self.end()