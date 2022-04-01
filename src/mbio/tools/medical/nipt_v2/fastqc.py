## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"
# last_modify:20170511
import zipfile
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class FastqcAgent(Agent):
    """
    产筛的生信shell部分功能
    处理fastq，得到bed2文件
    version v1.0
    author: moli.zhou
    lasted modified by hongdong 20171218
    """

    def __init__(self, parent):
        super(FastqcAgent, self).__init__(parent)
        options = [
            {"name": "sample_id", "type": "string"},
            {"name": "fastq_path", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "bam_file", "type": "infile", "format": "nipt.bed"}
        ]
        self.add_option(options)
        self.step.add_steps("fastqc")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.fastqc.start()
        self.step.update()

    def stepfinish(self):
        self.step.fastqc.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        # if not self.option('query_amino'):
        #     raise OptionError("必须输入氨基酸序列")
        if not self.option('sample_id'):
            raise OptionError("必须输入样本名")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = '13G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [".bed.2", "bed.2", "信息表"],
        ])
        super(FastqcAgent, self).end()


class FastqcTool(Tool):
    """
    蛋白质互作组预测tool
    """

    def __init__(self, config):
        super(FastqcTool, self).__init__(config)
        self._version = '1.0.1'

        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin')
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')

    def run_tf(self):
        os.mkdir(self.work_dir+'/temp')

        gz_fastqc = 'bioinfo/medical/FastQc/fastqc -t 10 -o {} {}_R1.fastq.gz'\
            .format(self.work_dir+'/temp', os.path.join(self.option('fastq_path').prop['path'],
                                                        self.option('sample_id')))
        self.logger.info(gz_fastqc)
        cmd = self.add_command("gz_fastqc", gz_fastqc).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("gz_fastqc质控成功")
        else:
            raise Exception("gz_fastqc质控出错")

        bam_fastqc = 'bioinfo/medical/FastQc/fastqc -t 10 -o {} {}'. \
            format(self.work_dir+'/temp', self.option('bam_file').prop['path'])
        self.logger.info(bam_fastqc)
        cmd = self.add_command("bam_fastqc", bam_fastqc).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("bam_fastqc质控成功")
        else:
            raise Exception("bam_fastqc质控出错")

        # unzip = 'unzip {}.map.valid_fastqc.zip'.format(self.work_dir+'/temp/' + self.option('sample_id'))
        # self.logger.info(unzip)
        # cmd = self.add_command("unzip", unzip).run()
        # self.wait(cmd)
        # if cmd.return_code == 0:
        #   self.logger.info("解压成功")
        # else:
        #   raise Exception("解压出错")

        file_zip = self.work_dir+'/temp/' + self.option('sample_id') + '.map.valid_fastqc.zip'
        file_zip1 = self.work_dir + '/temp/' + self.option('sample_id') + '_R1_fastqc.zip'
        file_dir = self.work_dir
        fz = zipfile.ZipFile(file_zip, 'r')
        for f in fz.namelist():
            fz.extract(f, file_dir)
        fz1 = zipfile.ZipFile(file_zip1, 'r')
        for f in fz1.namelist():
            fz1.extract(f, file_dir)

        file_ = '{}.map.valid_fastqc/fastqc_data.txt'.format(self.work_dir+'/' + self.option('sample_id'))
        file1 = '{}_R1_fastqc/fastqc_data.txt'.format(self.work_dir+'/' + self.option('sample_id'))
        adapt_percent = 0
        percent_data = []
        with open(file_, 'r+') as f, open(file1, "r+") as f1:
            for line1 in f1:
                line1 = line1.strip().split('\t')
                for i in range(len(line1)):
                    if re.match(r'TruSeq', line1[i]):
                        adapt_percent1 = line1[i - 1]
                        percent_data.append(float(adapt_percent1))
            if len(percent_data) >= 1:
                for m in percent_data:
                    adapt_percent += m
            else:
                adapt_percent = 0
            for line in f:
                line = line.strip()
                line = line.split('\t')
                if line[0] == '%GC':
                    with open('{}.gc'.format(self.option('sample_id')), 'w+') as gc:
                        gc.write(line[0]+'\t'+line[1]+'\n')
                        gc.write("adapt_percent" + '\t' + str(adapt_percent) + '\n')
                        gc.close()
                    break

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir + '/temp')
        for f in results:
            if re.search(r'.*html$', f):
                os.link(self.work_dir + '/temp/' + f, self.output_dir + '/' + f)
        results_gc = os.listdir(self.work_dir)
        for fgc in results_gc:
            if re.search(r'.*gc$', fgc):
                os.link(self.work_dir + '/' + fgc, self.output_dir + '/' + fgc)

        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(FastqcTool, self).run()
        self.run_tf()
        self.set_output()
        self.end()
