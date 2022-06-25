# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import os, re, subprocess, shutil, time
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class MapAgent(Agent):
    """
    先将输入的fq进行整合，然后用bowtie2进行比对，最后计算丰度
    author: guhaidong
    last_modified: 20180521
    """

    def __init__(self, parent):
        super(MapAgent, self).__init__(parent)
        options = [
            # {"name": "ref_index", "type": "string"},  不需要，tool里已经给出
            {"name": "query", "type": "infile", "format": "sequence.fastq"},
            {"name": "query_dir", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "fa_query", "type": "infile", "format": "sequence.fasta"},
            {"name": "method", "type": "string", "default": "global"},
            {"name": "abund", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.step.add_steps("map")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 20

    def stepstart(self):
        self.step.map.start()
        self.step.update()

    def stepfinish(self):
        self.step.map.finish()
        self.step.update()

    def check_option(self):
        if not self.option('query') or not self.option('query_dir') or not self.option("fa_query"):
            raise OptionError('必须输入文件')
        if self.option("method") not in ['local', 'global']:
            raise OptionError('比对参数错误： %s' % self.option("method"))
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = "30G"


class MapTool(Tool):
    def __init__(self, config):
        super(MapTool, self).__init__(config)
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.fq2fa_path = '/bioinfo/metaGenomic/MaxBin-2.2.5/auxiliary/idba-1.1.1/bin/fq2fa '
        self.perl_path = '/miniconda2/bin/perl '
        self.bowtie2_path = '/bioinfo/align/bowtie2-2.2.9/'
        self.abund_script = self.config.PACKAGE_DIR + '/hmdb/scripts/abund.pl '
        self.cgc_index = self.config.SOFTWARE_DIR + "/bioinfo/CGC_bowtie/index"
        if self.option('method') == 'local':
            self.method = "--very-sensitive-local"
        else:
            self.method = ""
        self.cgc_set = self.config.SOFTWARE_DIR + "/bioinfo/CGC_bowtie/CGC.fa"
        self.abund_output = self.output_dir + "/cgc_abund.txt"

    def run_query(self):
        """
        对输入序列的文件进行解压缩和合并
        :return:
        """
        if self.option('query').is_gz:
            self.treat_qz(self.option('query').path)
        else:
            self.treat_fq(self.option('query').path)

    def run_query_dir(self):
        """
        对输入序列的文件夹进行解压缩和合并
        :return:
        """
        # self.option("query_dir").get_full_info(os.path.join(self.work_dir, "tmp_dir"))
        if self.option('query_dir').has_list_file:
            # 默认压缩的文件为原始数据，所以此时没有单端，而肯定有双端序列
            if re.search(r'\.(fastq|fq)\.gz', self.option('query_dir').fastqs[0]):
                self.treat_pe_gz()  # 要读list.txt
            else:
                self.treat_pe_fq()  # 要读list.txt，考虑有没有单端，还是只有双端
        else:
            if re.search(r'\.(fastq|fq)\.gz', self.option('query_dir').fastqs[0]):
                self.treat_gz(self.option('query_dir').fastqs)
            else:
                self.treat_fq(self.option('query_dir').fastqs)

    def treat_gz(self, fastqs):
        mission_path = os.path.join(self.work_dir, "tmp_fastq")
        if os.path.exists(mission_path):
            shutil.rmtree(mission_path)
        os.mkdir(mission_path)
        if isinstance(fastqs, str):
            self.run_gz(fastqs, mission_path + '/tmp.s.fq')
        elif isinstance(fastqs, list):
            fastq = ' '.join(fastqs)
            self.run_gz(fastq, mission_path + '/tmp.s.fq')
        else:
            self.set_error("unknown fastqs type: %s" % type(fastqs))
            raise Exception("unknown fastqs type: %s " % type(fastqs))

    def run_gz(self, fastq, output_path):
        unzip_cmd = "zcat " + fastq + " > " + output_path
        self.logger.info("start unzip")
        self.logger.info("command: %s" % unzip_cmd)
        try:
            subprocess.check_output(unzip_cmd, shell=True)
            self.logger.info("unzip done")
        except subprocess.CalledProcessError:
            self.set_error("unzip error")
            raise Exception("unzip error")

    def treat_fq(self, fastqs):
        mission_path = os.path.join(self.work_dir, "tmp_fastq")
        if os.path.exists(mission_path):
            shutil.rmtree(mission_path)
        os.mkdir(mission_path)
        if isinstance(fastqs, unicode):
            fastqs = str(fastqs)
        if isinstance(fastqs, str):
            # self.run_cat(fastqs, mission_path + "/tmp.s.fq")  # modified by ghd @ 20181130 propertion fastqs is file name whichout dirpath
            self.run_cat(os.path.join(self.option('query_dir').path, fastqs), mission_path + "/tmp.s.fq")
        elif isinstance(fastqs, list):
            fastqs = " ".join(os.path.join(self.option('query_dir').path, i) for i in fastqs)
            self.run_cat(fastqs, mission_path + "/tmp.s.fq")
        else:
            self.set_error("unknown fastqs type: %s" % type(fastqs))
            raise Exception("unknown fastqs type: %s" % type(fastqs))

    def run_cat(self, fastq, output_path):
        cat_cmd = "cat " + fastq + " > " + output_path
        self.logger.info("start cat")
        self.logger.info("command: %s" % cat_cmd)
        try:
            subprocess.check_output(cat_cmd, shell=True)
            self.logger.info("cat done")
        except subprocess.CalledProcessError:
            self.set_error("cat error")
            raise Exception("cat error")

    def treat_pe_gz(self):
        mission_path = os.path.join(self.work_dir, "tmp_fastq")
        if os.path.exists(mission_path):
            shutil.rmtree(mission_path)
        os.mkdir(mission_path)
        sample_dict = self.get_list()
        for sample in sample_dict:
            for d in sample:
                self.run_gz(sample[d], mission_path + '/tmp.' + d + '.fq')

    def treat_pe_fq(self):
        mission_path = os.path.join(self.work_dir, "tmp_fastq")
        if os.path.exists(mission_path):
            shutil.rmtree(mission_path)
        os.mkdir(mission_path)
        sample_dict = self.get_list()
        for sample in sample_dict:
            for d in sample:
                self.run_cat(sample[d], mission_path + '/tmp.' + d + '.fq')

    def get_list(self):
        list_path = os.path.join(self.option('query_dir').prop['path'], "list.txt")
        sample = {}
        with open(list_path, "rb") as file:
            for line in file:
                line = line.strip().split()
                if len(line) == 3:
                    sample_path = os.path.join(self.option('fastq_dir').prop['path'], line[0])
                    if line[1] not in sample:
                        sample[line[1]] = {line[2]: sample_path}
                    elif line[2] in sample[line[1]].keys():
                        sample[line[1]][line[2]] += " " + sample_path
                    else:
                        sample[line[1]][line[2]] = sample_path
                else:
                    raise Exception("list.txt文件格式有错误")
        if len(sample) != 1:
            raise Exception("只能用一个样本，现在样本数为%s" % len(sample))
        return sample

    def convert_fa(self):
        mission_path = os.path.join(self.work_dir, "tmp_fasta")
        if os.path.exists(mission_path):
            shutil.rmtree(mission_path)
        os.mkdir(mission_path)
        resource = os.path.join(self.work_dir, "tmp_fastq")
        for fastq in ["tmp.l.fq", "tmp.r.fq", "tmp.s.fq"]:
            fq_file = os.path.join(resource, fastq)
            fasta = fastq[:-2] + "fa"
            fa_file = os.path.join(mission_path, fasta)
            if os.path.exists(fq_file):
                self.run_convert(fq_file, fa_file)

    def run_convert(self, fastq, fasta):
        self.logger.info("convert fq to fa")
        cmd = self.fq2fa_path + ' --paired %s %s' % (fastq, fasta)
        command = self.add_command("fq2fa", cmd).run()
        try:
            self.wait(command)
        except:
            self.logger.info("check whether fq2fa run too fast or not ?")
            # self.add_state("finish","isfinished")
            command._is_end = true
        if command.return_code == 0:
            self.logger.info("fq2fa success")
        else:
            self.set_error("fq2fa error")

    def run_fa_query(self):
        self.logger.info("run_fa_query")
        file_dir = os.path.join(self.work_dir, "tmp_fasta")
        if os.path.exists(file_dir):
            shutil.rmtree(file_dir)
        os.mkdir(file_dir)
        os.link(self.option("fa_query").path, os.path.join(self.work_dir, "tmp_fasta/tmp.s.fa"))

    def run_map(self):
        fa_list = os.listdir(os.path.join(self.work_dir, "tmp_fasta"))
        if len(fa_list) == 3:
            # 双端加单端比对
            self.run_pair_map()
            self.run_single_map()
        elif len(fa_list) == 2:
            # 只做双端比对
            self.run_pair_map()
        elif len(fa_list) == 1:
            # 只做单端比对
            self.run_single_map()
        else:
            raise Exception("需要检查路径是否有结果%s" % os.path.join(self.work_dir, "tmp_fasta"))

    def run_pair_map(self):
        """
        如果输入了fq_list，做双端比对，注意检查有没有单端的文件
        :return:
        """
        self.logger.info("运行bowtie2比对 pair reads")
        l_fa = os.path.join(self.work_dir, "tmp_fasta/tmp.l.fa")
        r_fa = os.path.join(self.work_dir, "tmp_fasta/tmp.r.fa")
        sam_file = os.path.join(self.work_dir, "pe.sam")
        cmd = "{}bowtie2 -p 10 --no-unal --no-head --no-sq {} -f -x {} -1 {} -2 {} -S {}".format(self.bowtie2_path,
                                                                                                 self.method,
                                                                                                 self.cgc_index, l_fa,
                                                                                                 r_fa, sam_file)
        self.logger.info("pe bowtie2 command: %s" % cmd)
        command = self.add_command("bowtie2_map_pair", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bowtie2_map_pair运行完成")
        else:
            self.set_error("bowtie2_map_pair运行出错")

    def run_single_map(self):
        """
        如果输入了fq文件，直接做单端比对, 多样性的序列是双端序列合在一起比对，实际上依旧是pe的序列
        :return:
        """
        self.logger.info("运行bowtie2比对 single read")
        s_fa = os.path.join(self.work_dir, "tmp_fasta/tmp.s.fa")
        if self.option('method') == 'local':
            sam_file = os.path.join(self.work_dir, 'pe.sam')
        else:
            sam_file = os.path.join(self.work_dir, 'se.sam')
        cmd = "{}bowtie2 -p 10 --no-unal --no-head --no-sq {} -f -x {} -U {} -S {}".format(self.bowtie2_path,
                                                                                           self.method, self.cgc_index,
                                                                                           s_fa, sam_file)
        self.logger.info("se bowtie2 command: %s" % cmd)
        command = self.add_command("bowtie2_map_single", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bowtie2_map_single运行完成")
        else:
            self.set_error("bowtie2_map_pair运行出错")

    def run_abundance(self):
        """
        计算丰度
        :return:
        """
        self.logger.info("运行丰度统计")
        pe_sam_path = os.path.join(self.work_dir, "pe.sam")
        se_sam_path = os.path.join(self.work_dir, "se.sam")
        cmd = "{} {} -set {} -o {} ".format(self.perl_path, self.abund_script, self.cgc_set, self.abund_output)
        if os.path.exists(pe_sam_path):
            cmd += " -pe %s" % pe_sam_path
        if os.path.exists(se_sam_path):
            cmd += " -se %s" % se_sam_path
        self.logger.info("abund command: %s" % cmd)
        command = self.add_command("abund", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("丰度计算完成")
        else:
            self.set_error("丰度计算失败")

    def set_output(self):
        self.logger.info("设置结果目录")
        # if os.path.exists(self.abund_output):
        #     self.logger.info("exists %s" % self.abund_output)
        # else:
        #     self.logger.info("not exists %s" % self.abund_output)
        self.option("abund").set_path(self.abund_output)
        self.logger.info("设置结果目录成功")

    def wait_file(self, file, count=0):
        count += 1
        if count > 10:
            self.set_error("没有文件%s" % file)
            raise Exception("没有文件%s" % file)
        if os.path.exists(file):
            return
        else:
            self.logger.info("没有，计数%s" % count)
            time.sleep(10)
            self.wait_file(file, count)

    def run(self):
        super(MapTool, self).run()
        if self.option('query_dir').is_set:
            self.run_query_dir()
            self.convert_fa()
        elif self.option("query").is_set:
            self.run_query()
            self.convert_fa()
        elif self.option("fa_query").is_set:
            self.run_fa_query()
        self.run_map()
        self.run_abundance()
        self.wait_file(self.abund_output)
        self.set_output()
        self.end()
