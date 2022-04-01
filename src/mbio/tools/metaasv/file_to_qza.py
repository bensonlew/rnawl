# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
from mbio.packages.metaasv.common_function import link_dir, link_file


class FileToQzaAgent(Agent):
    """
    qiime2 格式转换
    将fastq、fasta、taxonomy转为qza
    """
    def __init__(self, parent):
        super(FileToQzaAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "string"}, ##输入的fastq文件或者二维表格式的文件
            {"name": "type", "type": "string"}, ## 输入的文件类型
            {"name": "prefix", "type": "string"},  ## 文件的前缀
        ]
        self.add_option(options)
        self.step.add_steps('convert')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.convert.start()
        self.step.update()

    def step_end(self):
        self.step.convert.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('input_file'):
            raise OptionError('必须提供输入的fasta序列')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        input_size = float(os.path.getsize(self.option("input_file"))) /(1024*1024*100)
        if input_size >= 1:
            memory = 10 + int(10 * input_size)
        else:
            memory = 10
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(FileToQzaAgent, self).end()


class FileToQzaTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(FileToQzaTool, self).__init__(config)
        self.qiime_path = "program/Python/bin/python"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/file2qza.sh")
        self.fasta_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/fasta2qza.sh")
        self.taxonomy_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/taxonomy2qza.sh")
        self.miniconda3 = os.path.join(self.config.SOFTWARE_DIR, "program/miniconda3/bin")
        self.set_environ(PATH=self.miniconda3)
        self.biom_path = "program/Python/bin/"

    def get_manifest_table(self):
        """
        根据数据文件目录获取配置文件的路径，fq文件夹整理成manifest文件
        :return:
        """
        self.logger.info("开始生成qiime输入文件")
        qiime_input = os.path.join(self.work_dir, "qiime_dir")
        if os.path.exists(qiime_input):
            shutil.rmtree(qiime_input)
        link_dir(self.option("input_file"), qiime_input)
        list_table = os.path.join(self.work_dir, "list.txt")
        if os.path.exists(os.path.join(qiime_input, "list.txt")):
            link_file(os.path.join(qiime_input, "list.txt"), list_table)
            os.remove(os.path.join(qiime_input, "list.txt"))
        else:
            input_f = open(list_table, 'w')
            for file in os.listdir(qiime_input):
                file_na = file
                sample_na = file.strip(".fq")
                input_f.write("{}\t{}\n".format(file_na, sample_na))
            input_f.close()
        self.asv_table = os.path.join(self.work_dir, "manifest.tsv")
        n = 0
        with open(self.asv_table, 'w') as w, open(list_table, 'r') as f:
            w.write("sample-id\tabsolute-filepath\n")
            for line in f:
                line = line.strip().split("\t")
                file_name = line[0]
                sample_name = line[1]
                file_path = os.path.join(qiime_input, file_name)
                w.write("{}\t{}\n".format(sample_name, file_path))
                n += 1
        if n < 1:
            self.set_error("不存在正确的fastq文件，请检查")

    def file_to_biom(self):
        """
        将file文件转为biom格式的文件
        :return:
        """
        self.logger.info("开始生成biom文件格式")
        self.asv_table = os.path.join(self.work_dir, "input_asv.biom")
        cmd = '{}biom convert -i {} -o {} --table-type "OTU table" --to-hdf5'.format(self.biom_path, self.option("input_file").prop['path'], self.asv_table)
        self.logger.info(cmd)
        command = self.add_command('file2biom', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("file2biom转格式成功！")
        else:
            self.set_error("file2biom转格式失败！")

    def run_qiime2(self):
        """
        输入qiime数据--fastq数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行转格式！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        input_file = os.path.join(self.work_dir, "manifest.tsv")
        self.out_file = os.path.join(self.work_dir, self.option("prefix")+".qza" )
        cmd = '{} {} {} {} {}'.format(self.shell, self.shell_path, qiime2_env, input_file, self.out_file) #1
        self.logger.info(cmd)
        command = self.add_command('file2qza', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2转格式成功！")
        else:
            self.set_error("qiime2转格式失败失败！")

    def run_fasta(self):
        """
        输入fasta数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行转格式！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        input_file = self.option("input_file")
        self.out_file = os.path.join(self.work_dir, self.option("prefix")+".qza")
        cmd = '{} {} {} {} {}'.format(self.shell, self.fasta_path, qiime2_env, input_file, self.out_file) #1
        self.logger.info(cmd)
        command = self.add_command('file2qza', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2转格式成功！")
        else:
            self.set_error("qiime2转格式失败失败！")

    def run_taxonmy(self):
        """
        输入taxonomy
        :return:
        """
        self.logger.info("开始运行qiime2软件进行转格式！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        input_file = self.option("input_file")
        new_input_file = os.path.join(self.work_dir, "taxon_file.xls")
        with open(input_file, 'r') as f, open(new_input_file, "w") as w:
            w.write("Feature ID\tTaxon\n")
            for line in f:
                w.write(line)
        self.out_file = os.path.join(self.work_dir, self.option("prefix")+".qza")
        cmd = '{} {} {} {} {}'.format(self.shell, self.taxonomy_path, qiime2_env, new_input_file, self.out_file) #1
        self.logger.info(cmd)
        command = self.add_command('file2qza', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2转格式成功！")
        else:
            self.set_error("qiime2转格式失败失败！")


    def run(self):
        """
        运行
        :return:
        """
        super(FileToQzaTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            if self.option("type") in ['fq_dir']:
                self.get_manifest_table()
                self.run_qiime2()
            elif self.option("type") in ['fasta']:
                self.run_fasta()
            elif self.option("type") in ['taxonomy']:
                self.run_taxonmy()
            self.set_output()
            self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        prefix_file = os.path.basename(self.out_file)
        out_file = os.path.join(self.output_dir, prefix_file)
        if os.path.exists(out_file):
            os.remove(out_file)
        os.link(self.out_file, out_file)
        self.logger.info("设置结果文件目录成功！")