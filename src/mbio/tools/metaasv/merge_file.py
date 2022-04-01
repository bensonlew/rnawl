# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file


class MergeFileAgent(Agent):
    """
    功能： 将不同类型的文件进行合并
    """

    def __init__(self, parent):
        super(MergeFileAgent, self).__init__(parent)
        options = [
            {"name": "input_dir", "type": "infile", "format": "paternity_test.data_dir"},  # 输入合并的文件夹
            {"name": "type", "type": "string", "default": "table"}, #合并文件的类型
        ]
        self.add_option(options)
        self.step.add_steps('merge')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 30

    def step_start(self):
        self.step.merge.start()
        self.step.update()

    def step_end(self):
        self.step.merge.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("input_dir").is_set:
            raise OptionError("请传入input_dir路径!")
        if not self.option("type"):
            raise OptionError("必须设置参数文件的类型!")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(MergeFileAgent, self).end()

class MergeFileTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(MergeFileTool, self).__init__(config)
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.merge_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.python = "program/Python/bin/python"
        self.python_script = os.path.join(self.config.PACKAGE_DIR, 'metaasv/merge_table2.py')
        self.seqkit = self.config.SOFTWARE_DIR + "/bioinfo/meta/seqkit/seqkit"

    def cat_seq(self):
        """
        合并序列
        :return:
        """
        file_list=os.listdir(self.option('input_dir').prop['path'])
        cmd = self.sh_path + 'cat_seq.sh'
        for file in file_list:
            cmd += ' ' + self.option('input_dir').prop['path'] + '/' + file
        cmd += ' ' + self.work_dir + '/ASV_reps.fasta'
        self.logger.info('运行cat_seq，将sequence进行合并')
        command = self.add_command("cat_seq", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cat_seq合并完成")
        else:
            self.set_error("cat_seq合并失败！")

    def de_duplicate(self):
        """
        对fasta进行去重
        :return:
        """
        #outfile = os.path.join(self.work_dir, "ASV_reps_duplicate.fasta")
        outfile = os.path.join(self.output_dir, "ASV_reps.fasta")
        if os.path.exists(outfile):
            os.remove(outfile)
        inputfile = os.path.join(self.work_dir, "ASV_reps.fasta")
        cmd = '{} rmdup -s -w 0 --seq-type dna -j 2 {} > {}'.format(self.seqkit, inputfile, outfile)
        self.logger.info(cmd)
        command = self.add_command("merge_table", cmd, ignore_error=False, shell=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("merge_table运行完成！")
        else:
            self.set_error("merge_table运行出错！")

    def sort_fasta(self):
        """
        对fasta进行排序从长到短
        :return:
        """
        outfile = os.path.join(self.output_dir, "ASV_reps.fasta")
        if os.path.exists(outfile):
            os.remove(outfile)
        inputfile = os.path.join(self.work_dir, "ASV_reps_duplicate.fasta")
        with open(inputfile) as f,open(outfile,"w") as t:
            data = f.read()
            all = {}
            for i in data.split(">"):
                if i.strip():
                    all[i] = len(i)
            all_sort = sorted(all.items(),key=lambda x:x[1],reverse=True)
            for x in all_sort:
                t.write(">"+x[0])

    def merge_table(self):
        """
        对asv_table进行合并
        :return:
        """
        """
        运行脚本 rank_abundance.py 计算结果
        :return:
        """
        outfile = os.path.join(self.work_dir, "asv_table")
        if os.path.exists(outfile):
            shutil.rmtree(outfile)
        os.mkdir(outfile)
        cmd = '{} {} {} {}'.format(self.python, self.python_script, self.option("input_dir").prop['path'], outfile)
        self.logger.info(cmd)
        command = self.add_command("merge_table", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("merge_table运行完成！")
        else:
            self.set_error("merge_table运行出错！")

    def merge_stat(self):
        """
        对降噪结果进行统计
        :return:
        """
        input_file = os.path.join(self.work_dir, "DADA2_sequence_info.txt")
        with open(input_file, "w") as w:
            w.write("sample-id\tinput\tfiltered\tpercentage of input passed filter\tdenoised\tnon-chimeric\tpercentage of input non-chimeric\n")
            w.write("#q2:types\tnumeric\tnumeric\tnumeric\tnumeric\tnumeric\tnumeric\n")
            dir_path = self.option("input_dir").prop['path']
            list_dirs = os.listdir(dir_path)
            for file in list_dirs:
                file_path = os.path.join(dir_path, file)
                with open(file_path, 'r') as f:
                    for line in f:
                        line = line.strip().split("\t")
                        if (line[0] == "sample-id") or (line[0] == "#q2:types"):
                            pass
                        else:
                            w.write("\t".join(line) + "\n")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        if self.option("type") in ["table"]:
            if os.path.exists(os.path.join(os.path.join(self.work_dir, "asv_table"))):
                link_dir(os.path.join(os.path.join(self.work_dir, "asv_table")), self.output_dir)
        elif self.option("type") in ["fasta"]:
            if os.path.exists(os.path.join(os.path.join(self.output_dir, "ASV_reps.fasta"))):
                pass
            else:
                self.set_error("合并reads失败，请检查reads合并程序！")
        elif self.option("type") in ["stat"]:
            if os.path.exists(os.path.join(self.work_dir, "DADA2_sequence_info.txt")):
                link_file(os.path.join(self.work_dir, "DADA2_sequence_info.txt"), os.path.join(self.output_dir, "DADA2_sequence_info.txt"))

    def run(self):
        super(MergeFileTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            if self.option("type") in ["table"]:
                self.merge_table()
            elif self.option("type") in ["fasta"]:
                self.cat_seq()
                self.de_duplicate()
                #self.sort_fasta()
            elif self.option("type") in ["stat"]:
                self.merge_stat()
            self.set_output()
            self.end()