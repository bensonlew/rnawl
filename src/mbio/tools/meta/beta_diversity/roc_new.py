# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import subprocess
import re
import pandas as pd


class RocNewAgent(Agent):
    """
    调用analysis_lefse.py 脚本进行lefse分析 同时进行两组和随机森林分析
    version v2.0
    author: zhangpeng
    last_modify: 2017.05.09 zhouxuan
    """
    def __init__(self, parent):
        super(RocNewAgent, self).__init__(parent)
        options = [
            {"name": "lefse_input", "type": "infile", "format": "meta.otu.otu_table"},  # 输入文件，biom格式的otu表
            {"name": "lefse_group", "type": "infile", "format": "meta.otu.group_table"},  # 输入分组文件
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "lda_filter", "type": "float", "default": 2.0},
            {"name": "strict", "type": "int", "default": 0},
            {"name": "lefse_gname", "type": "string"},
            {"name": "start_level", "type": "int", "default": 3},
            {"name": "end_level", "type": "int", "default": 7},
            {"name": "method_cal", "type": "string", "default": "student"},
            {"name": "ci", "type": "float", "default": 0.99},
            {"name": "q_test", "type": "string", "default": "fdr"},
            {"name": "tree_number", "type": "int", "default": 500}

        ]
        self.add_option(options)
        self.step.add_steps("run_biom", "tacxon_stat", "plot_lefse")

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("lefse_input").is_set:
            raise OptionError("必须设置输入的otutable文件.", code="32703701")
        if not self.option("lefse_group").is_set:
            raise OptionError("必须提供分组信息文件", code="32703702")
        if not self.option("lefse_gname"):
            raise OptionError("必须设置lefse分组方案名称", code="32703703")
        if self.option("strict") not in [0, 1]:
            raise OptionError("所设严格性超出范围值", code="32703704")
        if self.option("lda_filter") > 4.0 or self.option("lda_filter") < -4.0:
            raise OptionError("所设阈值超出范围值", code="32703705")
        if len(self.option('lefse_gname').split(',')) >= 3:
            raise OptionError("lefse分析不支持大于2个的分组方案", code="32703706")
        for i in self.option('lefse_gname').split(','):
            gnum = self.option('lefse_group').group_num(i)
            if gnum < 2:
                raise OptionError("lefse分析分组类别必须大于2", code="32703707")
        if self.option('start_level') not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            raise OptionError('起始分类水平不在范围内', code="32703708")
        if self.option('end_level') not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            raise OptionError('结束分类水平不在范围内', code="32703709")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '10G'

    def biom_start_callback(self):
        self.step.run_biom.start()
        self.step.update()

    def biom_finish_callback(self):
        self.step.run_biom.finish()
        self.step.update()

    def sum_taxa_start_callback(self):
        self.step.tacxon_stat.start()
        self.step.update()

    def sum_taxa_finish_callback(self):
        self.step.tacxon_stat.finish()
        self.step.update()

    def lefse_start_callback(self):
        self.step.plot_lefse.start()
        self.step.update()

    def lefse_finish_callback(self):
        self.step.plot_lefse.finish()
        self.step.update()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "lefse分析结果输出目录"],
            ["./lefse_LDA.cladogram.png", "png", "lefse分析cladogram结果图片"]
        ])
        super(RocNewAgent, self).end()


class RocNewTool(Tool):
    """
    Lefse tool
    """
    def __init__(self, config):
        super(RocNewTool, self).__init__(config)
        self._version = '1.0.1'
        self.biom_path = "/miniconda2/bin/"
        self.python_path = "/miniconda2/bin/python"
        self.perl_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/perl"
        self.sum_taxa_path = "/miniconda2/bin/"
        self.script_path = "/bioinfo/taxon/scripts/"
        self.plot_lefse_path = self.config.SOFTWARE_DIR + "/bioinfo/statistical/lefse/"
        self._path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self._path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.end_level = 7
        self.start_level = 3
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.cmd_path_lefse_randomforest = self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/RandomForest_two_group.pl'

    def run_biom(self):
        self.add_state("biom_start", data="开始生成biom格式文件")
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR+"/gcc/5.1.0/lib64:$LD_LIBRARY_PATH")
        biom_cmd = self.biom_path + "biom convert -i %s -o otu_taxa_table.biom --table-type \"OTU table\" --process" \
                                    "-obs-metadata taxonomy  --to-hdf5" % self.option('lefse_input').prop["path"]
        self.logger.info("开始运行biom_cmd")
        biom_command = self.add_command("biom_cmd", biom_cmd).run()
        self.wait(biom_command)
        if biom_command.return_code == 0:
            self.logger.info("biom_cmd运行完成")
        else:
            self.set_error("biom_cmd运行出错!", code="32703701")
        self.add_state("biom_finish", data="生成biom格式文件")

    def run_script(self):
        self.add_state("sum_taxa_start", data="开始生成每一水平的物种统计文件")
        if self.option('end_level') >= self.option('start_level'):
            self.start_level = self.option('start_level')
            self.end_level = self.option('end_level')
        else:
            self.end_level = self.option('start_level')
            self.start_level = self.option('end_level')
        if self.end_level == 9:
            level = ','.join([str(i) for i in range(self.start_level, self.end_level)])
        else:
            level = ','.join([str(i) for i in range(self.start_level, self.end_level + 1)])
        self.logger.info(level)
        if self.start_level != 9:
            script_cmd = self.python_path + " %ssummarize_taxa.py -i otu_taxa_table.biom " \
                                            "-o tax_summary_a -L %s -a" % (self.config.SOFTWARE_DIR + self.script_path, level)
            self.logger.info("开始运行script_cmd")
            script_command = self.add_command("script_cmd", script_cmd).run()
            self.wait(script_command)
            if script_command.return_code == 0:
                self.logger.info("script_cmd运行完成")
            else:
                self.set_error("script_cmd运行出错!", code="32703702")
                raise Exception("script_cmd运行出错!")
        self.get_otu_taxon()
        self.remove_parent_otu(self.work_dir + '/tax_summary_a')

    def remove_parent_otu(self, tax_summary_a):
        files = os.listdir(tax_summary_a)
        for i in files:
            if re.search(r'txt$', i):
                _path = os.path.join(tax_summary_a, i)
                if i == 'otu_taxa_table_L9.txt':
                    df = pd.read_table(_path, index_col=0, skiprows=None)
                else:
                    df = pd.read_table(_path, index_col=0, skiprows=1)
                tmp = df.rename(index=lambda x: ';'.join(x.split(';')[self.start_level - 1:]))
                tmp.to_csv(_path, sep='\t')

    def get_otu_taxon(self):
        if self.end_level == 9:
            if not os.path.exists(self.work_dir + '/tax_summary_a'):
                os.mkdir(self.work_dir + '/tax_summary_a')
            otu_taxon_otu = os.path.join(self.work_dir + '/tax_summary_a', "otu_taxa_table_L9.txt")
            with open(self.option('lefse_input').prop['path'], 'r') as r:
                with open(otu_taxon_otu, 'w') as w:
                    line1 = r.next()
                    if re.search(r'Constructed from biom', line1):
                        line1 = r.next()
                    w.write(line1)
                    for line in r:
                        line = re.sub(r'\.0', '', line)
                        line = line.strip('\n').split('\t')
                        name = line[-1].split('; ')
                        name.append(line[0])
                        line[0] = ';'.join(name)
                        line = '\t'.join(line[0:-1]) + '\n'
                        w.write(line)
        else:
            pass

    def format_input(self):
        self.add_state("lefse_start", data="开始进行lefse分析")
        glist = self.option('lefse_gname').split(',')
        self.option('lefse_group').sub_group('./lefse_group', glist)
        plot_cmd = self.python_path + ' ' + self.plot_lefse_path + \
                   "lefse-input.py -i tax_summary_a -g ./lefse_group -o lefse_input.txt"
        self.logger.info("开始运行format_input_cmd")
        plot_command = self.add_command("format_input_cmd", plot_cmd).run()
        self.wait(plot_command)
        if plot_command.return_code == 0:
            self.logger.info("format_input_cmd运行完成")
        else:
            self.set_error("format_input_cmd运行出错!", code="32703703")

    def run_format(self):
        if len(self.option('lefse_gname').split(',')) == 1:
            format_cmd = self.python_path + " " + self.plot_lefse_path + 'format_input.py  lefse_input.txt  lefse_format.txt  -f  r -c 1 -u 2 -o 1000000'
        elif len(self.option('lefse_gname').split(',')) == 2:
            format_cmd = self.python_path + " " + self.plot_lefse_path + 'format_input.py  lefse_input.txt  lefse_format.txt  -f  r -c 1 -s 2 -u 3 -o 1000000'
        self.logger.info("开始运行format_cmd")
        format_command = self.add_command("format_cmd", format_cmd).run()
        self.wait(format_command)
        if format_command.return_code == 0:
            self.logger.info("format_cmd运行完成")
        else:
            self.set_error("format_cmd运行出错!", code="32703704")

    def run_lefse(self):
        cmd = self.python_path + ' %srun_lefse.py lefse_format.txt lefse_LDA.xls ' \
              '-l %s -y %s' % (self.plot_lefse_path, self.option("lda_filter"), self.option("strict"))
        self.logger.info("开始运行run_lefse_cmd")
        self.logger.info(cmd)
        command = self.add_command("run_lefse_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_lefse_cmd运行完成")
        else:
            self.logger.info("run_lefse_cmd运行出错")

    def run_lefse_cmd_check(self, command, line):
        if re.search(r"Error\sin\slda\Wdefault", line):
            command.kill()
            self.set_error("该分组方案的分组类别所含样本量小于3，lda分析出错", code="32703705")

    def plot_res(self):
        cmd = self.python_path + ' %splot_res.py lefse_LDA.xls lefse_LDA.png' \
              ' --dpi 300 --format png --width 20' % (self.plot_lefse_path)
        self.logger.info("开始运行plot_res_cmd")
        command = self.add_command("plot_res_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("plot_res_cmd运行完成")
        else:
            self.logger.info("plot_res_cmd运行出错")

    def run_level_otu(self):
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.config.SOFTWARE_DIR + '/bioinfo/meta/scripts/lefse_ana.pl'
        cmd += ' -o %s ' %(self.work_dir)
        self.logger.info('开始运行相关数据')
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 cmd.r 成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 失败')
            self.set_error('无法生成 cmd.r 文件', code="32703706")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --restore --no-save < %s/cmd.r' % (self.work_dir), shell=True)
            self.logger.info('计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('计算失败')
        self.logger.info('运行metagenomeSeq.pl程序进行metagenomeseq计算完成')

            

    def plot_cladogram(self):
        cmd = '%s %splot_cladogram.py lefse_LDA.xls ' \
              'lefse_LDA.cladogram.png' ' --format png' % (self.python_path ,self.plot_lefse_path)
        self.logger.info("开始运行plot_cladogram_cmd")
        command = self.add_command("plot_cladogram_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("plot_cladogram_cmd运行完成")
        else:
            self.logger.info("plot_cladogram_cmd运行出错")
        self.add_state("lefse_finish", data="lefse分析完成")

    def set_lefse_output(self):
        """
        将结果文件链接至output
        """
        os.system('cp %s %s' % (self.work_dir + '/lefse_LDA.xls', self.work_dir + '/lefse_lda_head.xls'))
        os.system('sed -i "1i\\taxon\tmean\tgroup\tlda\tpvalue" %s' % (self.work_dir + '/lefse_lda_head.xls'))
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        os.link(self.work_dir + '/lefse_lda_head.xls', self.output_dir + '/lefse_LDA.xls')
        file_list = os.listdir(self.work_dir + "/tax_summary_a")
        new_file_list = []
        for i in file_list:
            if re.match('otu_taxa_table_L\d\.txt', i):
                new_file_list.append(i)
        con_list = []
        for f in new_file_list:
            file_path = os.path.join(self.work_dir, "tax_summary_a", f)
            con = pd.read_table(file_path, header=0, sep="\t")
            if f == "otu_taxa_table_L9.txt":
                con = con.drop("taxonomy", axis=1)
            else:
                con = con.rename(columns={'#OTU ID': 'OTU ID'})
            con_list.append(con)
        finally_otu = pd.concat(con_list)
        self.logger.info("finally_otu:{}".format(finally_otu))
        sample = pd.read_table(self.option("grouptable").prop["path"], header=0, sep="\t")
        sample_list = [i for i in sample["#sample"]]
        sample_list.append("OTU ID")
        self.logger.info("sample_list:{}".format(sample_list))
        origin_name = [i for i in finally_otu.columns]
        self.logger.info("origin_name:{}".format(origin_name))
        for sam in origin_name:
            if sam not in sample_list:
                finally_otu = finally_otu.drop(sam, axis=1)
        self.logger.info('finally_otu:{}'.format(finally_otu))
        fin_path = os.path.join(self.output_dir, "lefse_otu_table")
        finally_otu.to_csv(fin_path, sep="\t", index=False)

    def formattable(self, tablepath):
        print "sometime"

    def randomforest_two_group(self):
        cmd = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl ' + self.cmd_path_lefse_randomforest 
        cmd += ' -i %s -o %s' % (self.option('otutable').prop['path'], self.work_dir + '/RandomForest')
        if self.option('lefse_group').is_set:
            cmd += ' -g %s -m %s' % (self.option('grouptable').prop['path'],self.option('grouptable').prop['path'])
        cmd += ' -type 2'  # 这个固定给脚本就可以
        cmd += ' -ntree {}'.format(self.option('tree_number'))  # str(self.option('ntree')) 可以写数个个数
        cmd += ' -method_cal {}'.format(self.option('method_cal'))
        # student,welch,wilcox  最后一个wilcox我写的代码是其他的都用这个
        cmd += ' -alt two.side'  # 固定给，后面如果加可以有two.sided greater less
        cmd += ' -ci {}'.format(self.option('ci'))  # 0-1之间
        cmd += ' -q_test {}'.format(self.option('q_test'))
        # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
        self.logger.info(cmd)
        self.logger.info('运行RandomForest_perl.pl程序进行RandomForest计算')
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 cmd.r 文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('无法生成 cmd.r 文件', code="32703707")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin/R --restore --no-save < %s/cmd_randomforest.r' % (self.work_dir + '/RandomForest'), shell=True)
            self.logger.info('RandomForest计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('RandomForest计算失败')
        self.logger.info('运行RandomForest_perl.pl程序进行RandomForest计算完成')

    def set_random_twogroup_output(self):
        try:
            os.link(self.work_dir + "/RandomForest/two_group_table.xls", self.output_dir + "/two_group_table.xls")
            os.link(self.work_dir + "/RandomForest/randomforest_vimp_table.xls", self.output_dir + "/Random_table.xls")
        except Exception as e:
            self.logger.error("随机森林两组比较结果link出错{}".format(e))
            self.set_error("随机森林两组比较结果link出错", code="32703708")
        else:
            self.logger.info("随机森林两组比较结果link完成")

    def run(self):
        super(RocNewTool, self).run()
        self.run_biom()
        self.run_script()
        self.format_input()
        self.run_format()
        self.run_lefse()
        self.run_level_otu()
        self.set_lefse_output()
        self.randomforest_two_group()
        self.set_random_twogroup_output()
        self.end()
