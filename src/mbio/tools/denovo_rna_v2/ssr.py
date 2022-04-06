# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import re
from biocluster.core.exceptions import OptionError
import pandas as pd
import os

class SsrAgent(Agent):
    """
    对无参转录组序列进行SSR的预测，以及引物设计,序列位置分布预测。
    检查tool的详细的错误可以去tool的结果里面一个err的文件查看
    """
    def __init__(self, parent):
        super(SsrAgent, self).__init__(parent)
        """
        unigene_fa这个有写好的文件格式检查，对于页面的具体参数传递过来的形式再讨论，再看怎么写到misa.ini里面,以及这个默认值是怎么设置
        """
        options = [
            {"name": "unigene_fa", "type": "infile", "format":"denovo_rna_v2.trinity_fasta"},  # 这是fa文件的格式检查
            {"name": "bed", "type": "infile", "format": "gene_structure.bed"},
            {"name": "rept_1", "type": "int", "default": 10},
            {"name": "rept_2", "type": "int", "default": 6},
            {"name": "rept_3", "type": "int", "default": 5},
            {"name": "rept_4", "type": "int", "default": 5},
            {"name": "rept_5", "type": "int", "default": 5},
            {"name": "rept_6", "type": "int", "default": 5},
            {"name": "ssr_distance", "type": "int", "default": 100},
        ]
        self.add_option(options)
        self._memory_increase_step = 30
        self.step.add_steps("SSR_predict")
        self.on("start", self.step_start)
        self.on("end", self.step_end)


    def step_start(self):
        self.step.SSR_predict.start()
        self.step.update()

    def step_end(self):
        self.step.SSR_predict.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        if not self.option("unigene_fa"):
            raise OptionError("必须设置输入文件:组装完成的unigene.fa文件", code = "32005701")

        if not self.option("bed"):
            raise OptionError("必须设置输入文件:比对SSR位置的参考序列bed文件", code = "32005702")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "SSR分析结果目录"],
        ])
        result_dir.add_regexp_rules([
            ["./ssr_analysis_details", " ", "SSR分析详情表"],
            ["./ssr_type", " ", "SSR统计表"],
            ["./ssr_repeats", " ", "SSR类型分布统计表"],
        ])
        super(SsrAgent, self).end()


class SsrTool(Tool):
    def __init__(self, config):
        super(SsrTool, self).__init__(config)
        # 变量名不能有.  这个表示从属关系，只能用下划线
        self._version = '1.0'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.3/bin/Rscript"
        self.python_path = "/miniconda2/bin/python"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        # 有些软件的路径不需要加self.config.SOFTWARE_DIR，这个需要自己测试下, 这个 + 后面的路径一定要加/开始，要不然会出错
        #self.primer3_core = "bioinfo/denovo_rna_v2/SSR/misa/primer3_core"
        #self.primer3_core = self.config.SOFTWARE_DIR + "/bioinfo/denovo_rna_v2/SSR/scripts/misa/primer3_core"
        self.p3_in = self.config.SOFTWARE_DIR + "/bioinfo/denovo_rna_v2/SSR/misa/p3_in.pl"
        self.p3_out = "/bioinfo/denovo_rna_v2/SSR/misa/p3_out.pl"
        #self.p3_settings_file = self.config.SOFTWARE_DIR + "/bioinfo/denovo_rna_v2/primer3-2.3.7/primer3_v1_1_4_default_settings.txt"
        # 因为原来的misa.ini文件一直固定在原先的路径里面，如果同时并且misa.pl调用的是一个死路径，现在作为参数传递
        self.misa_pl = self.config.SOFTWARE_DIR + "/bioinfo/denovo_rna_v2/SSR/misa/misa2.pl"
        self.misa_anno_pl = "/bioinfo/denovo_rna_v2/SSR/misa/misa_anno.pl"
        # misa_anno.pl这里面有个对于line[0]的匹配，主要是为了消除文件的空行
        self.primer3_core = Config().SOFTWARE_DIR + "/bioinfo/denovo_rna_v2/SSR/misa/primer3_core"
        self.script_path = 'bioinfo/rna/scripts'
        self.misa_class_pl = "/bioinfo/denovo_rna_v2/SSR/misa/misa_class.pl"
        self.query_path = self.config.PACKAGE_DIR + "/denovo_rna_v2/primer.sh"

    def make_misa_ini(self):
        """
        解析页面传过来的参数，然后写成标准格式的misa.ini文件,并且要确保它和misa.pl(刘老师修改版)在一个文件夹里
        字典解析，可以考虑分别打印出其中的值
        :param kwargs:
        :return:
        """
        with open(self.work_dir + '/misa.ini','w') as ini:
            param1 = '1' + '-' + str(self.option("rept_1"))
            param2 = '2' + '-' + str(self.option("rept_2"))
            param3 = '3' + '-' + str(self.option("rept_3"))
            param4 = '4' + '-' + str(self.option("rept_4"))
            param5 = '5' + '-' + str(self.option("rept_5"))
            param6 = '6' + '-' + str(self.option("rept_6"))
            param7 = str(self.option("ssr_distance"))
            line1_ideography = 'definition(unit_size,min_repeats):'
            line2_ideography = 'interruptions(max_difference_between_2_SSRs):'
            line1 = [line1_ideography, param1, param2, param3, param4, param5, param6]
            line2 = [line2_ideography, param7]
            ini.write(" ".join(tuple(line1)) + '\n')
            ini.write(" ".join(tuple(line2)) + '\n')

    # 一般只有给客户的文件我们最终才会上传到output.dir, 所以我们的中间文件的写入保存最好都在work_dir里面
    def unigene_fa2regular_fa(self):
        self.logger.info("111")
        # self.logger.info(self.option("unigene_fa").prop['path'])
        with open(self.option("unigene_fa").prop['path'], "r") as f1, open(self.work_dir + '/testt.fasta','w') as w:
            self.logger.info("打开unigene.fa文件，处理名字里面的多余描述")
            for line in f1:
                line = line.strip()
                line = re.sub(r'\t.*|\s.*', '', line)
                w.write(line + '\n')
        # 中间文件testt.fasta.misa和输入文件testt.fasta在同一个目录下面，所以我们还是从这里读取
        cmd1 = "{} {} {} {}".format(self.perl_path, self.misa_pl, self.work_dir + '/testt.fasta', self.work_dir + '/misa.ini')
        self.logger.info("运行misa.pl")
        self.logger.info(cmd1)
        cmd1_obj = self.add_command("cmd1", cmd1, ignore_error=True).run()
        self.wait(cmd1_obj)
        if cmd1_obj.return_code == 0:
            self.logger.info("运行misa.pl完成")
        elif cmd1_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd1_obj.rerun()
            self.wait(cmd1_obj)
            if cmd1_obj.return_code == 0:
                self.logger.info("运行misa.pl完成")
            else:
                self.set_error("运行misa.pl出错", code = "32005703")

        cmd2 = "{} {} {}".format(self.perl_path, self.p3_in, self.work_dir + '/testt.fasta.misa')
        self.logger.info("运行p3_in")
        self.logger.info(cmd2)
        cmd2_obj = self.add_command("cmd2", cmd2, ignore_error=True).run()
        self.wait(cmd2_obj)
        if cmd2_obj.return_code == 0:
            self.logger.info("运行p3_in完成")
        elif cmd2_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd2_obj.rerun()
            self.wait(cmd2_obj)
            if cmd2_obj.return_code == 0:
                self.logger.info("运行p3_in完成")
            else:
                self.set_error("运行p3_in出错", code = "32005704")

        # -output直接写在format里面，会说变量没有定义，这个时候需要加一个引号，变为"-output"
        # 这个软件需要修改一下权限，单独先在命令行运行一次，看是否会权限失败，chmod 755 your_software
        # sh脚本不加#!/bin/bash， 可以单独测试跑通，但是在框架里面一定要加上这个解释器行
        # cmd3 = "{} {} {} {}".format(self.script_path + '/primer.sh', self.primer3_core, self.work_dir + '/testt.fasta.p3in',  self.work_dir + '/testt.fasta.p3in.p3out')
        cmd3 = "{} {} {} {} {}".format("/sh", self.query_path, self.primer3_core, self.work_dir + '/testt.fasta.p3in', self.work_dir + '/testt.fasta.p3in.p3out')
        self.logger.info("运行primer3_core")
        self.logger.info(cmd3)
        cmd3_obj = self.add_command("cmd3", cmd3, ignore_error=True)
        cmd3_obj.software_dir = "/bin"
        cmd3_obj.run()
        self.wait(cmd3_obj)
        if cmd3_obj.return_code == 0:
            self.logger.info("运行primer3_core完成")
        elif cmd3_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd3_obj.rerun()
            self.wait(cmd3_obj)
            if cmd3_obj.return_code == 0:
                self.logger.info("运行primer3_core完成")
            else:
                self.set_error("运行primer3_core出错", code = "32005705")

        cmd4 = "{}  {}  {}".format(self.p3_out, self.work_dir + '/testt.fasta.p3in.p3out', self.work_dir + '/testt.fasta.misa')
        self.logger.info("运行p3_out")
        self.logger.info(cmd4)
        cmd4_obj = self.add_command("cmd4", cmd4, ignore_error=True).run()
        self.wait(cmd4_obj)
        if cmd4_obj.return_code == 0:
            self.logger.info("运行p3_out完成")
        elif cmd4_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd4_obj.rerun()
            self.wait(cmd4_obj)
            if cmd4_obj.return_code == 0:
                self.logger.info("运行p3_out完成")
            else:
                self.set_error("运行p3_out出错", code = "32005706")

        cmd5 = "{} -i {} -orf {}".format(self.misa_anno_pl,
                                         self.work_dir + '/testt.fasta.misa', self.option("bed").prop["path"])
        self.logger.info("运行misa_anno_pl")
        self.logger.info(cmd5)
        cmd5_obj = self.add_command("cmd5", cmd5, ignore_error=True).run()
        self.wait(cmd5_obj)
        if cmd5_obj.return_code == 0:
            self.logger.info("运行misa_anno_pl完成")
        elif cmd5_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd5_obj.rerun()
            self.wait(cmd5_obj)
            if cmd5_obj.return_code == 0:
                self.logger.info("运行misa_anno_pl完成")
            else:
                self.set_error("运行misa_anno_pl出错", code = "32005707")

        with open(self.work_dir + '/testt.fasta.statistics', "r") as f2, open(self.work_dir + '/testt.fasta.txt','w') as w2:
            self.logger.info("写入第二个repeats后面的SSR的分类统计信息")
            count = 0
            for line in f2:
                line = line.lstrip()
                if line.startswith("Repeats"):
                    count += 1
                if count == 2:
                    w2.write(line)

        cmd6 = "{} {}".format(self.misa_class_pl, self.work_dir + '/testt.fasta.txt')
        self.logger.info("运行misa_class.pl")
        self.logger.info(cmd6)
        cmd6_obj = self.add_command("cmd6", cmd6, ignore_error=True).run()
        self.wait(cmd6_obj)
        if cmd6_obj.return_code == 0:
            self.logger.info("运行misa_class.pl完成")
        elif cmd6_obj.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            cmd6_obj.rerun()
            self.wait(cmd6_obj)
            if cmd6_obj.return_code == 0:
                self.logger.info("运行misa_class.pl完成")
            else:
                self.set_error("运行misa_class.pl出错", code = "32005708")
        # 把unigene.fa.misa.txt（自己把脚本里这个由xls变为txt）的最后一列（也就是SSR的位置信息）加到testt.fasta.p3in.results的最后一列
        t1 = pd.read_table(self.work_dir + '/testt.fasta.p3in.results', header=0, sep='\t')
        t2 = pd.read_table(self.work_dir + '/unigene.fa.misa.txt', header=0, sep='\t')
        t2_last_col = t2['ssr_position']
        t3 = pd.concat([t1,t2_last_col], axis=1)
        t3.columns = ['Tm' if x.startswith('Tm') else x for x in t3.columns]
        t3.to_csv(self.work_dir + '/tmp.txt', sep='\t', index=False)

        # 对合并表格tmp.txt的SSR Type进行分类统计,并且计算testt.fasta.txt里面对重复数字1-5， 6-10， 10-15，>15的计算分类统计总数
        # 对unigene.fa.misa.txt的unigene个数和SSR个数大于1的基因进行统计
        # # -*- coding: utf-8 -*-
        # import os
        # import re
        #
        # os.chdir("G:/Single_test_ssr66499/Ssr")
        with open(self.work_dir + '/testt.fasta.txt') as f3, open(self.work_dir + '/unigene.fa.misa.txt') as f4, \
                open(self.work_dir + '/ssr_type.txt', "w") as w3, \
                open(self.work_dir + '/ssr_repeats_class.txt', "w") as w4:
            header_line = f3.readline()
            header_list = [int(x) for x in header_line.strip().split()[1:-1]]
            all_targets = list()
            target = [x for x in header_list if x <= 5]
            all_targets.append(target)
            target2 = [x for x in header_list if 6 <= x <= 10]
            all_targets.append(target2)
            target3 = [x for x in header_list if 11 <= x <= 15]
            all_targets.append(target3)
            target4 = [x for x in header_list if x > 15]
            all_targets.append(target4)
            # 去除空列表，这样就可以避免因为临界值导致的错误
            all_targets = [x for x in all_targets if x]

            num_Mono_nucleotide_repeats = 0
            num_Di_nucleotide_repeats = 0
            num_Tri_nucleotide_repeats = 0
            num_Tetra_nucleotide_repeats = 0
            num_Penta_nucleotide_repeats = 0
            num_Hexa_nucleotide_repeats = 0
            w4.write(
                "SSR type" + "\t" + "repaeat number" + '\n' + " " + "\t" + "1-5" + "\t" + "6-10" + "\t" + "11-15" + "\t" + ">15" + '\n')

            for line in f3:
                target_sum_list = list()
                line_list = [int(x) for x in line.strip().replace('-', '0').split()[1:-1]]
                d = dict(zip(header_list, line_list))
                for each in all_targets:
                    target_sum_list.append(sum(d[x] for x in each))
                # 写入SSR的按照1-5,6-10,11-15，>15的统计信息分类
                class1, class2, class3, class4 = target_sum_list[0], target_sum_list[1], target_sum_list[2], target_sum_list[3]
                w4.write(
                    line.split("\t")[0] + "\t" + str(class1) + "\t" + str(class2) + "\t" + str(class3) + "\t" + str(
                        class4) + '\n')

                line = line.strip().replace('-', '0').split('\t')
                # 用0替代-，这样才能方便统计个数，去除换行符
                # print line,不加strip还是会打印出\n
                line1 = [int(x) for x in line[1:]]

                if len(line[0].replace('/', "")) == 2:
                    num_Mono_nucleotide_repeats += int(line1[-1])

                if len(line[0].replace('/', "")) == 4:
                    num_Di_nucleotide_repeats += int(line1[-1])

                if len(line[0].replace('/', "")) == 6:
                    num_Tri_nucleotide_repeats += int(line1[-1])

                if len(line[0].replace('/', "")) == 8:
                    num_Tetra_nucleotide_repeats += int(line1[-1])

                if len(line[0].replace('/', "")) == 10:
                    num_Penta_nucleotide_repeats += int(line1[-1])

                if len(line[0].replace('/', "")) == 12:
                    num_Hexa_nucleotide_repeats += int(line1[-1])

            num_Total_number_of_identified_SSRs = num_Mono_nucleotide_repeats + num_Di_nucleotide_repeats + \
                                                  num_Tri_nucleotide_repeats + num_Tetra_nucleotide_repeats + \
                                                  num_Penta_nucleotide_repeats + num_Hexa_nucleotide_repeats

            Mono_percent = float(num_Mono_nucleotide_repeats)/num_Total_number_of_identified_SSRs
            Di_percent = float(num_Di_nucleotide_repeats)/num_Total_number_of_identified_SSRs
            Tri_percent = float(num_Tri_nucleotide_repeats)/num_Total_number_of_identified_SSRs
            Tetra_percent = float(num_Tetra_nucleotide_repeats)/num_Total_number_of_identified_SSRs
            Penta_percent = float(num_Penta_nucleotide_repeats)/num_Total_number_of_identified_SSRs
            Hexa_percent = float(num_Hexa_nucleotide_repeats)/num_Total_number_of_identified_SSRs

            header4 = f4.readline()
            unigene = []
            num_repeat = 0
            for i in f4:
                tmp = i.split()
                unigene.append(tmp[0])
            for item in set(unigene):
                if unigene.count(item) > 1:
                    num_repeat += 1
            Number_of_unigenes_containing_SSRs = len(set(unigene))
            Number_of_unigenes_containing_more_than_1_SSRs = num_repeat

            # 写入SSR统计表
            w3.write("ssr_information" + "\t" + "number" + "\t" + "percent" + "\t" + "\n" \
                                                           "total_ssrs" + "\t" + str(
                num_Total_number_of_identified_SSRs) + "\t" + "0" + "\t" + "\n" \
                                                       "num_uni_ssrs" + "\t" + str(
                Number_of_unigenes_containing_SSRs) + "\t" + "0" + "\t" + "\n" \
                                                      "num_uni_1_ssr" + "\t" +
                     str(Number_of_unigenes_containing_more_than_1_SSRs) + "\t" + "0" + "\t" + "\n" \
                                                                           "mono_rep" + "\t" + str(
                num_Mono_nucleotide_repeats) + "\t" + str(Mono_percent) + "\t" + "\n" \
                                               "di_rep" + "\t" + str(num_Di_nucleotide_repeats) + "\t" + str(Di_percent) + "\t" + "\n" \
                                                                                                                 "tri_rep" + "\t" + str(
                num_Tri_nucleotide_repeats) + "\t" + str(Tri_percent) + "\t" + "\n" \
                                              "tetra_rep" + "\t" + str(
                num_Tetra_nucleotide_repeats) + "\t" + str(Tetra_percent) + "\t" + "\n" \
                                                "penta_rep" + "\t" + str(
                num_Penta_nucleotide_repeats) + "\t" + str(Penta_percent) + "\t" + "\n" \
                                                "hexa_rep" + "\t" + str(
                num_Hexa_nucleotide_repeats) + "\t" + str(Hexa_percent) + "\t" + "\n")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        ssr_analysis_details = self.work_dir + '/tmp.txt'
        if os.path.exists(self.output_dir + '/ssr_analysis_details'):
            os.remove(self.output_dir + '/ssr_analysis_details')
        os.link(ssr_analysis_details, os.path.join(self.output_dir, 'ssr_analysis_details'))

        ssr_type = self.work_dir + '/ssr_type.txt'
        if os.path.exists(self.output_dir + '/ssr_type'):
            os.remove(self.output_dir + '/ssr_type')
        os.link(ssr_type, os.path.join(self.output_dir, 'ssr_type'))

        ssr_repeats_class = self.work_dir + '/ssr_repeats_class.txt'
        if os.path.exists(self.output_dir + '/ssr_repeats'):
            os.remove(self.output_dir + '/ssr_repeats')
        os.link(ssr_repeats_class, os.path.join(self.output_dir, 'ssr_repeats'))

        self.logger.info("设置结果目录")

    def run(self):
        super(SsrTool, self).run()
        self.make_misa_ini()
        self.unigene_fa2regular_fa()
        self.set_output()
        self.end()
