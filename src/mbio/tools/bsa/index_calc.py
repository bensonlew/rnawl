# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
from collections import defaultdict


class IndexCalcAgent(Agent):
    """
    bsa分析亲本子代深度过滤
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2018.02.22
    """
    def __init__(self, parent):
        super(IndexCalcAgent, self).__init__(parent)
        options = [
            {"name": "pop_vcf", "type": "string"},  # pop.final.vcf文件
            {"name": "wp", "type": "string"},  # 野生型亲本名称
            {"name": "mp", "type": "string"},  # 突变型亲本名称
            {"name": "wb", "type": "string"},  # 野生型混池名称
            {"name": "mb", "type": "string"},  # 突变型混池名称
            {"name": "pdep", "type": "int", "default": 10},  # 亲本测序深度，默认10X
            {"name": "bdep", "type": "int", "default": 10},  # 混池测序深度，默认10X
            {"name": "popt", "type": "string", "default": "F2"},  # 群体类型，默认F2（仅qtlseq.pl使用）
            {"name": "variant_type", "type": "string", "default": "ALL"},  # 只输出变异类型是该类型的结果
        ]
        self.add_option(options)
        self.step.add_steps("index_calc")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.index_calc.start()
        self.step.update()

    def stepfinish(self):
        self.step.index_calc.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('pop_vcf'):
            raise OptionError('必须输入pop.final.vcf', code="31500601")
        if not self.option('mb'):
            raise OptionError('必须输入突变型混池名称', code="31500602")
        if self.option('variant_type') not in ["ALL", "SNP", "INDEL"]:
            raise OptionError('输出的突变类型只能是"ALL,SNP,INDEL"中的一种', code="31500603")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(IndexCalcAgent, self).end()


class IndexCalcTool(Tool):
    def __init__(self, config):
        super(IndexCalcTool, self).__init__(config)
        self._version = "v1.0.1"
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/perl '
        self.mutmap_path = self.config.PACKAGE_DIR + '/bsa/mutmap.pl'
        self.qtlseq_path = self.config.PACKAGE_DIR + '/bsa/qtlseq.pl'
        self.gene_dict = defaultdict(set)

    def run(self):
        """
        运行
        :return:
        """
        super(IndexCalcTool, self).run()
        self.logger.info("wp:{}, mp:{}, wb:{}, mb:{}".format(self.option('wp'), self.option('mp'),
                                                             self.option('wb'), self.option('mb')))
        if self.option('wb'):  # 根据混池个数，调用不同的脚本运行
            self.run_qtlseq()
        else:
            self.run_mutmap()
        # self.get_gene_num()
        self.set_output()
        self.end()

    def run_mutmap(self):
        """
        运行mutmap.pl,计算位点index值
        """
        cmd = self.perl_path + self.mutmap_path + ' -vcf {} -out {} -mb {} -pdep {} -bdep {} -type {} -popt {}'.format(
            self.option('pop_vcf'), self.work_dir + '/index-calc.result', self.option('mb'),
            self.option('pdep'), self.option('bdep'), self.option("variant_type"), self.option('popt'))   # 最少一个突变型混池
        if self.option('wp'):
            cmd += ' -wp ' + self.option('wp')
        if self.option('mp'):
            cmd += ' -mp ' + self.option('mp')
        self.logger.info('运行mutmap.pl,计算位点index值')
        command = self.add_command("mutmap_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("mutmap运行完成")
        else:
            self.set_error("mutmap运行出错!", code="31500601")

    def run_qtlseq(self):
        """
        运行qtlseq.pl,计算位点index值
        """
        cmd = self.perl_path + self.qtlseq_path + ' -vcf {} -out {} -wb {}  -mb {}  -pdep {} -bdep {} -popt {} -type {}'.format(
            self.option('pop_vcf'), self.work_dir + '/index-calc.result', self.option('wb'),
            self.option('mb'), self.option('pdep'), self.option('bdep'), self.option('popt'), self.option("variant_type"))
        if self.option('wp'):
            cmd += ' -wp ' + self.option('wp')
        if self.option('mp'):
            cmd += ' -mp ' + self.option('mp')
        self.logger.info('运行qtlseq.pl,计算位点index值')
        command = self.add_command("qtlseq_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qtlseq运行完成")
        else:
            self.set_error("qtlseq运行出错!", code="31500602")

    # def get_gene_num(self):
    #     """
    #     通过pop.summary统计每条染色体的数目
    #     :return:
    #     """
    #     with open(self.option('pop_summary').prop['path'])as fr:
    #         lines = fr.readlines()
    #         for line in lines[2:]:  # 从第三行还是统计，忽略前两行注释及标题
    #             tmp = line.strip().split("\t")
    #             self.gene_dict[tmp[4]].add(tmp[0])  # 和黄总监确认过，统计转录本的个数tmp[0]，而不是统计基因个数tmp[1]
    #     with open(self.work_dir + "/gene.stat", "w+")as fw1:
    #         fw1.write("#chr\tgene\n")
    #         for key in self.gene_dict.keys():
    #             fw1.write(key + '\t' + str(len(self.gene_dict[key])) + '\n')
    #     gene_total = []
    #     snp_total = []
    #     indel_total = []
    #     with open(self.work_dir + "/index-calc.result.stat")as f, open(self.work_dir + "/index-calc.result.final.stat", "w+")as fw:
    #         lines = f.readlines()
    #         fw.write("#chr\tgene\tsnp\tindel\n")
    #         for line in lines[1:]:
    #             tmp = line.strip().split("\t")
    #             if tmp[0] not in self.gene_dict.keys():
    #                 # pass
    #                 self.set_error("基因统计信息中的染色体数与变异位点统计染色体数{}不匹配，请核实!".format(tmp[0]))
    #             else:
    #                 gene_total.append(len(self.gene_dict[tmp[0]]))
    #                 snp_total.append(int(tmp[1]))
    #                 indel_total.append(int(tmp[2]))
    #                 fw.write(tmp[0] + "\t" + str(len(self.gene_dict[tmp[0]])) + "\t" + tmp[1] + "\t" + tmp[2] + "\n")
    #         fw.write("Total\t" + str(sum(gene_total)) + "\t" + str(sum(snp_total)) + "\t" + str(sum(indel_total)) + "\n")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面  --这里加上output中文件已经存在要删除的判断 modified by HONGDONG
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        self.check_index(self.work_dir + "/index-calc.result.index")
        try:
            os.link(self.work_dir + "/index-calc.result.index", self.output_dir + "/index-calc.result.index")
            os.link(self.work_dir + "/index-calc.result.stat", self.output_dir + "/index-calc.result.final.stat")
            os.link(self.work_dir + "/index-calc.result.variant", self.output_dir + "/index-calc.result.variant")
            self.logger.info("设置index_calc分析结果目录成功")

        except Exception as e:
            self.logger.info("设置index_calc分析结果目录失败{}".format(e))
            self.set_error("设置index_calc分析结果目录失败%s", variables=(e), code="31500607")
            self.set_error("设置index_calc分析结果目录失败%s", variables=(e), code="31500603")

    def check_index(self, file_path):
        """
        用于检查index-calc.result.index中文件是不是空的，为空可能是因为样本名字不正确或者群体类型错了
        :param file_path:
        :return:
        """
        if len(open(file_path, 'rU').readlines()) == 1:
            self.set_error("给予您的分析策略，未获得有效标记，请重新调整！", code="31500604")
            # self.set_error("流程中群体类型选择错误，请重新设定群体类型参数再次计算，默认群体类型为F2！")
