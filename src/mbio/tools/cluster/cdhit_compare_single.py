# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class CdhitCompareSingleAgent(Agent):
    """
    cd-hit-est
    version v1.0
    author: zouxuan
    last modified:2018.7.17
    """

    def __init__(self, parent):
        super(CdhitCompareSingleAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入fasta文件
            {"name": "qunum", "type": "int", "default": 0},  # fasta编号
            {"name": "identity", "type": "float", "default": 0.95},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "memory_limit", "type": "int", "default": 10000},  # 内存大小，0为无限制
            {"name": "method", "type": "int", "default": 0},  # 1为全局比对，0为局部比对
            {"name": "direction", "type": "int", "default": 1},  # 1为双向比对，0为单向比对
            {"name": "num_threads", "type": "int", "default": 2},  # cpu数
            {"name": "select", "type": "int", "default": 1},  # 1为聚类到最相似的类中，0为聚类到第一个符合阈值的类
            {"name": "compare", "type": "string", "default": ""},  # 比对结果输出路径
            {"name": "pre", "type": "string", "default": "gene.geneset.tmp.fa.div-"},  # 文件前缀
            {"name": "output", "type": "outfile", "format": "sequence.fasta"}, # 输出fasta文件
            {"name": "ana_type", "type": "string", "default": "nucl"}  # 输入分析类型，是对核酸聚类还是随蛋白聚类
        ]
        self.add_option(options)
        self.step.add_steps('cdhitcomparesingle')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 50  # add memory limit error by hao.gao @ 20191030

    def step_start(self):
        self.step.cdhitcomparesingle.start()
        self.step.update()

    def step_end(self):
        self.step.cdhitcomparesingle.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31600201")
        if not self.option("compare").strip():
            raise OptionError("必须设置输出路径compare", code="31600202")
        if self.option("ana_type") == "nucl":
            if not 0.75 <= self.option("identity") <= 1:
                raise OptionError("identity必须在0.75，1之间", code="31600203")
        else:
            if not 0.7 <= self.option("identity") <= 1:
                raise OptionError("identity必须在0.7，1之间", code="31600205")
        if not 0 <= self.option("coverage") <= 1:
            raise OptionError("coverage必须在0,1之间", code="31600204")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = self.option("num_threads")
        self._memory = str(self.option("memory_limit") / 1000) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(CdhitCompareSingleAgent, self).end()


class CdhitCompareSingleTool(Tool):
    def __init__(self, config):
        super(CdhitCompareSingleTool, self).__init__(config)
        self._version = '1.0'
        self.cdhit_est_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit-est'
        self.cdhit_prot_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit'

    def run(self):
        super(CdhitCompareSingleTool, self).run()
        self.single_compare()
        self.set_output()

    def nuc_word_len(self):
        """
        核酸序列，根据identity判断确定种子长度
        """
        word_length = 11
        if self.option("identity") >= 0.95:
            word_length = 11
        elif 0.9 <= self.option("identity") < 0.95:
            word_length = 9
        elif 0.88 <= self.option("identity") < 0.9:
            word_length = 7
        elif 0.85 <= self.option("identity") < 0.88:
            word_length = 6
        elif 0.8 <= self.option("identity") < 0.85:
            word_length = 5
        elif 0.75 <= self.option("identity") < 0.8:
            word_length = 4
        return word_length

    def prot_word_len(self):
        """
        蛋白序列，根据identity值确定种子长度
        """
        word_length = 5
        if self.option("identity") >= 0.7:
            word_length = 5
        return word_length

    def single_compare(self):
        #length = self.word_len()
        compare_dir = os.path.join(self.work_dir, "compare_dir")
        if os.path.exists(compare_dir):
            shutil.rmtree(compare_dir)
        os.mkdir(compare_dir)
        if self.option("ana_type") == "nucl":
            word_length = self.nuc_word_len()
            cmd = '%s -i %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -r %s -g %s -T %s' % (
                self.cdhit_est_path, self.option("query").prop['path'],os.path.join(compare_dir, "o"), self.option("identity"),
                self.option("coverage"), word_length, 1, 0, 0,
                self.option("direction"), self.option("select"), self.option("num_threads"))
        else:
            word_length = self.prot_word_len()
            cmd = '%s -i %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -g %s -T %s' % (
                self.cdhit_prot_path, self.option("query").prop['path'], os.path.join(compare_dir, "o"), self.option("identity"),
                self.option("coverage"), word_length, self.option("method"), self.option("memory_limit"), 0, self.option("select"), self.option("num_threads"))
        self.logger.info(cmd)
        command1 = self.add_command('cmd_1', cmd, ignore_error=True)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            # f=open(out_dir + "/o",'r')
            # line=f.readline()
            # if line.startswith('>'):
            #     self.logger.info("compare single succeed")
            # 修改读文件的方法为，fasta检查 by ghd @ 20180717
            # try:
            #     self.option('output', out_dir + '/o')
            #     self.logger.info("compare single succeed at fist time")
            # except:
            # # else:
            #     command1.rerun()
            #     self.wait(command1)
            #     if command1.return_code == 0:
            #         self.option('output', out_dir + "/o")
            self.logger.info("compare single succeed at second time")
            #     else:
            #         self.set_error("compare single failed", code="31600201")
            #         raise Exception("compare single failed")
        elif command1.return_code == -6:
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by hao.gao @ 20191030
        else:
            self.set_error("compare single failed", code="31600202")
            raise Exception("compare single failed")



    def set_output(self):
        out_dir = self.option("compare") + '/' + self.option("pre") + str(self.option("qunum")) + "-"
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)

        destination_file = out_dir + "/o"
        if os.path.exists(destination_file):
            os.remove(destination_file)
        work_file = os.path.join(self.work_dir, "compare_dir", "o")
        if os.path.exists(work_file):
            os.link(work_file, destination_file)
            self.logger.info("链接文件成功！")
        else:
            self.logger.info("链接文件失败！")
            self.set_error("未能成功的生成结果文件！")
        # self.linkdir(self.option("compare") + '/gene.geneset.tmp.fa.div-' + str(self.option("qunum")) + "-",
        #              "gene.geneset.tmp.fa.div-" + str(self.option("qunum")) + "-")
        newdir = os.path.join(self.output_dir, self.option("pre") + str(self.option("qunum")) + "-")
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = os.path.join(self.option("compare") + '/' + self.option("pre") + str(self.option("qunum")) + "-", 'o')
        newfiles = os.path.join(newdir, 'o')
        if os.path.exists(newfiles):
            os.remove(newfiles)
        os.link(oldfiles,newfiles)
        self.end()

    # def linkdir(self, dirpath, dirname):
    #     """
    #     link一个文件夹下的所有文件到本module的output目录
    #     """
    #     allfiles = os.listdir(dirpath)
    #     newdir = os.path.join(self.output_dir, dirname)
    #     if not os.path.exists(newdir):
    #         os.mkdir(newdir)
    #     oldfiles = [os.path.join(dirpath, i) for i in allfiles]
    #     newfiles = [os.path.join(newdir, i) for i in allfiles]
    #     for newfile in newfiles:
    #         if os.path.exists(newfile):
    #             os.remove(newfile)
    #     for i in range(len(allfiles)):
    #         os.link(oldfiles[i], newfiles[i])
