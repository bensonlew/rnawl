# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
# __lastmodified__ = 'guhaidong' @ 20180522
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class CdhitCompareBetweenAgent(Agent):
    """
    cd-hit-est-2d
    version v1.0
    author: zouxuan
    last modified:2018.5.22
    """

    def __init__(self, parent):
        super(CdhitCompareBetweenAgent, self).__init__(parent)
        options = [
            {"name": "database", "type": "infile", "format": "sequence.fasta"},  # 输入fasta1文件
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入fasta2文件
            {"name": "dbnum", "type": "int", "default": 0},  # fasta1编号
            {"name": "qunum", "type": "int", "default": 1},  # fasta2编号
            {"name": "identity", "type": "float", "default": 0.95},  # 给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "memory_limit", "type": "int", "default": 10000},  # 内存大小，0为无限制
            {"name": "method", "type": "int", "default": 0},  # 1为全局比对，0为局部比对
            {"name": "direction", "type": "int", "default": 1},  # 1为双向比对，0为单向比对
            {"name": "num_threads", "type": "int", "default": 2},  # cpu数
            {"name": "select", "type": "int", "default": 1},  # 1为聚类到最相似的类中，0为聚类到第一个符合阈值的类
            {"name": "compare", "type": "string", "default": ""},  # 比对结果输出路径
            {"name": "output", "type": "outfile", "format": "sequence.fasta"},  # 输出fasta文件
            {"name": "ana_type", "type": "string", "default": "nucl"}  # 输入分析类型，是对核酸聚类还是随蛋白聚类
        ]
        self.add_option(options)
        self.step.add_steps('cdhitcomparebetween')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cdhitcomparebetween.start()
        self.step.update()

    def step_end(self):
        self.step.cdhitcomparebetween.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="31600101")
        if not self.option("database").is_set:
            raise OptionError("必须设置参数database", code="31600102")
        if not self.option("compare").strip():
            raise OptionError("必须设置输出路径compare", code="31600103")
        if self.option("ana_type") == "nucl":
            if not 0.75 <= self.option("identity") <= 1:
                raise OptionError("identity必须在0.75，1之间", code="31600203")
        else:
            if not 0.7 <= self.option("identity") <= 1:
                raise OptionError("identity必须在0.7，1之间", code="31600205")
        if not 0 <= self.option("coverage") <= 1:
            raise OptionError("coverage必须在0,1之间", code="31600105")
        if not self.option("qunum") > self.option("dbnum"):
            raise OptionError("fasta2编号必须大于fasta1编号", code="31600106")

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
        super(CdhitCompareBetweenAgent, self).end()


class CdhitCompareBetweenTool(Tool):
    def __init__(self, config):
        super(CdhitCompareBetweenTool, self).__init__(config)
        self._version = '1.0'
        self.cdhit_est_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit-est-2d'
        self.cdhit_prot_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit-2d'

    def run(self):
        super(CdhitCompareBetweenTool, self).run()
        self.between_compare()
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

    def between_compare(self):
        #length = self.word_len()
        compare_dir = os.path.join(self.work_dir, "compare_dir")
        if os.path.exists(compare_dir):
            shutil.rmtree(compare_dir)
        os.mkdir(compare_dir)
        if self.option("ana_type") == "nucl":
            word_length = self.nuc_word_len()
            cmd = '%s -i %s -i2 %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -r %s -g %s -T %s' % (self.cdhit_est_path, self.option("database").prop['path'], self.option("query").prop['path'],os.path.join(compare_dir, "vs." + str(self.option("dbnum"))), self.option("identity"), self.option("coverage"), word_length,1, self.option("memory_limit"), 0, self.option("direction"), self.option("select"),self.option("num_threads"))
        else:
            word_length = self.prot_word_len()
            cmd = '%s -i %s -i2 %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -g %s -T %s' % (self.cdhit_prot_path, self.option("database").prop['path'], self.option("query").prop['path'],os.path.join(compare_dir, "vs." + str(self.option("dbnum"))), self.option("identity"), self.option("coverage"), word_length,self.option("method"), self.option("memory_limit"), 0, self.option("select"),self.option("num_threads"))
        self.logger.info(cmd)
        command1 = self.add_command('cmd_1', cmd)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            # self.set_error("compare between failed")  # deleted @ 20180125
            # raise Exception("compare between failed")  # deleted @ 20180125
            # f=open(out_dir + "/vs." + str(self.option("dbnum")),'r')  # deleted @ 20180330
            # line=f.readline()
            # if line.startswith('>'):
            #     self.logger.info("compare between succeed")
            # else:
            # 生成乱码后，判断开头为'>'无效，现改为fa文件检查  by GHD @ 20180330
            # self.option('output').set_path(out_dir + "/vs." + str(self.option("dbnum")))
            # self.option('output').get_info()
            # if self.option('output').prop['file_format'] == 'FASTA':
            #     self.logger.info("first time compare between succeed")
            # else:
            # 生成乱码后，fa文件调get_info检查无效，改为option方法 by GHD @ 20180717
            # try:
            #     self.option('output', out_dir + "/vs." + str(self.option("dbnum")))
            # except:
            #     command1.rerun()
            #     self.wait(command1)
            #     if command1.return_code == 0:
            #         # self.option(output, out_dir + "/o")
            self.logger.info("second time compare between succeed")
            #     else:
            #         self.set_error("second time compare between failed", code="31600101")
            #         raise Exception("second time compare between failed")
        else:
            self.set_error("first time compare between failed", code="31600102")
            raise Exception("first time compare between failed")

    def set_output(self):
        out_dir = self.option("compare") + "/gene.geneset.tmp.fa.div-" + str(self.option("qunum")) + "-"
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        destination_file = out_dir + "/vs." + str(self.option("dbnum"))
        if os.path.exists(destination_file):
            os.remove(destination_file)
        work_file = os.path.join(self.work_dir, "compare_dir", "vs." + str(self.option("dbnum")))
        if os.path.exists(work_file):
            os.link(work_file, destination_file)
            self.logger.info("链接文件成功！")
        else:
            self.logger.info("链接文件失败！")
            self.set_error("未能成功的生成结果文件！")
        # self.linkdir(self.option("compare") + "/gene.geneset.tmp.fa.div-" + str(self.option("qunum")) + "-",
        #              "gene.geneset.tmp.fa.div-" + str(self.option("qunum")) + "-")
        newdir = os.path.join(self.output_dir, "gene.geneset.tmp.fa.div-" + str(self.option("qunum")) + "-")
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = os.path.join(self.option("compare") + '/gene.geneset.tmp.fa.div-' + str(self.option("qunum")) + "-",
                                "vs." + str(self.option("dbnum")))
        newfiles = os.path.join(newdir, "vs." + str(self.option("dbnum")))
        if os.path.exists(newfiles):  # 防止重运行时输出结果失败 by ghd @20180522
            os.remove(newfiles)
        os.link(oldfiles, newfiles)
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
