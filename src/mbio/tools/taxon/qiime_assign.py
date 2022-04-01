# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
# update_info:拆解了qiime的分类过程，分离库训练和序列分类，使用了重新包装的脚本，脚本直接放在了贝叶斯分类器rdp_classifier_2.11中文件位置固定.
# 或者在环境变量中添加RDP_JAR_PATH环境变量，指向classifier.jar

"""RDP taxon 物种分类工具"""

from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import subprocess


class QiimeAssignAgent(Agent):
    """
    classifier.jar: 2.11
    version v1.1
    author shenghe
    last_modified:2016.5.24
    """
    def __init__(self, parent=None):
        """
        """
        super(QiimeAssignAgent, self).__init__(parent)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 输入fasta文件
            {'name': 'revcomp', 'type': 'bool', 'default': False},  # 序列是否翻转
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            # {"name": "customer_mode", "type": "bool", "default": False},  # customer 自定义数据库
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'}  # 输出序列的分类信息文件
            ]
        self.add_option(options)
        self.step.add_steps('qiime_assign')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.qiime_assign.start()
        self.step.update()

    def step_end(self):
        self.step.qiime_assign.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("fasta").is_set:
            raise OptionError("必须设置参数fasta", code="34200301")
        if self.option("revcomp") not in [True, False]:
            raise OptionError("必须设置序列是否翻转", code="34200302")
        if self.option('database') == "custom_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件", code="34200303")
        else:
            if self.option("database") not in ['silva123/16s_bacteria', 'silva123/16s_archaea',
                                               'silva123/16s', 'silva123/18s_eukaryota', 'silva123',
                                               'silva119/16s_bacteria', 'silva119/16s_archaea',
                                               'silva119/16s', 'silva119/18s_eukaryota',
                                               'unite8.0/its_fungi','unite7.2/its_fungi','unite7.0/its_fungi',
                                               'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS',
                                               'fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA',
                                               'fgr/amoA_archaea', 'fgr/amoA_bacteria',
                                               'maarjam081/AM', 'Human_HOMD', 'Human_HOMD_v15.2', 'Human_HPB', 'Protist_PR2_v4.5',
                                               'silva128/16s_archaea', 'silva128/16s_bacteria',
                                               'silva128/18s_eukaryota', 'silva128/16s',
                                               'silva132/16s_archaea', 'silva132/16s_bacteria',
                                               'silva132/18s_eukaryota', 'silva132/16s',
                                               'silva138/16s_archaea', 'silva138/16s_bacteria',
                                               'silva138/18s_eukaryota', 'silva138/16s',
                                               'greengenes135/16s', 'greengenes135/16s_archaea',
                                               'greengenes135/16s_bacteria',
                                               'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea',
                                               'nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s',
                                               'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi',
                                               'nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria','nt_v20210917/16s',
                                               'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi',
                                               'fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012',
                                               'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012','fgr/nosZ_202012',
                                               'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012','fgr/nirK_202012',
                                               'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012','fgr/pmoA_202012', 'fgr/mmoX_202012']:
                                        # 王兆月 2016.11.14 增加数据库silva128 2016.11.23增加数据库mrcA 2016.11.28增加数据库greengenes135
                raise OptionError("数据库%s不被支持", variables=(self.option("database")), code="34200304")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["seqs_tax_assignments.txt", "xls", "OTU的分类学信息"]
            ])
        super(QiimeAssignAgent, self).end()

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 16
        if self.option('database') == 'custom_mode':
            fasta_size = self.option('ref_fasta').get_size() / 1024.00 / 1024.00 / 1024.00  # 单位为G
            if fasta_size > 3.5:## fix by qingchen.zhang @20200918 change 1.5 to 3.5 参考silva和rdp
                raise OptionError('提供的库fasta序列过大: %sG，暂不支持', variables=(fasta_size), code="34200305")
            self._memory = str(QiimeAssignAgent.max_memory_func(fasta_size)) + 'G'
        else:
            self._memory = '60G'
        self.logger.info('Memory:{}  CPU:{}'.format(self._memory, self._cpu))

    @staticmethod
    def max_memory_func(memory):
        """根据提供的fasta大小（单位G）来设定需要的内存大小（单位G）"""
        return int(round(70 * memory)) + 60


class QiimeAssignTool(Tool):
    """
    Qiime Taxon Classify tool
    """
    def __init__(self, config):
        super(QiimeAssignTool, self).__init__(config)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.train_taxon = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/taxon/rdp_classifier_2.11/train_taxon_by_RDP.py")
        self.RDP_classifier = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/taxon/rdp_classifier_2.11/RDP_classifier.py")

    def run_prepare(self):
        if self.option('revcomp'):
            self.logger.info("revcomp 输入的fasta文件")
            try:
                cmd = self.config.SOFTWARE_DIR + "/bioinfo/seq/fastx_toolkit_0.0.14/revcomp " + self.option('fasta').prop['path'] + " > seqs.fasta"
                subprocess.check_output(cmd, shell=True)
                self.logger.info("revcomp 输入的fasta文件 完成")
                return True
            except subprocess.CalledProcessError:
                self.logger.info("revcomp 出错")
                return False
        else:
            self.logger.info("链接输入文件到工作目录")
            if os.path.exists(self.work_dir + '/seqs.fasta'):
                os.remove(self.work_dir + '/seqs.fasta')
            os.link(self.option('fasta').prop['path'], self.work_dir + "/seqs.fasta")
            self.logger.info("链接输入文件到工作目录 完成")
            return True

    def run_assign(self):
        if self.option('database') == "custom_mode":
            ref_fas = self.option('ref_fasta').prop['path']
            os.system("dos2unix {}".format(self.option('ref_fasta').prop['path'])) ## add_by qingchen.zhang @20200813
            ref_tax = self.option('ref_taxon').prop['path']
            fasta_size = os.path.getsize(ref_fas) / 1024.00 / 1024.00 / 1024.00
            max_memory = QiimeAssignTool.max_memory_func(fasta_size)
            cmd = '/program/Python/bin/python ' + self.train_taxon + ' -t ' + ref_tax + ' -s ' + ref_fas +\
                  ' -o ' + self.work_dir + '/RDP_trained' + ' -m ' + str(max_memory)
            trainer = self.add_command('train', cmd)
            self.logger.info('开始对自定义分类库文件进行RDP训练')
            trainer.run()
            self.wait(trainer)
            if trainer.return_code == 0:
                self.logger.info('训练程序正确完成')
                prop_file = self.work_dir + '/RDP_trained/Classifier.properties'
            else:
                self.set_error('trainer运行出错', code="34200301")
                return
            if max_memory > 60:
                max_memory = 24  # 当训练的库文件增加，classifier适当增加，经验设定给后续classifier使用
        else:
            if self.option('database') in ['rdp11.5/16s', 'rdp11.5/16s_bacteria']:
                max_memory = 50  # 给classifier使用
            else:
                max_memory = 10  # 给classifier使用
            prop_dir = 'RDP_trained_' + '_'.join(self.option('database').split('/'))
            prop_file = os.path.join(self.config.SOFTWARE_DIR, 'database/taxon_db/train_RDP_taxon',

                                     prop_dir, 'Classifier.properties')
        cmd = '/program/Python/bin/python ' + self.RDP_classifier + ' -p ' + prop_file + ' -q seqs.fasta' + ' -o '\
              + self.work_dir + '/seqs_taxon' + ' -c ' + str(self.option('confidence')) + ' -m ' + str(max_memory)
        classifier = self.add_command('classifiy', cmd)
        self.logger.info('开始进行序列分类')
        classifier.run()
        self.wait(classifier)
        if classifier.return_code == 0:
            self.logger.info('分类程序正确完成')
            filename = self.work_dir + '/seqs_taxon/format_classified_with_confidence.txt'
            linkfile = self.output_dir + '/seqs_tax_assignments.txt'
            if os.path.exists(linkfile):
                os.remove(linkfile)
            os.link(filename, linkfile)
            self.option('taxon_file', linkfile)
        else:
            self.set_error("classifier运行出错", code="34200302")

    def run(self):
        super(QiimeAssignTool, self).run()
        if self.run_prepare():
            if len(os.listdir(self.output_dir)) == 0:
                self.run_assign()
            self.end()
        else:
            self.set_error("run_prepare运行出错!", code="34200303")

    @staticmethod
    def max_memory_func(memory):
        """根据提供的fasta大小（单位G）来设定需要的内存大小（单位G）"""
        #return int(round(60 * memory)) + 10
        return int(round(60 * memory)) + 50   #zouguanqing 20180814
