# -*- coding: utf-8 -*-
# __author__ = 'zouxuan
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.files.sequence.fasta import FastaFile
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil
from Bio import SeqIO

class CdhitMergeAgent(Agent):
    """
    version v1.0
    author: zouxuan
    last modified:2017.8.18
    last modified:20190424 by qingchen.zhang
    """

    def __init__(self, parent):
        super(CdhitMergeAgent, self).__init__(parent)
        options = [
            {"name": "compare_dir", "type": "infile", "format": "sequence.cdhit_cluster_dir"},  # 输入cd-hit比对后的文件夹
            {"name": "faa", "type": "outfile", "format": "sequence.fasta"},  # 非冗余基因集蛋白序列
            {"name": "fa", "type": "outfile", "format": "sequence.fasta"},  # 非冗余基因集核算序列
            {"name": "table", "type": "int", "default": 11},  # 给出transeq参数table，11为bacteria。
            {"name": "clstr", "type": "int", "default": 1},  # 是否生成cluster文件，0不生成，1生成。
            {"name": "num1", "type": "int", "default": 0},  # 单拼去冗余生成的.o文件个数。
            {"name": "num2", "type": "int", "default": 0},  # 混拼去冗余生成的.o文件个数。
            {"name": "is_trans", "type": "bool", "default": False},  # 是否需要翻译成蛋白文件。
            {"name": "gene_fa", "type": "infile", "format": "sequence.fasta"},  # 基因预测的单拼样本核酸序列
            {"name": "gene_fa_mix", "type": "infile", "format": "sequence.fasta"},  # 基因预测的混拼样本核酸序列
        ]
        self.add_option(options)
        self.step.add_steps('cdhitmerge')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 60  # 每次重运行增加内存60G by xieshichang @ 20200511

    def step_start(self):
        self.step.cdhitmerge.start()
        self.step.update()

    def step_end(self):
        self.step.cdhitmerge.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("compare_dir").is_set:
            raise OptionError("必须设置参数compare_dir", code="31600301")
        if not self.option("clstr") in [0, 1]:
            raise OptionError("clstr参数必须为0或1", code="31600302")
        if self.option("num1") + self.option("num2") < 1:
            raise OptionError("生成的.o文件必须大于1", code="31600303")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        #self._memory = str(len(os.listdir(self.option("compare_dir").prop['path'])) / 6 + 1) + 'G'
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(CdhitMergeAgent, self).end()


class CdhitMergeTool(Tool):
    def __init__(self, config):
        super(CdhitMergeTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = '/program/perl-5.24.0/bin/perl '
        self.python = '/program/Python/bin/python'
        self.merge_path = self.config.SOFTWARE_DIR + '/bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/clstr_merge.pl'
        self.sort_clstr_path = self.config.SOFTWARE_DIR + '/bioinfo/uniGene/scripts/sort_clstr.pl'
        self.trans_path = 'bioinfo/seq/emboss6.6/bin/transeq'
        if not os.path.isfile(os.path.join(self.config.SOFTWARE_DIR, self.trans_path)):
            self.trans_path = 'bioinfo/seq/EMBOSS-6.6.0/emboss/transeq'
        self.choose_script = self.config.PACKAGE_DIR + '/metagbin/'
        self.cat_sh = self.config.PACKAGE_DIR + '/sequence/scripts/cat_seq.sh'

    def run(self):
        super(CdhitMergeTool, self).run()
        self.fasta_merge()
        self.link_stat()
        self.end()

    def fasta_merge(self):
        if os.path.isfile(self.work_dir + '/gene.uniGeneset.bak.clstr'):
            os.remove(self.work_dir + '/gene.uniGeneset.bak.clstr')
        else:
            pass
        if 'clstr' not in self.get_option_object().keys() or self.option("clstr") == 1:
            # 为保证工作流不出错而增加的判断 ^^^^^^^^^^^^
            for i in range(0, self.option("num1")):
                cmd2 = self.config.SOFTWARE_DIR + self.perl_path + self.merge_path + ' '
                clstr = self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/o.clstr '
                if i < self.option("num1") - 1:
                    for j in range(i + 1, self.option("num1")):
                        clstr += self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(
                            j) + '-/vs.' + str(i) + '.clstr '
                cmd2 += clstr + '>> ' + self.work_dir + '/gene.uniGeneset.bak.clstr'
                try:
                    subprocess.check_output(cmd2, shell=True)
                    self.logger.info("clstr" + str(i) + "succeed")
                except subprocess.CalledProcessError:
                    self.set_error("clstr %s failed", variables=(i), code="31600301")
                    # if i == 0 :
                    #   cmd2 += self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/o.clstr '
                    # else:
                    #   cmd2 += self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/vs0.clstr '
        # cmd2 += '>> ' + self.work_dir + '/gene.uniGeneset.clstr'
        cmd1 = 'cat '
        if self.option("num1") > 0:
            for i in range(0, self.option("num1")):
                cmd1 += self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/o '
        if self.option("num2") > 0:
            for i in range(0, self.option("num2")):
                cmd1 += self.option("compare_dir").prop['path'] + '/mix.geneset.tmp.fa.div-' + str(i) + '-/o '
        cmd1 += ' > ' + self.work_dir + '/gene.uniGeneset.faa'
        self.logger.info(cmd1)
        try:
            subprocess.check_output(cmd1, shell=True)
            rm_dump_id = "perl -i.rm_dump_bak -lne '$id = $1 and $a{$id} += 1 if />(.*)/;  print if $a{$id} == 1;' "
            rm_dump_id += self.work_dir + '/gene.uniGeneset.faa'
            self.logger.info(rm_dump_id)
            if os.system(rm_dump_id) == 0:
                self.logger.info("基因集去重复id完成")
            else:
                self.set_error("基因组去重复id失败")
            self.successful('faa')
        except subprocess.CalledProcessError:
            self.set_error("fasta failed", code="31600302")
            raise Exception("fasta failed")
        """
        fastafile = FastaFile()
        file = self.work_dir + '/geneCatalog_stat.xls'
        fastafile.set_path(self.work_dir + '/gene.uniGeneset.fa')
        self.logger.info('成功生成fasta文件夹,开始非冗余基因集信息')
        with open(file, "w") as f:
            f.write("Catalog_genes\tCatalog_total_length(bp)\tCatalog_average_length(bp)\n")
            if fastafile.check():
                info_ = list()
                fastafile.get_info()
                info_.append(fastafile.prop["seq_number"])
                info_.append(fastafile.prop["bases"])
                avg = round(float(fastafile.prop["bases"]) / float(fastafile.prop["seq_number"]), 2)
                avg = str(avg)
                info_.append(avg)
                f.write("\t".join(info_) + "\n")
        self.logger.info('非冗余基因集信息统计完毕！')
        """
        #        self.logger.info(cmd2)
        #        try:
        #            subprocess.check_output(cmd2,shell=True)
        #            self.successful('clstr')
        #            self.end()
        #        except subprocess.CalledProcessError:
        #            self.set_error("clstr failed")
        #        self.new_clstr(self.work_dir + '/gene.uniGeneset.bak.clstr',self.work_dir + '/gene.uniGeneset.clstr')
        #        self.wait()
        if self.option("is_trans"):
            real_fasta = os.path.join(self.work_dir, 'gene.uniGeneset.faa')
            real_fastaa = os.path.join(self.work_dir, 'gene.uniGeneset.faa')
            cmd3 = '%s -sequence %s -table %s -trim -outseq %s' % (
                self.trans_path, real_fasta, self.option("table"), real_fastaa)
            self.logger.info(cmd3)
            command3 = self.add_command('cmd_3', cmd3)
            command3.run()
            self.wait(command3)
            if command3.return_code == 0:
                self.successful('faa')
            else:
                self.set_error("fastaa failed", code="31600303")
                raise Exception("fastaa failed")
        else:
            real_fastaa = os.path.join(self.work_dir, 'gene.uniGeneset.faa')
            if self.option('gene_fa_mix').is_set:
                before_fasta = self.work_dir + "/gene.total.uniGeneset.fa"
            else:
                before_fasta = self.option("gene_fa").prop['path']
            after_fasta = os.path.join(self.work_dir, 'gene.uniGeneset.fa')
            prot_name_list = os.path.join(self.work_dir, 'protein_name.list')
            with open(prot_name_list, 'w') as w:
                for seq_record in SeqIO.parse(real_fastaa, 'fasta'):
                    seq_id = seq_record.id
                    genes = seq_id.split("_")
                    newgene = "_".join(genes[0:len(genes) - 1])
                    w.write('{}\n'.format(newgene))
            if self.option('gene_fa_mix').is_set:
                if os.path.exists(before_fasta):
                    os.remove(before_fasta)
                cmd6 = '{} {} {}'.format(self.cat_sh, self.option("gene_fa").prop['path'], before_fasta)
                self.logger.info(cmd6)
                command6 = self.add_command("cat_reads_cmd1", cmd6, shell=True).run()
                self.wait(command6)
                if command6.return_code == 0:
                    self.logger.info("混拼序列单样本结果运行完成")
                else:
                    self.set_error("混拼序列单样本结果运行失败")
                self.logger.info('运行cat_reads，进行合并')
                cmd7 = '{} {} {}'.format(self.cat_sh, self.option("gene_fa_mix").prop['path'], before_fasta)
                self.logger.info(cmd7)
                command7 = self.add_command("cat_reads_cmd2", cmd7, shell=True).run()
                self.wait(command7)
                if command7.return_code == 0:
                    self.logger.info("混拼序列混拼样本结果运行完成")
                else:
                    self.set_error("混拼序列混拼样本结果运行完成")
            self.logger.info("正在进行提取核酸序列")
            #cmd5 = '{} {}choose_seqs.pl -f {} -l {} -o {}'.format(
            #    self.perl_path, self.choose_script, before_fasta, prot_name_list, after_fasta
            #)
            cmd5 = '{} {}seqchoose.py -i {} -l {} -o {}'.format(
                self.python, self.choose_script, before_fasta, prot_name_list, after_fasta
            )
            self.logger.info(cmd5)
            command5 = self.add_command('cmd_5', cmd5)
            command5.run()
            self.wait(command5)
            if command5.return_code == 0:
                self.successful('fa')
            else:
                self.set_error("fasta failed")
                raise Exception("fasta failed")
        fastafile = FastaFile()
        file = self.work_dir + '/geneCatalog_stat.xls'
        fastafile.set_path(self.work_dir + '/gene.uniGeneset.fa')
        self.logger.info('成功生成fasta文件夹,开始非冗余基因集信息')
        with open(file, "w") as f:
            f.write("Catalog_genes\tCatalog_total_length(bp)\tCatalog_average_length(bp)\n")
            if fastafile.check():
                info_ = list()
                fastafile.get_info()
                info_.append(fastafile.prop["seq_number"])
                info_.append(fastafile.prop["bases"])
                avg = round(float(fastafile.prop["bases"]) / float(fastafile.prop["seq_number"]), 2)
                avg = str(avg)
                info_.append(avg)
                f.write("\t".join(info_) + "\n")
        self.logger.info('非冗余基因集信息统计完毕！')

        if 'clstr' not in self.get_option_object().keys() or self.option("clstr") == 1:
            # 为保证现有工作流不出错而加的判断
            cmd4 = '%s %s %s %s' % (self.perl_path, self.sort_clstr_path, self.work_dir + '/gene.uniGeneset.bak.clstr',
                                    self.work_dir + '/gene.uniGeneset.clstr')
            command4 = self.add_command('cmd_4', cmd4)
            command4.run()
            self.wait(command4)
            if command4.return_code == 0:
                self.logger.info("clstr succeed")
                os.remove(self.work_dir + '/gene.uniGeneset.bak.clstr')
            else:
                self.set_error("clstr failed", code="31600304")

    def successful(self, type):
        self.logger.info(type + " succeed")
        name = "gene.uniGeneset." + type
        filename = os.path.join(self.work_dir, name)
        linkfile = os.path.join(self.output_dir, name)
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename, linkfile)
        if type in ["faa", "fa"]:
            self.option(type, linkfile)
        else:
            pass

    def link_stat(self):
        filename = os.path.join(self.work_dir, 'geneCatalog_stat.xls')
        linkfile = os.path.join(self.output_dir, 'geneCatalog_stat.xls')
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename, linkfile)
