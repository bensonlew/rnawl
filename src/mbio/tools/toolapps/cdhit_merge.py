# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.files.sequence.fasta import FastaFile
import os
import subprocess
from biocluster.core.exceptions import OptionError
from Bio import SeqIO
import numpy as np


class CdhitMergeAgent(Agent):
    """
    宏转录组 cdhit软件合并
    配置参数和资源
    """
    def __init__(self, parent):
        super(CdhitMergeAgent, self).__init__(parent)
        options = [
            {"name": "compare_dir", "type": "infile", "format": "sequence.cdhit_cluster_dir"},  # 输入cd-hit比对后的文件夹
            {"name": "table", "type": "int", "default": 11},  # 给出transeq参数table，11为bacteria
            {"name": "clstr", "type": "int", "default": 1},  # 是否生成cluster文件，0不生成，1生成
            {"name": "num1", "type": "int", "default": 0},  # 单拼去冗余生成的.o文件个数
            {"name": "num2", "type": "int", "default": 0},  # 混拼去冗余生成的.o文件个数
            {"name": "is_trans", "type": "string", "default": "false"},  # 是否需要翻译成蛋白文件
            {"name": "ana_type", "type": "string", "default": "nucl"},  # 输入分析类型，对蛋白还是对核酸聚类
            {"name": "gene_fa", "type": "infile", "format": "sequence.fasta"},  #输入的基因预测的样本核酸序列
            {"name": "gene_faa", "type": "infile", "format": "sequence.fasta"}  #输入的基因预测的样本蛋白序列
        ]
        self.add_option(options)
        self.step.add_steps('cdhitmerge')
        self._memory_increase_step = 50

    def check_options(self):
        """
        参数二次检查
        """
        if not self.option("compare_dir").is_set:
            raise OptionError("必须设置参数compare_dir", code="34401001")
        if not self.option("clstr") in [0, 1]:
            raise OptionError("是否需要中间的过程文件", code="34401002")
        if self.option("num1") < 1:
            raise OptionError("切分文件的个数必须要大于1", code="34401003")
        if self.option("ana_type") in ['nucl']:
            if not self.option("gene_fa").is_set:
                raise OptionError("必须输入核酸文件", code="34401004")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        if self.option("gene_fa").is_set:
            size = float(os.path.getsize(self.option("gene_fa").prop['path'])) / 1024*1024*10## 10M
        elif self.option("gene_faa").is_set:
            size = float(os.path.getsize(self.option("gene_faa").prop['path'])) / 1024*1024*10## 10M
        else:
            size = 20
        if int(size) < 20:
            total_memory = 20
        else:
            if int(size) > 50:
                total_memory = 100
            else:
                total_memory = 50 + int(size)
        self._memory = '{}G'.format(total_memory)

    def end(self):
        """
        运行结束
        """
        super(CdhitMergeAgent, self).end()


class CdhitMergeTool(Tool):
    """
    tool运行模块
    """
    def __init__(self, config):
        super(CdhitMergeTool, self).__init__(config)
        self._version = '1.0'
        self.perl_path = '/program/perl-5.24.0/bin/perl '
        self.merge_path = self.config.SOFTWARE_DIR + '/bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/clstr_merge.pl'
        self.sort_clstr_path = self.config.SOFTWARE_DIR + '/bioinfo/uniGene/scripts/sort_clstr.pl'
        self.trans_path = 'bioinfo/seq/emboss6.6/bin/transeq'
        if not os.path.isfile(os.path.join(self.config.SOFTWARE_DIR, self.trans_path)):
            self.trans_path = 'bioinfo/seq/EMBOSS-6.6.0/emboss/transeq'
        self.choose_script = self.config.PACKAGE_DIR + '/metagbin/'
        self.cat_sh = self.config.PACKAGE_DIR + '/sequence/scripts/cat_seq.sh'
        if self.option("gene_fa").is_set:
            self.gene_fa = self.option("gene_fa").prop['path']
        if self.option("gene_faa").is_set:
            self.gene_fa = self.option("gene_faa").prop['path']
        self.ana_type = self.option("ana_type")
        if self.ana_type in ['nucl']:
            self.out_path = os.path.join(self.work_dir, 'gene.uniGeneset.fa')
        else:
            self.out_path = os.path.join(self.work_dir, 'gene.uniGeneset.faa')

    def run(self):
        """
        运行
        具体逻辑是
        如果ana_type为nucl，则is_trans为true；
        否则ana_type为prot，则is_trans为false
        """
        super(CdhitMergeTool, self).run()
        self.logger.info("开始运行 tool！")
        if self.option("ana_type") in ['nucl']:
            self.run_merge_file()
            if self.option("is_trans") in ['true']:
                self.run_translate()
            else:
                if self.option("gene_faa").is_set:
                    self.run_pick_id_from_genefaa() # 从输入的蛋白序列中挑选出去冗余后的蛋白序列
                else:
                    self.set_error("必须进行翻译，is_trans参数错误", code="34401001")
            self.run_stat_gene()
        else:
            self.run_merge_file()
            if self.option("gene_fa").is_set:
                self.run_pick_id_from_genefa() ## 从输入的核酸序列挑选出去冗余后的核酸序列
                self.run_stat_gene()
            else:
                self.run_stat_gene_prot() ## 对蛋白序列文件进行统计
        if self.option("clstr") in [1]:
            self.run_clstr()
            self.run_sort_clstr()
        self.set_output()
        self.end()

    def run_merge_file(self):
        """
        合并去冗余后的序列
        """
        self.logger.info("正在对去冗余后的序列进行合并")
        cmd = 'cat '
        if self.option("num1") > 0:
            for i in range(0, self.option("num1")):
                cmd += self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/o '
        cmd += ' > ' + self.out_path
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info("cat success!")
        except subprocess.CalledProcessError:
            self.set_error("cat failed!", code="34401002")

    def run_translate(self):
        """
        对核酸序列进行翻译
        注意：这里如果进行翻译，那么必然是通过核酸进行去冗余
        """
        self.logger.info("开始对核酸文件进行翻译")
        real_fastaa = os.path.join(self.work_dir, 'gene.uniGeneset.faa')
        cmd = '%s -sequence %s -table %s -trim -outseq %s' % (
            self.trans_path, self.out_path, self.option("table"), real_fastaa)
        self.logger.info(cmd)
        command = self.add_command('translate', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("translate success!")
        else:
            self.set_error("translate failed", code="34401003")

    def run_clstr(self):
        """
        对去冗余中间文件进行合并
        """
        self.logger.info("开始对去冗余中间文件进行合并")
        for i in range(0, self.option("num1")):
            cmd2 = self.config.SOFTWARE_DIR + self.perl_path + self.merge_path + ' '
            clstr = self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(i) + '-/o.clstr'
            if os.path.exists(clstr):
                if i < self.option("num1") - 1:
                    for j in range(i + 1, self.option("num1")):
                        clstr += " " + self.option("compare_dir").prop['path'] + '/gene.geneset.tmp.fa.div-' + str(
                            j) + '-/vs.' + str(i) + '.clstr '
                cmd2 += clstr + '>> ' + os.path.join(self.work_dir, 'gene.uniGeneset.bak.clstr')
                try:
                    subprocess.check_output(cmd2, shell=True)
                    self.logger.info("clstr" + str(i) + "succeed")
                except subprocess.CalledProcessError:
                    self.set_error("clstr %s failed", code="34401004")
        self.logger.info("对所有去冗余的中间文件合并完成")

    def run_pick_id_from_genefa(self):
        """
        从原始的输入核酸文件中提取出去冗余后的核酸id，
        并根据提取后的核酸id挑选出核酸序列
        """
        self.logger.info("开始从原始基因序列文件中提取出去冗余后的核酸文件!")
        self.fa_path = os.path.join(self.work_dir, 'gene.uniGeneset.fa')
        if os.path.exists(self.fa_path):
            os.remove(self.fa_path)
        # prot_name_list = os.path.join(self.work_dir, 'protein_name.list')
        # prot_name_list = list(seq_record.id for seq_record in SeqIO.parse(self.out_path, 'fasta'))
        # origin_gene_dict = {}
        # for seq_record in SeqIO.parse(self.gene_fa, 'fasta'):
        #     gene_id = seq_record.id
        #     gene_seq = seq_record.seq
        #     origin_gene_dict[gene_id] = gene_seq
        #
        # # txpt_arr = np.intersect1d(prot_name_list, origin_gene_dict.keys())

        self.logger.info("正在进行提取核酸序列!")
        prot_name_list = os.path.join(self.work_dir, 'protein_name.list')
        if os.path.exists(prot_name_list):
            os.remove(prot_name_list)
        with open(prot_name_list, 'w') as w:
            for seq_record in SeqIO.parse(self.out_path, 'fasta'):
                seq_id = seq_record.id
                genes = seq_id.split("_")
                newgene = "_".join(genes[0:len(genes) - 1])
                w.write('{}\n'.format(newgene))

        cmd = '{} {}choose_seqs.pl -f {} -l {} -o {}'.format(
                self.perl_path, self.choose_script, self.gene_fa, prot_name_list, self.fa_path
            )
        self.logger.info(cmd)
        command5 = self.add_command('chooose', cmd)
        command5.run()
        self.wait(command5)
        if command5.return_code == 0:
            self.logger.info("succeed choose!")
        else:
            self.set_error("fasta failed")

        # with open(self.fa_path, 'w') as m:
        #     for gene in prot_name_list:
        #         newgene = gene.rstrip()[0:-2]
        #         if newgene in origin_gene_dict.keys():
        #             m.write(">{}\n".format(newgene))
        #             m.write("{}\n".format(origin_gene_dict[newgene]))
        #         else:
        #             self.logger.info("gene:{}不存在序列！".format(newgene))
        self.logger.info("提取核酸序列完成!")

    def run_pick_id_from_genefaa(self):
        """
        从原始的输入蛋白文件中提取出去冗余后的核酸id，
        并根据提取后的核酸id挑选出蛋白序列
        """
        self.logger.info("开始从原始蛋白文件中提取出去冗余后的蛋白文件!")
        self.faa_path = os.path.join(self.work_dir, 'gene.uniGeneset.faa')
        if os.path.exists(self.faa_path):
            os.remove(self.faa_path)
        # prot_name_list = list(seq_record.id for seq_record in SeqIO.parse(self.out_path, 'fasta'))
        # origin_gene_dict = {}
        # for seq_record in SeqIO.parse(self.gene_fa, 'fasta'):
        #     gene_id = seq_record.id
        #     gene_seq = seq_record.seq
        #     origin_gene_dict[gene_id] = gene_seq
        #
        # txpt_arr = np.intersect1d(prot_name_list, origin_gene_dict.keys())
        prot_name_list = os.path.join(self.work_dir, 'protein_name.list')
        if os.path.exists(prot_name_list):
            os.remove(prot_name_list)
        with open(prot_name_list, 'w') as w:
            for seq_record in SeqIO.parse(self.out_path, 'fasta'):
                seq_id = seq_record.id
                genes = seq_id.split("_")
                newgene = "_".join(genes[0:len(genes) - 1])
                w.write('{}\n'.format(newgene))

        cmd = '{} {}choose_seqs.pl -f {} -l {} -o {}'.format(
                self.perl_path, self.choose_script, self.gene_fa, prot_name_list, self.faa_path
            )
        self.logger.info(cmd)
        command5 = self.add_command('chooose', cmd)
        command5.run()
        self.wait(command5)
        if command5.return_code == 0:
            self.logger.info("succeed choose!")
        else:
            self.set_error("fasta failed")
        # self.logger.info("正在进行提取蛋白序列!")
        # with open(self.faa_path, 'w') as m:
        #     for gene in txpt_arr:
        #         newgene = gene.rstrip("_1")
        #         if newgene in origin_gene_dict.keys():
        #             m.write(">{}\n".format(gene))
        #             m.write("{}\n".format(origin_gene_dict[gene]))
        #         else:
        #             self.logger.info("gene:{}不存在序列！".format(gene))
        self.logger.info("提取蛋白序列完成!")

    def run_stat_gene(self):
        """
        对核酸序列进行信息统计
        """
        self.logger.info("开始对核酸序列进行信息统计")
        fastafile = FastaFile()
        file = os.path.join(self.work_dir, 'geneCatalog_stat.xls')
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

    def run_stat_gene_prot(self):
        """
        对核酸序列进行信息统计
        """
        self.logger.info("开始对核酸序列进行信息统计")
        fastafile = FastaFile()
        file = os.path.join(self.work_dir, 'geneCatalog_stat.xls')
        fastafile.set_path(self.work_dir + '/gene.uniGeneset.faa')
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

    def run_sort_clstr(self):
        """
        对中间文件进行排序
        """
        self.logger.info("开始对中间过程文件进行排序")
        if os.path.exists(self.work_dir + '/gene.uniGeneset.bak.clstr'):
            cmd = '%s %s %s %s' % (self.perl_path, self.sort_clstr_path, self.work_dir + '/gene.uniGeneset.bak.clstr',
                                        self.work_dir + '/gene.uniGeneset.clstr')
            command = self.add_command('sort_clstr', cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("clstr succeed")
                os.remove(self.work_dir + '/gene.uniGeneset.bak.clstr')
                if os.path.exists(os.path.join(self.output_dir, 'gene.uniGeneset.clstr')):
                    os.remove(os.path.join(self.output_dir, 'gene.uniGeneset.clstr'))
                os.link(os.path.join(self.work_dir, 'gene.uniGeneset.clstr'), os.path.join(self.output_dir, 'gene.uniGeneset.clstr'))
            else:
                self.set_error("clstr failed", code="34401005")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        filename = os.path.join(self.work_dir, 'geneCatalog_stat.xls')
        linkfile = os.path.join(self.output_dir, 'geneCatalog_stat.xls')
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(filename, linkfile)

        fa_path = os.path.join(self.output_dir, "gene.uniGeneset.fa")
        if os.path.exists(fa_path):
            os.remove(fa_path)
            if os.path.exists(os.path.join(self.work_dir, "gene.uniGeneset.fa")):
                os.link(os.path.join(self.work_dir, "gene.uniGeneset.fa"), fa_path)
        else:
            if os.path.exists(os.path.join(self.work_dir, "gene.uniGeneset.fa")):
                os.link(os.path.join(self.work_dir, "gene.uniGeneset.fa"), fa_path)

        faa_path = os.path.join(self.output_dir, "gene.uniGeneset.faa")
        if os.path.exists(faa_path):
            os.remove(faa_path)
            if os.path.exists(os.path.join(self.work_dir, "gene.uniGeneset.faa")):
                os.link(os.path.join(self.work_dir, "gene.uniGeneset.faa"), faa_path)
        else:
            if os.path.exists(os.path.join(self.work_dir, "gene.uniGeneset.faa")):
                os.link(os.path.join(self.work_dir, "gene.uniGeneset.faa"), faa_path)

