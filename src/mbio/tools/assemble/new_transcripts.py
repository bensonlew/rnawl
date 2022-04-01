# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.trans_step import merged_add_code
import os
import shutil
import re


class NewTranscriptsAgent(Agent):
    """
    新转录本预测
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.14
    """
    def __init__(self, parent):
        super(NewTranscriptsAgent, self).__init__(parent)
        options = [
            {"name": "tmap", "type": "infile", "format": "assembly.tmap"},  # compare后的tmap文件
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考序列的注释文件
            {"name": "merged_gtf", "type": "string"},  # 输入文件，拼接后的注释文件
            {"name": "new_trans_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 新转录本注释文件
            {"name": "new_genes_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 新基因gtf文件
            {"name": "new_trans_fa", "type": "outfile", "format": "sequence.fasta"},  # 新转录本注释文件
            {"name": "new_genes_fa", "type": "outfile", "format": "sequence.fasta"},  # 新基因注释文件
            {"name": "add_code_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 增加class_code后的注释文件
            {"name": "change_id_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 转换ID后的注释文件
            {"name": "change_id_fa", "type": "outfile", "format": "sequence.fasta"}  # 转换ID后的序列文件文件
        ]
        self.add_option(options)
        self.step.add_steps("newtranscripts")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.newtranscripts.start()
        self.step.update()

    def stepfinish(self):
        self.step.newtranscripts.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('tmap'):
            raise OptionError('必须输入compare后tmap文件')
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa')
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列的注释文件')
        if not self.option('merged_gtf'):
            raise OptionError('必须输入参考序列merged_gtf')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "3G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["new_transcripts.gtf", "gtf", "新转录本注释文件"],
            ["new_genes.gtf", "gtf", "新基因注释文件"],
            ["new_transcripts.fa", "fa", "新转录本序列文件"],
            ["new_genes.fa", "fa", "新基因序列文件"],
        ])
        super(NewTranscriptsAgent, self).end()


class NewTranscriptsTool(Tool):
    def __init__(self, config):
        super(NewTranscriptsTool, self).__init__(config)
        self._version = "v1.0.1"
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        self.Python_path = 'program/Python/bin/python '
        self.newtranscripts_gtf_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/scripts/assembly_stat.py'
        self.change_id_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/scripts/gtfmerge.pl'
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/"

    def run(self):
        """
        运行
        :return:
        """
        super(NewTranscriptsTool, self).run()
        self.run_change_id()
        # self.run_change_id_fa()
        self.run_newtranscripts_gtf()
        self.run_change_newtranscripts_gtf()
        # self.run_newtranscripts_fa()
        # self.run_gene_fa()
        # self.run_ref_add()
        self.set_output()
        self.end()

    def run_change_id(self):
        """
        运行perl脚本，将class_code为“=”全部替换掉
        """
        merged_gtf = self.work_dir + '/add_code_merged.gtf'
        merged_add_code(self.option('merged_gtf'), self.option('tmap').prop['path'], merged_gtf)
        cmd = self.perl_path + self.change_id_path + " -i %s -compare %s -ref %s -o %schange_id_merged.gtf" % (
            merged_gtf, self.option('tmap').prop['path'], self.option('ref_gtf').prop['path'], self.work_dir + "/")
        self.logger.info('运行perl，转换id')
        command = self.add_command("change_id_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("perl运行完成")
        else:
            self.set_error("perl运行出错!")

    def run_change_id_fa(self):
        """
        运行
        :return:
        """
        cmd = self.gffread_path + "gffread %s -g %s -w change_id_merged.fa" % (
            self.work_dir + "/" + "change_id_merged.gtf", self.option('ref_fa').prop['path'])
        self.logger.info('运行gtf_to_fasta，形成转换ID后的fasta文件')
        command = self.add_command("change_id_merged_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

    def run_newtranscripts_gtf(self):
        """
        运行python，挑出新转录本gtf文件(merged.gtf)
        """
        olddir = os.path.join(self.work_dir, 'merged')
        if not os.path.exists(olddir):
            os.mkdir(olddir)
        cmd = self.Python_path + self.newtranscripts_gtf_path \
            + " -tmapfile %s -transcript_file %s -out_new_trans %snew_transcripts.gtf -out_new_genes %snew_genes.gtf -out_old_trans %sold_trans.gtf -out_old_genes %sold_genes.gtf" % (
                self.option('tmap').prop['path'], self.work_dir + "/add_code_merged.gtf",
                self.work_dir+"/merged/", self.work_dir+"/merged/", self.work_dir+"/merged/", self.work_dir+"/merged/")
        self.logger.info('运行python，挑出新转录本gtf文件')
        command = self.add_command("merged_newtranscripts_gtf_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("python运行完成")
        else:
            self.set_error("python运行出错!")

    def run_change_newtranscripts_gtf(self):
        """
        运行python，挑出新转录本gtf文件(change_id_merged.gtf)
        """
        newdir = os.path.join(self.work_dir, 'changed')
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        cmd = self.Python_path + self.newtranscripts_gtf_path \
            + " -tmapfile %s -transcript_file %s -out_new_trans %snew_transcripts.gtf -out_new_genes %snew_genes.gtf -out_old_trans %sold_trans.gtf -out_old_genes %sold_genes.gtf" % (
                self.option('tmap').prop['path'], self.work_dir + "/change_id_merged.gtf",self.work_dir+"/changed/",
                self.work_dir+"/changed/", self.work_dir+"/changed/", self.work_dir+"/changed/")
        self.logger.info('运行python，挑出新转录本gtf文件')
        command = self.add_command("change_newtranscripts_gtf_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("python运行完成")
        else:
            self.set_error("python运行出错!")

    def run_newtranscripts_fa(self):
        """
        运行python，新转录本fa文件
        """
        cmd = self.gffread_path + "gffread %s -g %s -w new_transcripts.fa" % (
            self.work_dir + "/changed/" + "new_transcripts.gtf", self.option('ref_fa').prop['path'])
        self.logger.info('运行gtf_to_fasta，形成新转录本fasta文件')
        command = self.add_command("newtranscripts_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

    def run_gene_fa(self):
        """
        运行python，新基因fa文件
        """
        cmd = self.gffread_path + "gffread %s -g %s -w new_genes.fa" % (
            self.work_dir + "/changed/" + "new_genes.gtf", self.option('ref_fa').prop['path'])
        self.logger.info('运行gtf_to_fasta，形成新基因fasta文件')
        command = self.add_command("genes_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

    def run_ref_add(self):
        """
        将class_code为u的序列追加到ref.gtf中
        :return:
        """
        ref_file = self.option('ref_gtf').prop['path']
        new_ref_file = self.work_dir + "/" + "changed_ref.gtf"
        shutil.copy2(ref_file, new_ref_file)
        u_gtf_file = self.work_dir + "/" + "new_genes.gtf"
        with open(new_ref_file, 'a')as fw, open(new_ref_file, 'r')as f, open(u_gtf_file, 'r')as fr:
            first_line = f.readline()
            ref_nine_element = first_line.strip().split("\t")[8]
            for line in fr:
                line_split = line.strip().split("\t")
                eight_elements = "\t".join(line_split[0:7])
                # print line_split[8]
                if line_split[2] == 'transcript':
                    m = re.search(r'gene_id\s*\"(\S+)\";\s*transcript_id\s*\"(\S+)\"', line_split[8])
                    if m:
                        gene_id = m.group(1)
                        transcript_id = m.group(2)
                        n1 = re.search(r'gene_id\s*\"(\S+)\";\s*transcript_id\s*\"(\S+)\"', ref_nine_element)
                        n2 = re.search(r'transcript_id\s*\"(\S+)\";\s*gene_id\s*\"(\S+)\"', ref_nine_element)
                        if n1:
                            new_nine_element = 'gene_id \"' + gene_id + '\"; transcript_id \"' + transcript_id + '\";'
                        elif n2:
                            new_nine_element = 'transcript_id \"' + transcript_id + '\"; gene_id \"' + gene_id + '\";'
                        else:
                            print 'the type of ref.gtf error!'
                    else:
                        pass
                else:
                    m = re.search(r'gene_id\s*\"(\S+)\";\s*transcript_id\s*\"(\S+)\";\s*exon_number\s*\"(\S+)\"',
                                  line_split[8])
                    if m:
                        gene_id = m.group(1)
                        transcript_id = m.group(2)
                        exon_number = m.group(3)
                        n1 = re.search(r'gene_id\s*\"(\S+)\";\s*transcript_id\s*\"(\S+)\"', ref_nine_element)
                        n2 = re.search(r'transcript_id\s*\"(\S+)\";\s*gene_id\s*\"(\S+)\"', ref_nine_element)
                        if n1:
                            new_nine_element = 'gene_id \"' + gene_id + '\"; transcript_id \"' + transcript_id + '\"; exon_number \"' + exon_number + '\";'
                        elif n2:
                            new_nine_element = 'transcript_id \"' + transcript_id + '\"; gene_id \"' + gene_id + '\"; exon_number \"' + exon_number + '\";'
                        else:
                            print 'the type of ref.gtf error!'
                new_line = eight_elements + '\t' + new_nine_element + '\n'
                fw.write(new_line)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            if os.path.exists(self.output_dir + "/new_transcripts.gtf"):
                os.remove(self.output_dir + "/new_transcripts.gtf")
            if os.path.exists(self.output_dir + "/new_genes.gtf"):
                os.remove(self.output_dir + "/new_genes.gtf")
            if os.path.exists(self.output_dir + "/old_trans.gtf"):
                os.remove(self.output_dir + "/old_trans.gtf")
            if os.path.exists(self.output_dir + "/old_genes.gtf"):
                os.remove(self.output_dir + "/old_genes.gtf")
            if os.path.exists(self.output_dir + "/change_id_merged.gtf"):
                os.remove(self.output_dir + "/change_id_merged.gtf")
            if os.path.exists(self.output_dir + "/add_code_merged.gtf"):
                os.remove(self.output_dir + "/add_code_merged.gtf")
            os.link(self.work_dir + "/changed/new_transcripts.gtf", self.output_dir + "/new_transcripts.gtf")
            os.link(self.work_dir + "/changed/new_genes.gtf", self.output_dir + "/new_genes.gtf")
            os.link(self.work_dir + "/merged/old_trans.gtf", self.output_dir + "/old_trans.gtf")
            os.link(self.work_dir + "/merged/old_genes.gtf", self.output_dir + "/old_genes.gtf")
            os.link(self.work_dir + "/change_id_merged.gtf", self.output_dir + "/change_id_merged.gtf")
            os.link(self.work_dir + "/add_code_merged.gtf", self.output_dir + "/add_code_merged.gtf")
            self.option('new_trans_gtf').set_path(self.output_dir + "/new_transcripts.gtf")
            self.option('new_genes_gtf').set_path(self.output_dir + "/new_genes.gtf")
            self.option('add_code_gtf').set_path(self.output_dir + "/add_code_merged.gtf")
            self.option('change_id_gtf').set_path(self.output_dir + "/change_id_merged.gtf")
            self.logger.info("设置拼接比较结果目录成功")

        except Exception as e:
            self.logger.info("设置拼接比较分析结果目录失败{}".format(e))
            self.set_error("设置拼接比较分析结果目录失败{}".format(e))
